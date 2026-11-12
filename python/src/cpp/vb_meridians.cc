// vb_meridians: for every vertex x of a voltage graph (the shell of a defect cluster, frames fixed along a spanning
// forest), the length of the shortest closed walk through x of each nonzero, non-wrapping voltage class (up to sign),
// keeping every class within MARGIN of the shortest.
//
// Exact by meet in the middle: a closed walk of length l through x has a midpoint m with both halves of length
// <= ceil(l/2), so a breadth-first search from x in the covering graph (states = vertex + accumulated voltage) to
// radius ceil(l/2) reaches m on two sheets whose voltages differ by the walk's voltage, and their depths sum to <= l.
// Every pair of states over one vertex is a closed walk through x, so the minimum over pairs is exact once the
// radius reaches ceil((best + MARGIN) / 2).
//
// Termination: the caller marks as inactive every vertex of a component that has no nonzero non-wrapping class (an
// exact test on the component's cycle voltages). In an active component every vertex has such a closed walk, so each
// search finds one at finite depth and then stops; no radius limit is needed. BUDGET > 0 caps the states of one
// search (the vertex is then reported unresolved); 0 means no cap.
//
// Certification (EXTRA > 0): if x and y are k steps apart, d_c(y) <= d_c(x) + 2k for every class c (walk to x, around,
// and back: the two halves cancel). So if the shortest class at x beats every other by a margin M, every vertex within
// k < M/4 steps of x has that class as its unique shortest one. Phase A searches seeds (no two within SPACING steps)
// with margin EXTRA; a class not found is longer than 2 x (the radius reached), which bounds M from below. Each seed
// certifies the ball of radius floor((M - 1) / 4). Phase B searches every uncertified active vertex (margin MARGIN).
//
// Input (binary, stdin): int32 n, m, MARGIN, EXTRA, SPACING, BUDGET; int32 offs[n+1] (CSR); int32 dst[m];
//   int32 volt[m][6] (lattice coordinates x 12, then box image; the reverse edge carries the negative); int8 active[n].
// Output (binary, stdout): int32 k per vertex (-1: unresolved), then every vertex's k records int32 [6 voltage,
//   length] in vertex order; class sign: first nonzero lattice component positive. A certified vertex has one record, with length -1 (its class, not its length); an
//   unresolved vertex (budget reached) has k = -1. Then int32 [seeds, certified, searched in phase B, unresolved,
//   most sheets in one search].
// Build: g++ -O2 -fopenmp -std=c++17 vb_meridians.cc -o vb_meridians

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <vector>

typedef std::array<int32_t, 6> V6;
struct V6Hash {
    size_t operator()(const V6& v) const {
        uint64_t h = 1469598103934665603ULL;
        for (int32_t x : v) { h ^= (uint32_t)x; h *= 1099511628211ULL; }
        return (size_t)h;
    }
};

static int32_t n, m, MARGIN, EXTRA, SPACING, BUDGET;
static std::vector<int32_t> offs, dst;
static std::vector<V6> volt;
static std::vector<int8_t> active, nonzero;      // nonzero: the edge's voltage is not zero

struct Found { V6 v; int32_t len; };

// Open-addressing map uint64 -> int32 (linear probing). It grows with the entries, not with the graph, and clearing
// touches only the slots used.
struct FlatMap {
    std::vector<uint64_t> key;                     // ~0: empty
    std::vector<int32_t> val;
    std::vector<uint32_t> used;
    uint64_t mask = 0;
    FlatMap() { init(1 << 10); }
    void init(size_t cap) { key.assign(cap, ~0ULL); val.resize(cap); mask = cap - 1; used.clear(); }
    void clear() { for (uint32_t i : used) key[i] = ~0ULL; used.clear(); }
    static uint64_t mix(uint64_t k) { k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33; return k; }
    int32_t* find(uint64_t k) {
        for (uint64_t i = mix(k) & mask; key[i] != ~0ULL; i = (i + 1) & mask) if (key[i] == k) return &val[i];
        return nullptr;
    }
    // the value at k, inserting v if k is absent (then fresh = true)
    int32_t& get(uint64_t k, int32_t v, bool& fresh) {
        if (2 * (used.size() + 1) > key.size()) grow();
        uint64_t i = mix(k) & mask;
        for (; key[i] != ~0ULL; i = (i + 1) & mask) if (key[i] == k) { fresh = false; return val[i]; }
        key[i] = k; val[i] = v; used.push_back((uint32_t)i); fresh = true; return val[i];
    }
    void grow() {
        std::vector<uint64_t> ok; std::vector<int32_t> ov; std::vector<uint32_t> ou;
        ok.swap(key); ov.swap(val); ou.swap(used);
        init(ok.size() * 2);
        for (uint32_t i : ou) {
            uint64_t j = mix(ok[i]) & mask;
            while (key[j] != ~0ULL) j = (j + 1) & mask;
            key[j] = ok[i]; val[j] = ov[i]; used.push_back((uint32_t)j);
        }
    }
};

struct State { int32_t v, s, d, prev; };          // vertex, sheet, depth, previous state over the same vertex (-1)

struct Search {
    FlatMap state_of;                             // sheet << 32 | vertex -> index in q
    FlatMap last_at;                              // vertex -> the last state reached over it
    std::vector<V6> sheets;
    std::unordered_map<V6, int32_t, V6Hash> sheet_of;
    std::vector<State> q;
    std::vector<Found> out;
    int32_t most_sheets = 0;

    // returns the radius reached (every class of length <= 2 x radius is found), 1 << 29 if the search was exhausted,
    // or -1 if the budget was reached
    int32_t run(int32_t x, int32_t margin) {
        state_of.clear(); last_at.clear(); sheets.clear(); sheet_of.clear(); q.clear(); out.clear();
        sheets.push_back(V6{0, 0, 0, 0, 0, 0}); sheet_of[sheets[0]] = 0;
        bool fresh;
        state_of.get((uint64_t)(uint32_t)x, 0, fresh); last_at.get((uint64_t)(uint32_t)x, 0, fresh);
        q.push_back({x, 0, 0, -1});
        int32_t best = 1 << 30, stop = 1 << 30, reached = 1 << 29;
        for (size_t head = 0; head < q.size(); head++) {
            int32_t u = q[head].v, s = q[head].s, du = q[head].d;
            if (du >= stop) { reached = stop; break; }            // states come in order of depth
            if (BUDGET > 0 && (int64_t)q.size() > BUDGET) return -1;
            for (int32_t e = offs[u]; e < offs[u + 1]; e++) {
                int32_t w = dst[e];
                int32_t s2 = s;
                if (nonzero[e]) {                                // (most edges: frames fixed along a spanning forest)
                    V6 vw;
                    for (int k = 0; k < 6; k++) vw[k] = sheets[s][k] + volt[e][k];
                    auto it = sheet_of.find(vw);
                    if (it == sheet_of.end()) { s2 = (int32_t)sheets.size(); sheets.push_back(vw); sheet_of[vw] = s2; }
                    else s2 = it->second;
                }
                int32_t me = (int32_t)q.size();
                state_of.get(((uint64_t)s2 << 32) | (uint32_t)w, me, fresh);
                if (!fresh) continue;
                int32_t& la = last_at.get((uint64_t)(uint32_t)w, me, fresh);
                int32_t p = fresh ? -1 : la;
                la = me;
                q.push_back({w, s2, du + 1, p});
                for (int32_t t = p; t >= 0; t = q[t].prev) {     // pair the new state with every other sheet at w
                    V6 d;
                    for (int k = 0; k < 6; k++) d[k] = sheets[s2][k] - sheets[q[t].s][k];
                    if (d[3] || d[4] || d[5]) continue;          // wraps the box
                    int k0 = d[0] ? 0 : (d[1] ? 1 : (d[2] ? 2 : -1));
                    if (k0 < 0) continue;                        // trivial
                    if (d[k0] < 0) for (int k = 0; k < 6; k++) d[k] = -d[k];
                    int32_t len = du + 1 + q[t].d;
                    bool seen = false;
                    for (auto& f : out)
                        if (f.v == d) { seen = true; if (len < f.len) f.len = len; break; }
                    if (!seen) out.push_back(Found{d, len});
                    if (len < best) { best = len; stop = std::min(stop, (best + margin + 1) / 2); }
                }
            }
        }
        most_sheets = std::max(most_sheets, (int32_t)sheets.size());
        if (!out.empty()) {
            int32_t mn = 1 << 30;
            for (auto& f : out) mn = std::min(mn, f.len);
            out.erase(std::remove_if(out.begin(), out.end(), [&](const Found& f) { return f.len > mn + margin; }), out.end());
        }
        return reached;
    }
};

static void ball(int32_t x, int32_t r, std::vector<int32_t>& o, std::vector<char>& mark) {   // vertices within r steps
    o.clear(); o.push_back(x); mark[x] = 1;
    std::vector<std::pair<int32_t, int32_t>> bq{{x, 0}};
    for (size_t h = 0; h < bq.size(); h++) {
        auto [u, d] = bq[h];
        if (d >= r) continue;
        for (int32_t e = offs[u]; e < offs[u + 1]; e++)
            if (!mark[dst[e]]) { mark[dst[e]] = 1; o.push_back(dst[e]); bq.push_back({dst[e], d + 1}); }
    }
    for (int32_t v : o) mark[v] = 0;
}

int main() {
    int32_t hdr[6];
    if (fread(hdr, 4, 6, stdin) != 6) return 1;
    n = hdr[0]; m = hdr[1]; MARGIN = hdr[2]; EXTRA = hdr[3]; SPACING = hdr[4]; BUDGET = hdr[5];
    offs.resize(n + 1); dst.resize(m); volt.resize(m); active.resize(n);
    if (fread(offs.data(), 4, n + 1, stdin) != (size_t)n + 1) return 1;
    if (fread(dst.data(), 4, m, stdin) != (size_t)m) return 1;
    if (fread(volt.data(), 4, (size_t)m * 6, stdin) != (size_t)m * 6) return 1;
    if (fread(active.data(), 1, n, stdin) != (size_t)n) return 1;
    nonzero.assign(m, 0);
    for (int32_t e = 0; e < m; e++) for (int k = 0; k < 6; k++) if (volt[e][k]) nonzero[e] = 1;
    std::vector<std::vector<Found>> res(n);
    std::vector<char> done(n, 0), unresolved(n, 0);
    int32_t nseed = 0, ncert = 0, nb = 0, nun = 0, most = 0;
    std::vector<char> mark(n, 0);
    if (EXTRA > 0) {
        // phase A: seeds among active vertices, no two within SPACING steps (greedy in index order)
        std::vector<int32_t> seeds, tmp;
        std::vector<char> blocked(n, 0);
        for (int32_t x = 0; x < n; x++) {
            if (blocked[x] || !active[x]) continue;
            seeds.push_back(x); ball(x, SPACING - 1, tmp, mark);
            for (int32_t v : tmp) blocked[v] = 1;
        }
        nseed = (int32_t)seeds.size();
        std::vector<int32_t> rad(seeds.size(), -1);
        std::vector<std::vector<Found>> sres(seeds.size());
#pragma omp parallel
        {
            Search S;
#pragma omp for schedule(dynamic, 4)
            for (size_t i = 0; i < seeds.size(); i++) {
                int32_t reached = S.run(seeds[i], EXTRA);
                if (reached < 0 || S.out.empty()) continue;          // budget reached: left to phase B
                sres[i] = S.out;
                std::vector<int32_t> L;
                for (auto& f : S.out) L.push_back(f.len);
                std::sort(L.begin(), L.end());
                int32_t second = L.size() > 1 ? L[1] : (reached >= (1 << 29) ? (1 << 28) : 2 * reached + 1);
                rad[i] = std::max(0, (second - L[0] - 1) / 4);
            }
#pragma omp critical
            most = std::max(most, S.most_sheets);
        }
        for (size_t i = 0; i < seeds.size(); i++) {
            if (rad[i] < 0) continue;
            Found best = *std::min_element(sres[i].begin(), sres[i].end(),
                                           [](const Found& a, const Found& b) { return a.len < b.len; });
            int32_t mn = best.len; int nmin = 0;
            for (auto& f : sres[i]) nmin += f.len == mn;
            if (nmin > 1) continue;                                   // a tie at the seed: searched in phase B
            ball(seeds[i], rad[i], tmp, mark);
            for (int32_t v : tmp) {
                if (done[v]) {
                    if (res[v][0].v != best.v) { fprintf(stderr, "certification conflict at %d\n", v); return 2; }
                    continue;
                }
                done[v] = 1; res[v] = {Found{best.v, v == seeds[i] ? mn : -1}}; ncert++;
            }
        }
    }
    std::vector<int32_t> rest;
    for (int32_t x = 0; x < n; x++) if (!done[x] && active[x]) rest.push_back(x);
    nb = (int32_t)rest.size();
#pragma omp parallel
    {
        Search S;
#pragma omp for schedule(dynamic, 16)
        for (size_t i = 0; i < rest.size(); i++) {
            if (S.run(rest[i], MARGIN) < 0) unresolved[rest[i]] = 1;
            else res[rest[i]] = S.out;
        }
#pragma omp critical
        most = std::max(most, S.most_sheets);
    }
    // counts per vertex first (-1: unresolved), then all records (so the caller parses them as arrays)
    for (int32_t x = 0; x < n; x++) {
        int32_t k = unresolved[x] ? -1 : (int32_t)res[x].size();
        nun += unresolved[x];
        fwrite(&k, 4, 1, stdout);
    }
    for (int32_t x = 0; x < n; x++)
        if (!unresolved[x]) for (auto& f : res[x]) { fwrite(f.v.data(), 4, 6, stdout); fwrite(&f.len, 4, 1, stdout); }
    int32_t st[5] = {nseed, ncert, nb, nun, most};
    fwrite(st, 4, 5, stdout);
    return 0;
}
