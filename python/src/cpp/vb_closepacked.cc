// vb_closepacked: the close-packed route (fcc, hcp, and their stacking faults) in C++, writing the same output as
// vb_engine so that the Python side is a FastConnection (closepacked.CPConnection, which writes the parameters and
// the table below):
//   (a) each cell's readings (sets of faces joining first neighbours) from every isomorphism of its topology with a
//       table entry; agreement: a face joins first neighbours for both cells or for neither, readings that no
//       neighbour supports are removed until nothing changes (arc consistency, whose result does not depend on the
//       order), and cells left with none drop out and agreement restarts without them;
//   (b) a decided cell with twelve first neighbours, every Delaunay edge among them decided by a decided endpoint, is
//       labelled by an orientation-preserving isomorphism of its shell with an ideal shell (template); shells whose
//       combinatorial symmetry exceeds their geometric one (hcp-like) are labelled by consistency with labelled
//       neighbours, growing from the cell with the most labelled neighbours;
//   (c) every Delaunay edge between labelled cells must carry a gluing that is a crystal symmetry (first neighbours:
//       from their common neighbours; opposite corners of an octahedral hole: from its middle square).
// The Delaunay complex is the dual of the Voronoi tessellation (a Voronoi vertex shared by four cells is a
// tetrahedron), read from the same voro++ cells.
//
// Input: argv[1] parameters (text): G, then G point-group matrices (row major); K kinds, each "name ambiguous
//   perfect", 12 ideal vectors, then a 12 x 12 adjacency matrix (0/1); the kinds in the order they are tried.
//   argv[2] table (binary int32): entries E; per entry: nf, then per face: n, n vertex ids; then nsets, per set: m, m
//   face ids. argv[3] sample (binary): int32 N; double L[3]; int32 periodic[3]; double x[N][3].
//   argv[4] output prefix: .cells (kept, kind, region, rotation index), .nbrs, .glue, .faces (as vb_engine), .cp
//   (int8 status per atom: 0 ok, 1 unseen, 2 ambiguous, 3 nocell, 4 disagree, 5 undecided shell, 6 shell, 7 glue).
// Build: g++ -O2 -fopenmp -std=c++17 -I<voro>/src vb_closepacked.cc <voro>/libvoro++.a -o vb_closepacked

#include "voro++.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <functional>
#include <unistd.h>
#include <map>
#include <omp.h>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace voro;
typedef std::array<double, 3> V3;
typedef std::array<double, 9> M3;

enum { OK = 0, UNSEEN, AMBIG, NOCELL, DISAGREE, UNDECIDED_SHELL, SHELL, GLUE };

static int N;
static double Lb[3];
static int per[3];
static std::vector<V3> X;
static std::vector<M3> PG;

static V3 minimage(V3 x) {
    for (int k = 0; k < 3; k++) if (per[k]) x[k] -= Lb[k] * std::round(x[k] / Lb[k]);
    return x;
}
static V3 sub(int a, int b) { return minimage(V3{X[a][0] - X[b][0], X[a][1] - X[b][1], X[a][2] - X[b][2]}); }
static double det3(const V3& a, const V3& b, const V3& c) {
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
}
static V3 mvec(const M3& m, const V3& v) {
    return V3{m[0] * v[0] + m[1] * v[1] + m[2] * v[2], m[3] * v[0] + m[4] * v[1] + m[5] * v[2], m[6] * v[0] + m[7] * v[1] + m[8] * v[2]};
}
static M3 mul(const M3& a, const M3& b) {
    M3 c{};
    for (int r = 0; r < 3; r++) for (int s = 0; s < 3; s++) for (int k = 0; k < 3; k++) c[3 * r + s] += a[3 * r + k] * b[3 * k + s];
    return c;
}

// ------------------------------------------------------------------ cells (faces) of the tessellation
struct CellF {
    std::vector<int> nbr;               // neighbour across face f
    std::vector<int> v, o{0};           // face cycles, flat (local vertex ids)
    bool ok = false;
    int nf() const { return (int)nbr.size(); }
    int fs(int f) const { return o[f + 1] - o[f]; }
    const int* fc(int f) const { return v.data() + o[f]; }
};
static std::vector<CellF> cells;

// ------------------------------------------------------------------ planar maps and their isomorphisms
// darts (u, v) of a cell's faces, numbered face by face: next dart within the face, face of the dart, reverse dart
// (-1 if absent), and the sizes of its left and right faces. Built in time linear in the darts, from a list of the
// darts leaving each vertex. Faces come as one flat list of vertex ids (fv) with offsets (fo).
struct Map {
    int nd = 0;
    std::vector<int> nxt, face, inv, rev;            // per dart
};
static void make_map(const int* fv, const int* fo, int nf, Map& M) {
    int nd = fo[nf], nv = 0;
    for (int k = 0; k < nd; k++) nv = std::max(nv, fv[k] + 1);
    M.nd = nd; M.nxt.resize(nd); M.face.resize(nd); M.inv.resize(nd); M.rev.assign(nd, -1);
    std::vector<int> head(nv + 1, 0), out(nd);       // darts by their first vertex (CSR)
    for (int f = 0; f < nf; f++)
        for (int k = fo[f]; k < fo[f + 1]; k++) {
            int k1 = k + 1 == fo[f + 1] ? fo[f] : k + 1;
            M.nxt[k] = k1; M.face[k] = f; head[fv[k] + 1]++;
        }
    for (int v = 0; v < nv; v++) head[v + 1] += head[v];
    { std::vector<int> fill(head.begin(), head.end() - 1); for (int k = 0; k < nd; k++) out[fill[fv[k]]++] = k; }
    for (int k = 0; k < nd; k++) {                   // reverse of u -> w: the dart leaving w that ends at u
        int u = fv[k], w = fv[M.nxt[k]];
        for (int t = head[w]; t < head[w + 1]; t++) if (fv[M.nxt[out[t]]] == u) { M.rev[k] = out[t]; break; }
    }
    for (int k = 0; k < nd; k++) {
        int L = fo[M.face[k] + 1] - fo[M.face[k]], R = M.rev[k] >= 0 ? fo[M.face[M.rev[k]] + 1] - fo[M.face[M.rev[k]]] : 0;
        M.inv[k] = L * 1000 + R;                     // (left, right) face sizes
    }
}
// key: sorted face sizes and the sorted dart invariants (equal for isomorphic cells and mirror images)
static std::vector<int> map_key(const int* fo, int nf, const Map& M) {
    std::vector<int> k; k.reserve(nf + 1 + M.nd);
    for (int f = 0; f < nf; f++) k.push_back(fo[f + 1] - fo[f]);
    std::sort(k.begin(), k.end());
    k.push_back(-1);
    size_t a = k.size(); k.insert(k.end(), M.inv.begin(), M.inv.end()); std::sort(k.begin() + a, k.end());
    return k;
}
struct KeyHash { size_t operator()(const std::vector<int>& v) const { uint64_t h = 1469598103934665603ULL; for (int x : v) { h ^= (uint32_t)x; h *= 1099511628211ULL; } return (size_t)h; } };
// all orientation-preserving isomorphisms a -> b as face maps (from a dart of a with the rarest invariant, tried
// only against darts of b with the same invariant, which every isomorphism must map it to); an isomorphism is fixed
// by the image of one dart, so this finds each exactly once
static void isos(const Map& A, int nfa, const Map& B, int nfb, std::vector<std::vector<int>>& out) {
    out.clear();
    if (A.nd != B.nd || nfa != nfb) return;
    int D = A.nd;
    std::vector<int> sv(A.inv); std::sort(sv.begin(), sv.end());
    auto cnt = [&](int x) { return (int)(std::upper_bound(sv.begin(), sv.end(), x) - std::lower_bound(sv.begin(), sv.end(), x)); };
    int d0 = 0, c0 = cnt(A.inv[0]);
    for (int d = 1; d < D; d++) { int c = cnt(A.inv[d]); if (c < c0) { c0 = c; d0 = d; } }
    std::vector<int> m(D), mb(D), st;
    for (int e0 = 0; e0 < D; e0++) {
        if (B.inv[e0] != A.inv[d0]) continue;
        std::fill(m.begin(), m.end(), -1); std::fill(mb.begin(), mb.end(), -1);
        st.assign(1, d0); m[d0] = e0; mb[e0] = d0; bool ok = true; int nm = 1;
        while (!st.empty() && ok) {
            int d = st.back(); st.pop_back(); int e = m[d];
            int pairs[2][2] = {{A.nxt[d], B.nxt[e]}, {A.rev[d], B.rev[e]}};
            for (auto& pr : pairs) {
                int dd = pr[0], ee = pr[1];
                if (dd < 0) { if (ee >= 0) ok = false; if (!ok) break; continue; }
                if (m[dd] >= 0) { if (m[dd] != ee) { ok = false; break; } }
                else if (ee < 0 || mb[ee] >= 0) { ok = false; break; }
                else { m[dd] = ee; mb[ee] = dd; nm++; st.push_back(dd); }
            }
        }
        if (!ok || nm != D) continue;
        std::vector<int> fm(nfa, -1);
        for (int d = 0; d < D; d++) fm[A.face[d]] = B.face[m[d]];
        if (std::find(out.begin(), out.end(), fm) == out.end()) out.push_back(fm);
    }
}

struct Entry { std::vector<int> fv, fo; int nf = 0; Map M; std::vector<uint64_t> sets; };   // sets: bitmasks over faces
static std::vector<Entry> table;
static std::unordered_map<std::vector<int>, std::vector<int>, KeyHash> by_key;

// readings of cell i (bitmasks over its faces), or unseen (topology not in the table)
static bool options(int i, std::vector<uint64_t>& out) {
    const CellF& c = cells[i];
    out.clear();
    int nf = c.nf();
    if (nf > 64) return false;
    Map M; make_map(c.v.data(), c.o.data(), nf, M);
    auto it = by_key.find(map_key(c.o.data(), nf, M));
    if (it == by_key.end()) return false;
    bool seen = false, mm = false; Map Mm;
    std::vector<int> rv, ro;                         // mirror image: every face reversed
    std::vector<std::vector<int>> I;
    for (int e : it->second) {
        const Entry& E = table[e];
        isos(M, nf, E.M, E.nf, I);
        if (I.empty()) {
            if (!mm) {
                rv.resize(c.v.size()); ro.assign(c.o.begin(), c.o.end());
                for (int f = 0; f < nf; f++) std::reverse_copy(c.fc(f), c.fc(f) + c.fs(f), rv.begin() + c.o[f]);
                make_map(rv.data(), ro.data(), nf, Mm); mm = true;
            }
            isos(Mm, nf, E.M, E.nf, I);
        }
        if (I.empty()) continue;
        seen = true;
        for (auto& fm : I)
            for (uint64_t S : E.sets) {
                uint64_t mask = 0;
                for (int f = 0; f < nf; f++) if ((S >> fm[f]) & 1) mask |= 1ULL << f;
                out.push_back(mask);
            }
    }
    std::sort(out.begin(), out.end()); out.erase(std::unique(out.begin(), out.end()), out.end());
    return seen;
}

// ------------------------------------------------------------------ Delaunay complex (dual of the tessellation)
static std::vector<std::array<int, 4>> tets;           // sorted ids, sorted list
static std::vector<std::vector<int>> dnb;              // Delaunay neighbours (sorted)
static int orient(int v, int a, int b, int c) {         // +1/-1 if (v,a,b,c) is a tetrahedron, else 0
    std::array<int, 4> k = {v, a, b, c}; std::sort(k.begin(), k.end());
    if (!std::binary_search(tets.begin(), tets.end(), k)) return 0;
    double d = det3(sub(a, v), sub(b, v), sub(c, v));
    return d > 0 ? 1 : (d < 0 ? -1 : 0);
}
static bool has(const std::vector<int>& s, int x) { return std::binary_search(s.begin(), s.end(), x); }

// ------------------------------------------------------------------ templates and labelling
struct Kind {
    std::string name; bool ambiguous, perfect; std::vector<V3> V; bool adj[12][12]; std::vector<std::array<int, 3>> tri;
    std::vector<int> ehist;                     // histogram over edges of the number of common neighbours of the ends
    std::vector<std::array<int, 12>> aut;       // automorphisms of the template graph
    std::vector<int> acls;                      // per automorphism: its class, where a1 ~ a2 if a point-group element
                                                // R has V[a2(x)] = R V[a1(x)] for every x (the first of the class)
};
static std::vector<Kind> kinds;
struct Lab { int kind = -1; std::vector<std::pair<int, int>> m; int cls = -1; };   // m: (neighbour atom, template index),
                                                                     // sorted by atom; cls: class of a candidate
static std::vector<Lab> lab;                                       // labelling of each atom (kind -1: none)
// candidates of an ambiguous shell, compactly: per kind found, the first isomorphism as template indices in the order
// of the shell atoms (base, 12 per kind), and per candidate the automorphism ai of the template it is composed with
struct Cand { int16_t kind, ai; int32_t cls; int8_t b; };          // b: which base (of 12) in the set
struct CandSet { std::vector<int8_t> base; std::vector<Cand> c; };
static std::vector<std::vector<int>> lnb;                          // first neighbours (Delaunay minus diagonals)

static const V3& vecof(const Lab& L, int atom, bool& found) {
    static V3 zero{0, 0, 0};
    auto it = std::lower_bound(L.m.begin(), L.m.end(), std::make_pair(atom, -1));
    found = it != L.m.end() && it->first == atom;
    return found ? kinds[L.kind].V[it->second] : zero;
}

// backtracking extension of a partial map shell -> template (fcc_topo._extend)
static bool extend(std::vector<std::pair<int, int>>& mp, const std::vector<int>& sh, const std::vector<std::vector<int>>& adj, const Kind& K) {
    auto mapped_of = [&](int u) { for (auto& p : mp) if (p.first == u) return p.second; return -1; };
    std::vector<int> todo; for (int u : sh) if (mapped_of(u) < 0) todo.push_back(u);
    std::vector<int> order; std::set<int> seen; for (auto& p : mp) seen.insert(p.first);
    while (order.size() < todo.size()) {
        int best = -1, bs = -1;
        for (int u : todo) if (!seen.count(u)) {
            int s = 0; for (int z : adj[std::find(sh.begin(), sh.end(), u) - sh.begin()]) s += seen.count(z);
            if (s > bs) { bs = s; best = u; }
        }
        order.push_back(best); seen.insert(best);
    }
    std::vector<char> used(12, 0); for (auto& p : mp) used[p.second] = 1;
    std::function<bool(size_t)> rec = [&](size_t k) -> bool {
        if (k == order.size()) return true;
        int u = order[k]; const auto& au = adj[std::find(sh.begin(), sh.end(), u) - sh.begin()];
        for (int c = 0; c < 12; c++) {
            if (used[c]) continue;
            bool ok = true;
            for (int z : au) { int q = mapped_of(z); if (q >= 0 && !K.adj[q][c]) { ok = false; break; } }
            if (!ok) continue;
            for (auto& p : mp) if (std::find(au.begin(), au.end(), p.first) == au.end() && K.adj[p.second][c]) { ok = false; break; }
            if (!ok) continue;
            mp.push_back({u, c}); used[c] = 1;
            if (rec(k + 1)) return true;
            mp.pop_back(); used[c] = 0;
        }
        return false;
    };
    return rec(0);
}

// orientation-preserving shell isomorphisms of atom v onto kind K (all of them, or the first)
// what the labelling of atom v's shell needs whatever the kind: the shell graph (each of the 12 atoms with its 4
// neighbours in the shell, else bad), the histogram of common neighbours over its edges, and the orientation
// triangle (computed when first needed: o = -2)
struct ShellCtx {
    bool bad = false;
    std::vector<std::vector<int>> adj;
    std::vector<int> hist;
    int a0 = -1, b0 = -1, c0 = -1, o = -2;
};
static void shell_ctx(int v, ShellCtx& S) {
    const std::vector<int>& sh = lnb[v];
    S.adj.assign(12, {});
    for (int a = 0; a < 12; a++) { for (int b : sh) if (has(lnb[sh[a]], b)) S.adj[a].push_back(b); if (S.adj[a].size() != 4) { S.bad = true; return; } }
    S.hist.assign(13, 0);
    for (int a = 0; a < 12; a++)
        for (int b : S.adj[a]) {
            int ib = (int)(std::find(sh.begin(), sh.end(), b) - sh.begin());
            if (ib <= a) continue;
            int cm = 0; for (int z : S.adj[a]) cm += std::find(S.adj[ib].begin(), S.adj[ib].end(), z) != S.adj[ib].end();
            S.hist[cm]++;
        }
}
static std::vector<Lab> shell_isos(int v, int kind, bool all, ShellCtx& S, CandSet* cs = nullptr) {
    std::vector<Lab> out;
    const Kind& K = kinds[kind];
    const std::vector<int>& sh = lnb[v];
    if (S.bad) return out;
    // the shell graph and the template must have the same histogram of common neighbours over edges (an isomorphism
    // invariant): if not, no isomorphism exists and none of the searches below can succeed
    if (S.hist != K.ehist) return out;
    const auto& adj = S.adj;
    // the orientation is fixed by one triangle of the shell that is a Delaunay tetrahedron with v (the first, in
    // order of the neighbours; a triangle of the shell graph need not be one)
    if (S.o == -2) {
        S.o = 0;
        for (int ia = 0; ia < 12 && !S.o; ia++)
            for (int b : adj[ia]) {
                const auto& ab = adj[std::find(sh.begin(), sh.end(), b) - sh.begin()];
                for (int c : adj[ia]) if (std::find(ab.begin(), ab.end(), c) != ab.end()) {
                    int oo = orient(v, sh[ia], b, c);
                    if (oo) { S.a0 = sh[ia]; S.b0 = b; S.c0 = c; S.o = oo; break; }
                }
                if (S.o) break;
            }
    }
    int a0 = S.a0, b0 = S.b0, c0 = S.c0, o = S.o;
    if (!o) return out;
    // every isomorphism shell -> template is a.phi for one automorphism a of the template, phi the first one found;
    // each is fixed by the image of (a0, b0, c0), so the search from triangle t, order q succeeds exactly for the
    // images of a.phi with the right orientation, and gives a.phi. All of them, in the order of (t, q):
    auto all_from = [&](const Lab& phi) {
        std::vector<int> ph(12, -1);
        for (auto& pr : phi.m) ph[std::find(sh.begin(), sh.end(), pr.first) - sh.begin()] = pr.second;
        int ia0 = (int)(std::find(sh.begin(), sh.end(), a0) - sh.begin()), ib0 = (int)(std::find(sh.begin(), sh.end(), b0) - sh.begin()),
            ic0 = (int)(std::find(sh.begin(), sh.end(), c0) - sh.begin());
        const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
        std::vector<std::pair<int, int>> res;          // ((t, q) key, automorphism)
        for (size_t ai = 0; ai < K.aut.size(); ai++) {
            const auto& al = K.aut[ai];
            int A1 = al[ph[ia0]], B1 = al[ph[ib0]], C1 = al[ph[ic0]];
            double d = det3(K.V[A1], K.V[B1], K.V[C1]);
            int sg = d > 0 ? 1 : (d < 0 ? -1 : 0);
            if (sg != o) continue;
            std::array<int, 3> srt = {A1, B1, C1}; std::sort(srt.begin(), srt.end());
            int t = (int)(std::lower_bound(K.tri.begin(), K.tri.end(), srt) - K.tri.begin());
            int q = 0; for (; q < 6; q++) if (srt[perms[q][0]] == A1 && srt[perms[q][1]] == B1 && srt[perms[q][2]] == C1) break;
            res.push_back({t * 6 + q, (int)ai});
        }
        std::sort(res.begin(), res.end());              // keys are distinct (one automorphism per image of a dart)
        int8_t b = (int8_t)(cs->base.size() / 12);
        for (int a = 0; a < 12; a++) cs->base.push_back((int8_t)ph[a]);
        for (auto& r : res) cs->c.push_back(Cand{(int16_t)kind, (int16_t)r.second, kind * 4096 + K.acls[r.second], b});
    };
    for (auto& t : K.tri) {
        int p[3] = {t[0], t[1], t[2]};
        const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
        for (auto& q : perms) {
            int A1 = p[q[0]], B1 = p[q[1]], C1 = p[q[2]];
            double d = det3(K.V[A1], K.V[B1], K.V[C1]);
            int s = d > 0 ? 1 : (d < 0 ? -1 : 0);
            if (s != o) continue;
            std::vector<std::pair<int, int>> mp{{a0, A1}, {b0, B1}, {c0, C1}};
            if (extend(mp, sh, adj, K)) {
                Lab L; L.kind = kind; L.m = mp; std::sort(L.m.begin(), L.m.end());
                if (all) { all_from(L); return out; }
                out.push_back(L);
                return out;
            }
        }
    }
    return out;
}

// candidate c of atom v as a labelling (the shell atoms, lnb[v], are sorted, so m comes out sorted by atom)
static void make_lab(int v, const CandSet& cs, const Cand& c, Lab& L) {
    const auto& sh = lnb[v]; const auto& al = kinds[c.kind].aut[c.ai];
    L.kind = c.kind; L.cls = c.cls; L.m.resize(12);
    for (int a = 0; a < 12; a++) L.m[a] = {sh[a], al[cs.base[c.b * 12 + a]]};
}

// gluing x_i = R x_j + d between labelled atoms: first neighbours from their common neighbours, else (opposite corners
// of an octahedral hole) from the middle square; the first point-group element that maps the points (fcc_topo.g)
static int gluing(int i, int j, const Lab& Li, const Lab& Lj, V3& dout) {
    bool fij, fji;
    const V3& dij = vecof(Li, j, fij); const V3& dji = vecof(Lj, i, fji);
    V3 Pi[13], Pj[13], d; size_t np = 0;              // at most 12 common neighbours, and one more point
    if (fij && fji) {
        for (auto& p : Li.m) { bool f; const V3& w = vecof(Lj, p.first, f); if (f) { Pi[np] = kinds[Li.kind].V[p.second]; Pj[np++] = w; } }
        Pi[np] = V3{0, 0, 0}; Pj[np++] = dji; d = dij;
    } else if (!fij) {
        for (auto& p : Li.m) { bool f; const V3& w = vecof(Lj, p.first, f); if (f) { Pi[np] = kinds[Li.kind].V[p.second]; Pj[np++] = w; } }
        if (np != 4) return -1;
        V3 ci{0, 0, 0}, cj{0, 0, 0};
        for (int k = 0; k < 4; k++) for (int r = 0; r < 3; r++) { ci[r] += Pi[k][r] / 4; cj[r] += Pj[k][r] / 4; }
        d = V3{2 * ci[0], 2 * ci[1], 2 * ci[2]};
        Pi[np] = V3{0, 0, 0}; Pj[np++] = V3{2 * cj[0], 2 * cj[1], 2 * cj[2]};
    } else return -1;
    double mx = 1.0; for (size_t k = 0; k < np; k++) for (int r = 0; r < 3; r++) mx = std::max(mx, std::fabs(Pi[k][r]));
    double tol = 1e-6 * mx;
    for (int g = 0; g < (int)PG.size(); g++) {
        bool ok = true;
        for (size_t k = 0; k < np && ok; k++) {
            V3 y = mvec(PG[g], Pj[k]);
            for (int r = 0; r < 3; r++) if (!(std::fabs(y[r] + d[r] - Pi[k][r]) < tol)) { ok = false; break; }
        }
        if (ok) { dout = d; return g; }
    }
    return -1;
}

int main(int argc, char** argv) {
    if (argc < 5) { fprintf(stderr, "usage: vb_closepacked params table sample prefix\n"); return 1; }
    double T0 = omp_get_wtime(), Tl = T0;
    auto tick = [&](const char* w) { if (getenv("VB_TIMING")) { double t = omp_get_wtime(); fprintf(stderr, "%-18s %.3f s\n", w, t - Tl); Tl = t; } };
    // ---------------- parameters
    FILE* fp = fopen(argv[1], "r");
    int G, K;
    if (fscanf(fp, "%d", &G) != 1) return 2;
    PG.resize(G);
    for (int g = 0; g < G; g++) for (int k = 0; k < 9; k++) if (fscanf(fp, "%lf", &PG[g][k]) != 1) return 2;
    if (fscanf(fp, "%d", &K) != 1) return 2;
    kinds.resize(K);
    for (int k = 0; k < K; k++) {
        char nm[64]; int amb, prf;
        if (fscanf(fp, "%63s %d %d", nm, &amb, &prf) != 3) return 2;
        kinds[k].name = nm; kinds[k].ambiguous = amb; kinds[k].perfect = prf; kinds[k].V.resize(12);
        for (int a = 0; a < 12; a++) for (int r = 0; r < 3; r++) if (fscanf(fp, "%lf", &kinds[k].V[a][r]) != 1) return 2;
        for (int a = 0; a < 12; a++) for (int b = 0; b < 12; b++) { int x; if (fscanf(fp, "%d", &x) != 1) return 2; kinds[k].adj[a][b] = x; }
        for (int a = 0; a < 12; a++) for (int b = a + 1; b < 12; b++) for (int c = b + 1; c < 12; c++)
            if (kinds[k].adj[a][b] && kinds[k].adj[a][c] && kinds[k].adj[b][c]) kinds[k].tri.push_back({a, b, c});
        Kind& Kk = kinds[k];
        Kk.ehist.assign(13, 0);
        for (int a = 0; a < 12; a++) for (int b = a + 1; b < 12; b++) if (Kk.adj[a][b]) {
            int cm = 0; for (int z = 0; z < 12; z++) cm += Kk.adj[a][z] && Kk.adj[b][z];
            Kk.ehist[cm]++;
        }
        // automorphisms of the template graph (backtracking over the images of 0, 1, ..., 11)
        std::array<int, 12> pm; std::array<char, 12> usd{};
        std::function<void(int)> rec = [&](int a) {
            if (a == 12) { Kk.aut.push_back(pm); return; }
            for (int c = 0; c < 12; c++) {
                if (usd[c]) continue;
                bool okk = true;
                for (int z = 0; z < a && okk; z++) if (Kk.adj[a][z] != Kk.adj[c][pm[z]]) okk = false;
                if (!okk) continue;
                pm[a] = c; usd[c] = 1; rec(a + 1); usd[c] = 0;
            }
        };
        rec(0);
        Kk.acls.assign(Kk.aut.size(), -1);
        for (size_t a2 = 0; a2 < Kk.aut.size(); a2++)
            for (size_t a1 = 0; a1 <= a2 && Kk.acls[a2] < 0; a1++) {
                if (Kk.acls[a1] != (int)a1 && a1 != a2) continue;           // compare with the first of each class
                for (auto& R : PG) {
                    bool ok = true;
                    for (int x = 0; x < 12 && ok; x++) { V3 y = mvec(R, Kk.V[Kk.aut[a1][x]]); for (int r = 0; r < 3; r++) if (std::fabs(y[r] - Kk.V[Kk.aut[a2][x]][r]) > 1e-6) ok = false; }
                    if (ok) { Kk.acls[a2] = (int)a1; break; }
                }
            }
    }
    fclose(fp);
    // multiplication table of the point group (products are elements)
    std::vector<std::vector<int>> gmul(G, std::vector<int>(G, -1));
    for (int a = 0; a < G; a++) for (int b = 0; b < G; b++) {
        M3 p = mul(PG[a], PG[b]);
        for (int c = 0; c < G; c++) { double e = 0; for (int k = 0; k < 9; k++) e = std::max(e, std::fabs(p[k] - PG[c][k])); if (e < 1e-6) { gmul[a][b] = c; break; } }
    }
    int gid = -1;
    for (int c = 0; c < G; c++) { double e = 0; for (int k = 0; k < 9; k++) e = std::max(e, std::fabs(PG[c][k] - (k % 4 == 0 ? 1.0 : 0.0))); if (e < 1e-9) gid = c; }
    // ---------------- table
    {
        FILE* ft = fopen(argv[2], "rb");
        int32_t E; if (fread(&E, 4, 1, ft) != 1) return 2;
        table.resize(E);
        for (int e = 0; e < E; e++) {
            int32_t nf; if (fread(&nf, 4, 1, ft) != 1) return 2;
            Entry& T = table[e];
            T.nf = nf; T.fo.assign(1, 0);
            for (int f = 0; f < nf; f++) {
                int32_t n; if (fread(&n, 4, 1, ft) != 1) return 2;
                size_t a = T.fv.size(); T.fv.resize(a + n); if (fread(T.fv.data() + a, 4, n, ft) != (size_t)n) return 2;
                T.fo.push_back((int)T.fv.size());
            }
            if (nf > 64) { fprintf(stderr, "table entry with more than 64 faces\n"); return 2; }
            int32_t ns; if (fread(&ns, 4, 1, ft) != 1) return 2;
            T.sets.assign(ns, 0);
            for (int s = 0; s < ns; s++) {
                int32_t m; if (fread(&m, 4, 1, ft) != 1) return 2;
                std::vector<int32_t> S(m); if (fread(S.data(), 4, m, ft) != (size_t)m) return 2;
                for (int x : S) T.sets[s] |= 1ULL << x;
            }
            make_map(T.fv.data(), T.fo.data(), nf, T.M);
            by_key[map_key(T.fo.data(), nf, T.M)].push_back(e);
        }
        fclose(ft);
    }
    tick("read_table");
    // ---------------- sample and tessellation
    {
        FILE* fs = fopen(argv[3], "rb");
        int32_t n; if (fread(&n, 4, 1, fs) != 1) return 3;
        N = n;
        if (fread(Lb, 8, 3, fs) != 3) return 3;
        int32_t pp[3]; if (fread(pp, 4, 3, fs) != 3) return 3;
        for (int k = 0; k < 3; k++) per[k] = pp[k];
        X.resize(N);
        if (fread(X.data(), 8, (size_t)3 * N, fs) != (size_t)3 * N) return 3;
        fclose(fs);
    }
    int threads = omp_get_max_threads();
    cells.resize(N);
    {
        double dens = N / (Lb[0] * Lb[1] * Lb[2]);
        int nx = std::max(1, int(Lb[0] * pow(dens / 5, 1. / 3))), ny = std::max(1, int(Lb[1] * pow(dens / 5, 1. / 3))),
            nz = std::max(1, int(Lb[2] * pow(dens / 5, 1. / 3)));
        container_3d con(0, Lb[0], 0, Lb[1], 0, Lb[2], nx, ny, nz, per[0], per[1], per[2], 8, threads);
        for (int i = 0; i < N; i++) con.put(i, X[i][0], X[i][1], X[i][2]);
        #pragma omp parallel num_threads(threads)
        {
            voronoicell_neighbor_3d vc(con);
            std::vector<int> nb, fv;
            #pragma omp for schedule(dynamic, 256)
            for (auto cli = con.begin(); cli < con.end(); ++cli)
                if (con.compute_cell(vc, cli)) {
                    int pid = con.id[cli->ijk][cli->q];
                    vc.neighbors(nb); vc.face_vertices(fv);
                    CellF& c = cells[pid];
                    c.nbr = nb; c.v.clear(); c.o.assign(1, 0);
                    size_t k = 0;
                    for (size_t f = 0; f < nb.size(); f++) { int m = fv[k]; c.v.insert(c.v.end(), fv.begin() + k + 1, fv.begin() + k + 1 + m); c.o.push_back((int)c.v.size()); k += m + 1; }
                    c.ok = true;
                }
        }
    }
    tick("tessellate");
    // Delaunay tetrahedra: the Voronoi vertices shared by four cells, each found from its smallest atom (whose cell
    // lists them in order of the atom, then sorted: the whole list comes out sorted, with no global sort)
    {
        auto cell_tets = [&](int i, std::vector<std::array<int, 4>>& out) {
            out.clear();
            const CellF& c = cells[i];
            if (!c.ok) return;
            int nv = 0; for (int x : c.v) nv = std::max(nv, x + 1);
            std::vector<int> cnt(nv, 0); std::vector<std::array<int, 3>> nb3(nv);
            for (int f = 0; f < c.nf(); f++) for (int k = 0; k < c.fs(f); k++) { int x = c.fc(f)[k]; if (cnt[x] < 3) nb3[x][cnt[x]] = c.nbr[f]; cnt[x]++; }
            for (int x = 0; x < nv; x++) {
                if (cnt[x] != 3) continue;
                std::array<int, 4> t = {i, nb3[x][0], nb3[x][1], nb3[x][2]};
                if (t[1] < 0 || t[2] < 0 || t[3] < 0) continue;
                std::sort(t.begin(), t.end());
                if (t[0] == t[1] || t[1] == t[2] || t[2] == t[3]) continue;
                if (t[0] == i) out.push_back(t);
            }
            std::sort(out.begin(), out.end()); out.erase(std::unique(out.begin(), out.end()), out.end());
        };
        std::vector<int64_t> to(N + 1, 0);
        #pragma omp parallel
        {
            std::vector<std::array<int, 4>> buf;
            #pragma omp for schedule(dynamic, 256)
            for (int i = 0; i < N; i++) { cell_tets(i, buf); to[i + 1] = (int64_t)buf.size(); }
        }
        for (int i = 0; i < N; i++) to[i + 1] += to[i];
        tets.resize(to[N]);
        #pragma omp parallel
        {
            std::vector<std::array<int, 4>> buf;
            #pragma omp for schedule(dynamic, 256)
            for (int i = 0; i < N; i++) { cell_tets(i, buf); std::copy(buf.begin(), buf.end(), tets.begin() + to[i]); }
        }
        // Delaunay neighbours: the other three atoms of each tetrahedron, sorted, without repeats
        std::vector<int64_t> no(N + 1, 0);
        for (auto& t : tets) for (int a = 0; a < 4; a++) no[t[a] + 1] += 3;
        for (int i = 0; i < N; i++) no[i + 1] += no[i];
        std::vector<int> nall(no[N]);
        { std::vector<int64_t> fill(no.begin(), no.end() - 1);
          for (auto& t : tets) for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) if (a != b) nall[fill[t[a]]++] = t[b]; }
        dnb.assign(N, {});
        #pragma omp parallel for schedule(dynamic, 1024)
        for (int i = 0; i < N; i++) {
            std::sort(nall.begin() + no[i], nall.begin() + no[i + 1]);
            auto e = std::unique(nall.begin() + no[i], nall.begin() + no[i + 1]);
            dnb[i].assign(nall.begin() + no[i], e);
        }
    }
    tick("delaunay");
    // ---------------- (a) readings and agreement
    std::vector<int8_t> st(N, OK);
    std::vector<char> incell(N, 0), seenc(N, 0);
    std::vector<std::vector<uint64_t>> dom0(N);
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < N; i++) {
        const CellF& c = cells[i];
        bool bad = !c.ok;
        if (!bad) {
            std::vector<int> s = c.nbr;
            for (int x : s) if (x < 0) bad = true;
            std::sort(s.begin(), s.end());
            if (!bad && std::adjacent_find(s.begin(), s.end()) != s.end()) bad = true;
            if (!bad && s != dnb[i]) bad = true;
        }
        if (bad) { st[i] = NOCELL; continue; }
        incell[i] = 1;
        if (options(i, dom0[i])) seenc[i] = 1; else st[i] = UNSEEN;
    }
    tick("readings");
    auto face_of = [&](int i, int j) { const auto& nb = cells[i].nbr; for (size_t f = 0; f < nb.size(); f++) if (nb[f] == j) return (int)f; return -1; };
    std::vector<char> dead(N, 0), mark(N, 0);
    std::vector<std::vector<uint64_t>> dom(N);
    // ff: the first face of i with the neighbour across face f (gluings and agreement use it); rv: the first face of
    // that neighbour with i
    std::vector<int64_t> aoff(N + 1, 0);
    for (int i = 0; i < N; i++) aoff[i + 1] = aoff[i] + cells[i].nf();
    std::vector<int> ff(aoff[N]), rv(aoff[N]);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < N; i++)
        for (int f = 0; f < cells[i].nf(); f++) {
            int j = cells[i].nbr[f];
            ff[aoff[i] + f] = face_of(i, j);
            rv[aoff[i] + f] = j >= 0 ? face_of(j, i) : -1;
        }
    // arc consistency in rounds, each taking all its decisions from the readings at its start (so the result depends
    // on no order, and the cells of a round are independent of one another); a cell left with none no longer
    // constrains its neighbours in that pass, and once the pass ends it drops out and agreement starts again
    while (true) {
        std::vector<int> act, next;
        for (int i = 0; i < N; i++) { dom[i].clear(); if (seenc[i] && !dead[i]) { dom[i] = dom0[i]; act.push_back(i); } }
        while (!act.empty()) {
            std::vector<std::vector<uint64_t>> nd(act.size());
            std::vector<char> changed(act.size(), 0);
            #pragma omp parallel
            {
                std::vector<uint64_t> keep;
                #pragma omp for schedule(dynamic, 256)
                for (size_t t = 0; t < act.size(); t++) {
                    int i = act[t];
                    keep.assign(dom[i].begin(), dom[i].end());
                    for (int f = 0; f < cells[i].nf() && !keep.empty(); f++) {
                        int j = cells[i].nbr[f];
                        if (j < 0 || dom[j].empty() || ff[aoff[i] + f] != f) continue;
                        int fji = rv[aoff[i] + f];
                        bool sup[2] = {false, false};
                        for (uint64_t r : dom[j]) sup[(r >> fji) & 1] = true;
                        size_t w = 0;
                        for (uint64_t r : keep) if (sup[(r >> f) & 1]) keep[w++] = r;
                        keep.resize(w);
                    }
                    if (keep.size() != dom[i].size()) { changed[t] = 1; nd[t] = keep; }
                }
            }
            next.clear();
            for (size_t t = 0; t < act.size(); t++) if (changed[t]) dom[act[t]].swap(nd[t]);
            for (size_t t = 0; t < act.size(); t++) if (changed[t])
                for (int k : cells[act[t]].nbr) if (k >= 0 && !dom[k].empty() && !mark[k]) { mark[k] = 1; next.push_back(k); }
            for (int k : next) mark[k] = 0;
            act.swap(next);
        }
        bool any = false;
        for (int i = 0; i < N; i++) if (seenc[i] && !dead[i] && dom[i].empty()) { dead[i] = 1; any = true; }
        if (!any) break;
    }
    std::vector<int>().swap(ff); std::vector<int>().swap(rv);
    for (int i = 0; i < N; i++) if (dead[i]) dom[i].clear();
    tick("agreement");
    // decided cells (one reading) and their first neighbours
    std::vector<std::vector<int>> prim(N);
    std::vector<char> isprim(N, 0);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < N; i++) {
        if (!seenc[i]) continue;
        if (dom[i].size() == 1) {
            for (int f = 0; f < cells[i].nf(); f++) if ((dom[i][0] >> f) & 1) prim[i].push_back(cells[i].nbr[f]);
            std::sort(prim[i].begin(), prim[i].end()); isprim[i] = 1;
        } else st[i] = dom[i].empty() ? DISAGREE : AMBIG;
    }
    // a decided pair that disagrees about their face marks both; the test is symmetric, so each cell marks itself
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < N; i++) if (isprim[i])
        for (int j : dnb[i]) if (isprim[j] && has(prim[i], j) != has(prim[j], i)) { st[i] = DISAGREE; break; }
    for (int i = 0; i < N; i++) if (isprim[i] && st[i] != OK) isprim[i] = 0;
    // first neighbours: Delaunay neighbours minus the diagonals (edges a decided endpoint excludes)
    lnb.assign(N, {});
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int a = 0; a < N; a++)
        for (int b : dnb[a]) {
            bool diag = (isprim[a] && !has(prim[a], b)) || (isprim[b] && !has(prim[b], a));
            if (!diag) lnb[a].push_back(b);
        }
    // (b) candidates: twelve neighbours, every Delaunay edge among them decided by a decided endpoint
    std::vector<char> cand(N, 0);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < N; i++) {
        if (!isprim[i]) continue;
        const auto& S = prim[i];
        if (S.size() != 12) { st[i] = SHELL; continue; }
        bool und = false;
        for (size_t x = 0; x < S.size() && !und; x++) for (size_t y = x + 1; y < S.size() && !und; y++) {
            int a = S[x], b = S[y];
            if (has(dnb[a], b) && !isprim[a] && !isprim[b]) und = true;
        }
        if (und) { st[i] = UNDECIDED_SHELL; continue; }
        cand[i] = 1;
    }
    std::vector<std::vector<uint64_t>>().swap(dom0); std::vector<std::vector<uint64_t>>().swap(dom);
    std::vector<std::vector<int>>().swap(prim);
    // labelling: the first unambiguous kind that fits; ambiguous kinds collect every isomorphism
    lab.assign(N, Lab());
    std::vector<CandSet> cands(N);
    #pragma omp parallel for schedule(dynamic, 64)
    for (int v = 0; v < N; v++) {
        if (!cand[v] || lnb[v].size() != 12) continue;
        CandSet c_;
        ShellCtx S; shell_ctx(v, S);
        for (int k = 0; k < K; k++) {
            if (kinds[k].ambiguous) { shell_isos(v, k, true, S, &c_); continue; }
            auto one = shell_isos(v, k, false, S);
            if (!one.empty()) { lab[v] = one[0]; c_ = CandSet(); break; }
        }
        if (!c_.c.empty() && lab[v].kind < 0) { c_.base.shrink_to_fit(); c_.c.shrink_to_fit(); cands[v] = std::move(c_); }
    }
    tick("shells");
    // ambiguous shells: grow from the cell with the most labelled neighbours (ties: smallest index), picking the
    // candidate that glues to all of them
    {
        // cells with candidates by (labelled first neighbours, -index), largest first: an indexed binary heap whose
        // keys only increase (the same order as a heap of entries re-pushed when outdated). First-neighbour lists are
        // symmetric, so labelling v adds one to the count of each of its first neighbours
        std::vector<int> nfx(N, 0), hpos(N, -1), hp;
        #pragma omp parallel for schedule(dynamic, 1024)
        for (int v = 0; v < N; v++) if (!cands[v].c.empty()) { int s2 = 0; for (int w : lnb[v]) s2 += lab[w].kind >= 0; nfx[v] = s2; }
        auto key = [&](int v) { return ((int64_t)nfx[v] << 32) | (uint32_t)(0x7fffffff - v); };
        auto up = [&](int t) {
            int j = hp[t]; int64_t kj = key(j);
            while (t > 0) { int u = (t - 1) / 2; if (key(hp[u]) >= kj) break; hp[t] = hp[u]; hpos[hp[t]] = t; t = u; }
            hp[t] = j; hpos[j] = t;
        };
        auto pop = [&]() {
            int top = hp[0]; hpos[top] = -1;
            int last = hp.back(); hp.pop_back();
            if (!hp.empty()) {
                int t = 0, n = (int)hp.size(); int64_t kl = key(last);
                while (true) {
                    int c = 2 * t + 1; if (c >= n) break;
                    if (c + 1 < n && key(hp[c + 1]) > key(hp[c])) c++;
                    if (key(hp[c]) <= kl) break;
                    hp[t] = hp[c]; hpos[hp[t]] = t; t = c;
                }
                hp[t] = last; hpos[last] = t;
            }
            return top;
        };
        for (int v = 0; v < N; v++) if (!cands[v].c.empty()) { hp.push_back(v); hpos[v] = (int)hp.size() - 1; up(hpos[v]); }
        std::vector<int> tried;                      // classes of candidates already tested
        Lab Lt, Lb;
        while (!hp.empty()) {
            int v = pop();
            std::vector<int> nbl; for (int w : lnb[v]) if (lab[w].kind >= 0) nbl.push_back(w);
            int best = -1, bok = -1;
            tried.clear();
            for (size_t t = 0; t < cands[v].c.size(); t++) {
                // a candidate R-related to one already tested (R in the point group) glues to exactly the same
                // neighbours (with R g in place of g), so its count is equal and it cannot be picked: skipped
                int cl = cands[v].c[t].cls;
                if (cl >= 0) { if (std::find(tried.begin(), tried.end(), cl) != tried.end()) continue; tried.push_back(cl); }
                make_lab(v, cands[v], cands[v].c[t], Lt);
                int ok = 0; V3 d;
                for (int w : nbl) {
                    bool f1, f2; vecof(Lt, w, f1); vecof(lab[w], v, f2);
                    if (f1 && f2 && gluing(v, w, Lt, lab[w], d) >= 0) ok++;
                }
                if (best < 0 || ok > bok) { best = (int)t; bok = ok; Lb = Lt; }
                if (ok == (int)nbl.size() && !nbl.empty()) break;
            }
            if (best < 0 || (!nbl.empty() && bok < (int)nbl.size())) continue;
            lab[v] = Lb;
            for (int w : lnb[v]) { nfx[w]++; if (hpos[w] >= 0) up(hpos[w]); }
        }
    }
    if (getenv("VB_DEBUG")) {
        long nc = 0, n12 = 0, nwith = 0, nlab = 0; std::map<int, long> byk, ncand;
        for (int v = 0; v < N; v++) { if (!cand[v]) continue; nc++; n12 += lnb[v].size() == 12; nwith += !cands[v].c.empty() || lab[v].kind >= 0; nlab += lab[v].kind >= 0; if (lab[v].kind >= 0) byk[lab[v].kind]++; ncand[(int)cands[v].c.size()]++; }
        fprintf(stderr, "cand %ld, with 12 neighbours %ld, with a labelling or candidates %ld, labelled %ld; by kind:", nc, n12, nwith, nlab);
        for (auto& kv : byk) fprintf(stderr, " %s:%ld", kinds[kv.first].name.c_str(), kv.second);
        fprintf(stderr, "; candidates per atom:"); for (auto& kv : ncand) fprintf(stderr, " %d:%ld", kv.first, kv.second); fprintf(stderr, "\n");
    }
    std::vector<CandSet>().swap(cands);
    for (int i = 0; i < N; i++) if (cand[i] && lab[i].kind < 0) st[i] = SHELL;
    tick("labelling");
    // (c) every Delaunay edge between labelled cells carries a gluing that is a crystal symmetry
    {
        std::vector<char> bad(N, 0);
        #pragma omp parallel for schedule(dynamic, 1024)
        for (int a = 0; a < N; a++) {
            if (lab[a].kind < 0) continue;
            for (int b : dnb[a]) if (b > a && lab[b].kind >= 0) {
                V3 d; if (gluing(a, b, lab[a], lab[b], d) < 0) { bad[a] = 1; bad[b] = 1; }
            }
        }
        for (int i = 0; i < N; i++) if (bad[i]) { st[i] = GLUE; lab[i] = Lab(); }
    }
    std::vector<char> kept(N, 0);
    int nkept = 0;
    for (int i = 0; i < N; i++) { kept[i] = lab[i].kind >= 0 && kinds[lab[i].kind].perfect; nkept += kept[i]; if (lab[i].kind < 0 && st[i] == OK) st[i] = SHELL; }
    tick("gluing_check");
    // ---------------- gluings between kept cells (Delaunay edges), rotation transport, regions
    std::vector<std::vector<int>> nbl(N);              // labelled first neighbours for kept atoms, else Delaunay
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < N; i++) {
        if (kept[i]) { for (auto& p : lab[i].m) nbl[i].push_back(p.first); }
        else nbl[i] = lnb[i].empty() ? dnb[i] : lnb[i];
    }
    std::vector<std::array<int, 3>> gl; std::vector<V3> gd;
    {
        std::vector<std::vector<std::array<int, 3>>> lg(threads); std::vector<std::vector<V3>> ld(threads);
        #pragma omp parallel num_threads(threads)
        {
            int t = omp_get_thread_num();
            #pragma omp for schedule(static)
            for (int a = 0; a < N; a++) {
                if (!kept[a]) continue;
                for (int b : dnb[a]) if (b > a && kept[b]) {
                    V3 d; int g = gluing(a, b, lab[a], lab[b], d);
                    if (g >= 0) { lg[t].push_back({a, b, g}); ld[t].push_back(d); }
                }
            }
        }
        size_t tot = 0; for (auto& x : lg) tot += x.size();
        gl.reserve(tot); gd.reserve(tot);
        for (int t = 0; t < threads; t++) {
            gl.insert(gl.end(), lg[t].begin(), lg[t].end()); gd.insert(gd.end(), ld[t].begin(), ld[t].end());
            std::vector<std::array<int, 3>>().swap(lg[t]); std::vector<V3>().swap(ld[t]);
        }
    }
    // lookup of the gluing of (i, j) from the sorted list (i < j stored; reverse = inverse)
    // the list is in order of (a, b) (static schedule: each thread a block of a, in order; b ascending): the
    // gluings of a are gs[a] .. gs[a + 1]
    std::vector<int64_t> gs(N + 1, 0);
    for (auto& e : gl) gs[e[0] + 1]++;
    for (int i = 0; i < N; i++) gs[i + 1] += gs[i];
    std::vector<int> ginv(G, -1);
    for (int a = 0; a < G; a++) for (int b = 0; b < G; b++) if (gmul[a][b] == gid) ginv[a] = b;
    auto rot_of = [&](int i, int j) -> int {             // rotation index of g(i, j), or -1
        int a = std::min(i, j), b = std::max(i, j);
        for (int64_t e = gs[a]; e < gs[a + 1]; e++) if (gl[e][1] == b) return i < j ? gl[e][2] : ginv[gl[e][2]];
        return -1;
    };
    std::vector<int> region(N, -1), rglob(N, -1);
    int nreg = 0, conflicts = 0;
    for (int root = 0; root < N; root++) {
        if (!kept[root] || rglob[root] >= 0) continue;
        rglob[root] = gid; region[root] = nreg;
        std::deque<int> dq{root};
        while (!dq.empty()) {
            int i = dq.front(); dq.pop_front();
            for (int j : nbl[i]) {
                if (!kept[j]) continue;
                int g = rot_of(i, j);
                if (g < 0) continue;
                int Rj = gmul[rglob[i]][g];
                if (rglob[j] >= 0) { if (rglob[j] != Rj) conflicts++; continue; }
                rglob[j] = Rj; region[j] = nreg; dq.push_back(j);
            }
        }
        nreg++;
    }
    tick("transport");
    // ---------------- output (as vb_engine)
    std::string pre(argv[4]);
    {
        FILE* fc = fopen((pre + ".cells").c_str(), "wb");
        for (int i = 0; i < N; i++) { int32_t row[4] = {(int32_t)kept[i], lab[i].kind, region[i], rglob[i]}; fwrite(row, 4, 4, fc); }
        fclose(fc);
        FILE* fn = fopen((pre + ".nbrs").c_str(), "wb");
        int32_t off = 0; fwrite(&off, 4, 1, fn);
        for (int i = 0; i < N; i++) { off += (int32_t)nbl[i].size(); fwrite(&off, 4, 1, fn); }
        for (int i = 0; i < N; i++) for (int j : nbl[i]) { int32_t jj = j; fwrite(&jj, 4, 1, fn); }
        fclose(fn);
        // Voronoi (face-sharing) neighbours of every atom: defect clusters and their shells are defined by them
        FILE* fv = fopen((pre + ".vnbrs").c_str(), "wb");
        off = 0; fwrite(&off, 4, 1, fv);
        for (int i = 0; i < N; i++) { off += (int32_t)dnb[i].size(); fwrite(&off, 4, 1, fv); }
        for (int i = 0; i < N; i++) for (int j : dnb[i]) { int32_t jj = j; fwrite(&jj, 4, 1, fv); }
        fclose(fv);
        FILE* fg = fopen((pre + ".glue").c_str(), "wb");
        for (size_t e = 0; e < gl.size(); e++) { int32_t r3[3] = {gl[e][0], gl[e][1], gl[e][2]}; fwrite(r3, 4, 3, fg); fwrite(gd[e].data(), 8, 3, fg); }
        fclose(fg);
        FILE* ff = fopen((pre + ".faces").c_str(), "wb");
        for (int i = 0; i < N; i++) {
            if (kept[i]) continue;
            int32_t hd[2] = {i, cells[i].ok ? (int32_t)cells[i].nf() : -1};
            fwrite(hd, 4, 2, ff);
            if (!cells[i].ok) continue;
            for (int f = 0; f < cells[i].nf(); f++) {
                int32_t fh[2] = {cells[i].nbr[f], (int32_t)cells[i].fs(f)};
                fwrite(fh, 4, 2, ff);
                fwrite(cells[i].fc(f), 4, cells[i].fs(f), ff);
            }
        }
        fclose(ff);
        FILE* fs = fopen((pre + ".cp").c_str(), "wb");
        fwrite(st.data(), 1, N, fs);
        fclose(fs);
    }
    int cnt[8] = {0}; for (int i = 0; i < N; i++) cnt[st[i]]++;
    printf("stats kept %d regions %d gluings %zu conflicts %d ok %d unseen %d ambiguous %d nocell %d disagree %d undecided_shell %d shell %d glue %d\n",
           nkept, nreg, gl.size(), conflicts, cnt[0], cnt[1], cnt[2], cnt[3], cnt[4], cnt[5], cnt[6], cnt[7]);
    tick("output");
    fflush(stdout); fflush(stderr);
    _exit(0);
}
