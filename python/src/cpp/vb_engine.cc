#include <unistd.h>
// vb_engine: the combinatorial connection of VoroBurgers in C++ (voro++ dev branch).
//
// Input 1 (reference file, written by voroburgers/engine.py):
//   line 1: NREF G tol
//   G lines: 9 numbers, the point-group matrices (row major)
//   per reference orbit:
//     "ref" species nfaces nverts mode        (mode: 0 rigid, 1 rigid with labels, 2 needs neighbours; +4: glue by
//                                             centers and common neighbours, for degenerate references)
//     nfaces lines:  sighash species n v_0 ... v_{n-1}   (neighbour label, 0 = any + face cycle, local vertex ids)
//     nfaces lines:  dx dy dz                             (ideal neighbour vector of each face)
//     nverts lines:  px py pz                             (ideal vertex positions)
// Input 2 (sample file): line 1 "N Lx Ly Lz px py pz use_species", then N lines "x y z [species]".
// Output (stdout), text:
//   "cells N"; then N lines: i kept orbit region rglob_index npseudo  (rglob into the point group)
//   "gluings M"; then M lines: i j Rindex dx dy dz   (x_i = R x_j + d), one per kept face pair i<j
//   "stats ..." summary line.
// Build: g++ -O2 -fopenmp -I<voro>/src vb_engine.cc -L<voro>/src -lvoro++ -o vb_engine

#include "voro++.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <map>
#include <atomic>
#include <queue>
#include <deque>
#include <set>
#include <string>
#include <vector>
#include <omp.h>
#include <parallel/algorithm>
using namespace voro;

typedef std::array<double, 3> V3;
typedef std::array<double, 9> M3;

// ------------------------------------------------------------------ cells
// half-edge table of a cell: (u, v) -> (face, position), as a sorted vector (a std::map per cell costs a few KB per
// atom); lookups by binary search, iteration in key order, and the last insertion wins for a repeated key, as in a map
struct DartMap {
    typedef std::pair<std::pair<int, int>, std::pair<int, int>> KV;
    std::vector<KV> v;
    typedef std::vector<KV>::const_iterator It;
    void clear() { v.clear(); }
    void put(std::pair<int, int> k, std::pair<int, int> val) { v.push_back({k, val}); }
    void finish() {
        std::stable_sort(v.begin(), v.end(), [](const KV& a, const KV& b) { return a.first < b.first; });
        std::vector<KV> w;
        for (size_t i = 0; i < v.size(); i++)
            if (i + 1 == v.size() || v[i + 1].first != v[i].first) w.push_back(v[i]);
        w.shrink_to_fit(); v.swap(w);
    }
    It find(const std::pair<int, int>& k) const {
        It it = std::lower_bound(v.begin(), v.end(), k, [](const KV& a, const std::pair<int, int>& b) { return a.first < b; });
        return (it != v.end() && it->first == k) ? it : v.end();
    }
    It begin() const { return v.begin(); }
    It end() const { return v.end(); }
};

// face cycles of a cell in one flat array (vertex ids) with offsets: a vector per face costs a heap block each
struct FaceView {
    const int* p; int n;
    size_t size() const { return (size_t)n; }
    int operator[](size_t k) const { return p[k]; }
    const int* begin() const { return p; }
    const int* end() const { return p + n; }
};
struct Cycles {
    std::vector<int> v, o{0};
    size_t size() const { return o.size() - 1; }
    FaceView operator[](size_t f) const { return FaceView{v.data() + o[f], o[f + 1] - o[f]}; }
    void clear() { v.clear(); o.assign(1, 0); }
    void add(const int* a, int n) { v.insert(v.end(), a, a + n); o.push_back((int)v.size()); }
    void set(const std::vector<std::vector<int>>& c) { clear(); for (auto& f : c) add(f.data(), (int)f.size()); }
    void shrink() { v.shrink_to_fit(); o.shrink_to_fit(); }
};

struct Cell {
    std::vector<int> nbr;                 // neighbour across face f
    Cycles cyc;                           // face cycles (local vertex ids)
    int nv = 0;
    bool simple = false, ok = false;
    std::vector<std::array<int, 4>> vkey; // global vertex identity (sorted cell ids), simple cells; freed once the
                                          // faces are aligned
    std::vector<int16_t> al;              // per face: position in the neighbour's copy of the face of this copy's first
                                          // vertex, + 128 if reversed (-1: none; >= 1000: an explicit map)
    // half-edge tables: dart = (face f, position k): u = cyc[f][k] -> v = cyc[f][k+1]
    DartMap he;                                               // (u,v) -> (f,k)
    uint64_t sig = 0;
};

static uint64_t sig_hash(const Cycles& cyc) {
    std::vector<int> s;
    for (size_t f = 0; f < cyc.size(); f++) s.push_back((int)cyc[f].size());
    std::sort(s.begin(), s.end());
    uint64_t h = 1469598103934665603ULL;
    for (int x : s) { h ^= (uint64_t)x; h *= 1099511628211ULL; }
    h ^= s.size(); h *= 1099511628211ULL;
    return h;
}

// half-edge table of a cell: built for reference cells once, and for a sample cell only while its labellings are
// found (then freed), so that the tables of all cells never exist at once
static void build_he(Cell& c) {
    c.he.clear();
    for (size_t f = 0; f < c.cyc.size(); f++)
        for (size_t k = 0; k < c.cyc[f].size(); k++)
            c.he.put({c.cyc[f][k], c.cyc[f][(k + 1) % c.cyc[f].size()]}, {(int)f, (int)k});
    c.he.finish();
}

static void free_he(Cell& c) { std::vector<DartMap::KV>().swap(c.he.v); }

static void build(Cell& c, int self, bool with_he) {
    if (with_he) build_he(c); else free_he(c);
    int maxv = -1;
    for (size_t f = 0; f < c.cyc.size(); f++)
        for (size_t k = 0; k < c.cyc[f].size(); k++) maxv = std::max(maxv, c.cyc[f][k]);
    c.nv = maxv + 1;
    // faces at each vertex (the first three, in face order) and their count
    std::vector<std::array<int, 3>> inc(c.nv); std::vector<int> ninc(c.nv, 0);
    for (size_t f = 0; f < c.cyc.size(); f++) for (int v : c.cyc[f]) { if (ninc[v] < 3) inc[v][ninc[v]] = (int)f; ninc[v]++; }
    c.simple = true;
    c.vkey.assign(c.nv, {-1, -1, -1, -1});
    for (int v = 0; v < c.nv; v++) {
        if (ninc[v] != 3) { c.simple = false; continue; }
        std::array<int, 4> k = {self, c.nbr[inc[v][0]], c.nbr[inc[v][1]], c.nbr[inc[v][2]]};
        std::sort(k.begin(), k.end());
        c.vkey[v] = k;
    }
    c.sig = sig_hash(c.cyc);
}

// ------------------------------------------------------------------ references
struct Ref {
    int species, mode;
    bool nbrglue = false;          // degenerate reference: glue by centers and common neighbors, not by vertices
    Cell cell, mirror;
    std::vector<uint64_t> flab_sig; std::vector<int> flab_sp;
    std::vector<V3> d, P;
};

// isomorphism cell -> target (orientation preserving), first half-edge of cell to target dart s
struct IsoFull { int orbit; std::vector<int16_t> vmap, fmap; };   // reference vertex and face ids (< 32767, checked)

static bool isomorphism(const Cell& a, const Cell& b, std::pair<int, int> H0, IsoFull& out,
                        const std::vector<uint64_t>* fl_sig, const std::vector<int>* fl_sp,
                        const Ref* r, bool use_species) {
    out.vmap.assign(a.nv, -1);
    out.fmap.assign(a.cyc.size(), -1);
    std::vector<int> finv(b.cyc.size(), -1);
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> st;
    std::pair<int, int> h0 = {a.cyc[0][0], a.cyc[0][1 % a.cyc[0].size()]};
    st.push_back({h0, H0});
    while (!st.empty()) {
        auto [h, H] = st.back(); st.pop_back();
        auto ia = a.he.find(h), ib = b.he.find(H);
        if (ia == a.he.end() || ib == b.he.end()) return false;
        int f = ia->second.first, F = ib->second.first;
        if (out.vmap[h.first] == -1) out.vmap[h.first] = H.first; else if (out.vmap[h.first] != H.first) return false;
        if (out.vmap[h.second] == -1) out.vmap[h.second] = H.second; else if (out.vmap[h.second] != H.second) return false;
        if (out.fmap[f] != -1) { if (out.fmap[f] != F) return false; continue; }
        if (finv[F] != -1) return false;
        if (a.cyc[f].size() != b.cyc[F].size()) return false;
        if (fl_sig) {
            if (r->flab_sig[F] != 0 && (*fl_sig)[f] != r->flab_sig[F]) return false;   // label 0: any neighbor
            if (use_species && (*fl_sp)[f] != r->flab_sp[F]) return false;
        }
        out.fmap[f] = F; finv[F] = f;
        int n = (int)a.cyc[f].size();
        int ka = ia->second.second, kb = ib->second.second;
        for (int t = 0; t < n; t++) {
            int u = a.cyc[f][(ka + t) % n], v = a.cyc[f][(ka + t + 1) % n];
            int U = b.cyc[F][(kb + t) % n], W = b.cyc[F][(kb + t + 1) % n];
            if (out.vmap[u] == -1) out.vmap[u] = U; else if (out.vmap[u] != U) return false;
            st.push_back({{v, u}, {W, U}});
        }
    }
    for (int x : out.fmap) if (x < 0) return false;
    return true;
}

// ------------------------------------------------------------------ small linear algebra
static bool inv3(const M3& m, M3& r) {
    double det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
    if (std::fabs(det) < 1e-12) return false;
    double id = 1.0 / det;
    r = {(m[4] * m[8] - m[5] * m[7]) * id, (m[2] * m[7] - m[1] * m[8]) * id, (m[1] * m[5] - m[2] * m[4]) * id,
         (m[5] * m[6] - m[3] * m[8]) * id, (m[0] * m[8] - m[2] * m[6]) * id, (m[2] * m[3] - m[0] * m[5]) * id,
         (m[3] * m[7] - m[4] * m[6]) * id, (m[1] * m[6] - m[0] * m[7]) * id, (m[0] * m[4] - m[1] * m[3]) * id};
    return true;
}
static M3 mul(const M3& a, const M3& b) {
    M3 c{};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) { double s = 0; for (int k = 0; k < 3; k++) s += a[3 * i + k] * b[3 * k + j]; c[3 * i + j] = s; }
    return c;
}
static V3 mvec(const M3& a, const V3& x) {
    return {a[0] * x[0] + a[1] * x[1] + a[2] * x[2], a[3] * x[0] + a[4] * x[1] + a[5] * x[2], a[6] * x[0] + a[7] * x[1] + a[8] * x[2]};
}

// ------------------------------------------------------------------ globals
static std::vector<Ref> refs;
static std::vector<M3> PG;
static double TOL;
static std::vector<Cell> cells;
static std::vector<std::vector<int>> special_maps;     // face alignments that are not a cyclic shift (none expected)
static uint64_t hmix(uint64_t h, uint64_t x) { h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); return h; }

// canonical form of a cell, independent of how voro++ numbered it (which can vary from run to run with threads): faces
// ordered by neighbour, vertices numbered in order of the sorted neighbours of their faces, each face's cycle started
// at its smallest vertex (orientation kept). Input: neighbours nb and voro++'s face list fv (per face: count, then
// vertices); output into nbr and cyc. Flat arrays throughout (a vector per vertex or face costs a heap block each).
static void canonical(const std::vector<int>& nb, const std::vector<int>& fv, std::vector<int>& nbr, Cycles& cyc) {
    int nf = (int)nb.size(), nv = 0;
    std::vector<int> fo_(nf + 1);                              // offsets of the faces' vertex lists in fv
    { size_t k = 0; for (int f = 0; f < nf; f++) { fo_[f] = (int)k + 1; for (int t = 0; t < fv[k]; t++) nv = std::max(nv, fv[k + 1 + t] + 1); k += fv[k] + 1; } fo_[nf] = (int)k + 1; }
    auto flen = [&](int f) { return fv[fo_[f] - 1]; };
    // a vertex's key: the sorted neighbours of its faces (flat, with offsets)
    std::vector<int> ko(nv + 1, 0), kv;
    for (int f = 0; f < nf; f++) for (int t = 0; t < flen(f); t++) ko[fv[fo_[f] + t] + 1]++;
    for (int v = 0; v < nv; v++) ko[v + 1] += ko[v];
    kv.resize(ko[nv]);
    { std::vector<int> pos(ko.begin(), ko.end() - 1);
      for (int f = 0; f < nf; f++) for (int t = 0; t < flen(f); t++) kv[pos[fv[fo_[f] + t]]++] = nb[f]; }
    for (int v = 0; v < nv; v++) std::sort(kv.begin() + ko[v], kv.begin() + ko[v + 1]);
    auto kless = [&](int a, int b) { return std::lexicographical_compare(kv.begin() + ko[a], kv.begin() + ko[a + 1], kv.begin() + ko[b], kv.begin() + ko[b + 1]); };
    std::vector<int> fo(nf); for (int f = 0; f < nf; f++) fo[f] = f;
    auto fkey = [&](int f) { std::vector<int> k; for (int t = 0; t < flen(f); t++) { int v = fv[fo_[f] + t]; k.insert(k.end(), kv.begin() + ko[v], kv.begin() + ko[v + 1]); } std::sort(k.begin(), k.end()); return k; };
    std::stable_sort(fo.begin(), fo.end(), [&](int a, int b) { return nb[a] != nb[b] ? nb[a] < nb[b] : fkey(a) < fkey(b); });
    std::vector<int> vo(nv); for (int v = 0; v < nv; v++) vo[v] = v;
    std::stable_sort(vo.begin(), vo.end(), kless);
    std::vector<int> newid(nv); for (int r = 0; r < nv; r++) newid[vo[r]] = r;
    nbr.resize(nf); cyc.clear(); cyc.v.reserve(fv.size() - nf); cyc.o.reserve(nf + 1);
    int c[64]; std::vector<int> big;
    for (int r = 0; r < nf; r++) {
        int f = fo[r], n = flen(f); nbr[r] = nb[f];
        int* p = c; if (n > 64) { big.resize(n); p = big.data(); }
        for (int t = 0; t < n; t++) p[t] = newid[fv[fo_[f] + t]];
        std::rotate(p, std::min_element(p, p + n), p + n);
        cyc.add(p, n);
    }
}
// a labelling of cell i, compactly: base maps (reference vertex of each vertex, then reference face of each face) at
// lpool[i][off...], composed with automorphism aut of the reference (autall[orbit][aut]; -1: none). The labellings of
// a cell that keeps all of them share one base map, so each costs 8 bytes, not a copy of the maps
struct Iso { int16_t orbit, aut; int32_t off; };
static std::vector<std::vector<Iso>> isos;
static std::vector<std::vector<int16_t>> lpool;
static std::vector<std::vector<IsoFull>> autall;
static std::vector<V3> X;
static double Lb[3]; static int per[3];
// triclinic box (fully periodic): a = (Lb0,0,0), b = (xy,Lb1,0), c = (xz,yz,Lb2), LAMMPS restricted-triclinic form
static bool tri = false; static double txy = 0, txz = 0, tyz = 0;

static int nearest_pg(const M3& R) {
    int best = -1; double bs = -1e300;
    for (size_t g = 0; g < PG.size(); g++) {
        double s = 0; for (int k = 0; k < 9; k++) s += PG[g][k] * R[k];
        // tr(P^T R) = 3 at a match; two distinct crystallographic elements have tr(P^T Q) <= 2, so a score
        // above 2.5 cannot be beaten by another element close enough to pass the gluing tolerance
        if (s > 2.5) return (int)g;
        if (s > bs) { bs = s; best = (int)g; }
    }
    return best;
}

// frame data of cell i under labelling a
static inline int vm(int i, const Iso& m, int v) { int x = lpool[i][m.off + v]; return m.aut < 0 ? x : autall[m.orbit][m.aut].vmap[x]; }
static inline int fm(int i, const Iso& m, int f) { int x = lpool[i][m.off + cells[i].nv + f]; return m.aut < 0 ? x : autall[m.orbit][m.aut].fmap[x]; }
static inline const V3& Pof(int i, const Iso& m, int v) { return refs[m.orbit].P[vm(i, m, v)]; }
static inline const V3& Dof(int i, const Iso& m, int f) { return refs[m.orbit].d[fm(i, m, f)]; }

// two labellings of one cell are equivalent if a crystal symmetry maps one onto the other (the same reading in another
// frame): choosing between them is a choice of gauge and cannot change any result. Every frame relation and every fit
// is then the same for both, up to that symmetry: if labelling b is R applied to labelling a, every ideal position b
// gives is R times the one a gives, so a frame relation g with a neighbor becomes R g, with the same misfit, and the
// misfit of a cell (after the best orthogonal map) is the same. Labellings are therefore tested one class at a time.
static bool equiv_lab(int j, const Iso& a, const Iso& b) {
    if (a.orbit != b.orbit) return false;
    size_t nf = cells[j].nbr.size();
    for (const M3& R : PG) {
        bool same = true;
        for (size_t f = 0; f < nf && same; f++) {
            V3 x = mvec(R, Dof(j, a, (int)f)); const V3& y = Dof(j, b, (int)f);
            double e = 0, s2 = 0; for (int k = 0; k < 3; k++) { e += (x[k] - y[k]) * (x[k] - y[k]); s2 += y[k] * y[k]; }
            same = e <= 1e-12 * (1 + s2);
        }
        if (same) return true;
    }
    return false;
}

// gluing of (i with labelling ma) and (j with labelling mb) across their shared face, fitted from i's side: R index
// or -1
// the gluing from the two centers and the common neighbors alone (degenerate references, whose vertex positions are
// not fixed by the labeling): the unique point-group element g with x_i = g x_j + d on all of them, or -1
static int glue_nbr(int i, const Iso& ma, int j, const Iso& mb, int fi, int fj, V3* dout) {
    const Cell& ci = cells[i]; const Cell& cj = cells[j];
    V3 d = Dof(i, ma, fi), dj = Dof(j, mb, fj);
    std::vector<std::pair<V3, V3>> cor;                  // (position in i's frame, in j's frame)
    cor.push_back({V3{0, 0, 0}, dj});
    for (size_t fa = 0; fa < ci.nbr.size(); fa++) {
        int kk = ci.nbr[fa];
        if (kk == j) continue;
        for (size_t fb = 0; fb < cj.nbr.size(); fb++) if (cj.nbr[fb] == kk) { cor.push_back({Dof(i, ma, (int)fa), Dof(j, mb, (int)fb)}); break; }
    }
    int found = -1;
    for (int g = 0; g < (int)PG.size(); g++) {
        double e = 0;
        for (size_t k = 0; k < cor.size() && e <= TOL; k++) {
            V3 y = mvec(PG[g], cor[k].second);
            for (int r = 0; r < 3; r++) e = std::max(e, std::fabs(y[r] + d[r] - cor[k].first[r]));
        }
        if (e <= TOL) { if (found >= 0) return -1; found = g; }
    }
    if (found >= 0 && dout) *dout = d;
    return found;
}

static int glue_one(int i, const Iso& ma, int j, const Iso& mb, V3* dout) {
    const Cell& ci = cells[i]; const Cell& cj = cells[j];
    int fi = -1, fj = -1;
    for (size_t f = 0; f < ci.nbr.size(); f++) if (ci.nbr[f] == j) { fi = (int)f; break; }
    for (size_t f = 0; f < cj.nbr.size(); f++) if (cj.nbr[f] == i) { fj = (int)f; break; }
    if (fi < 0 || fj < 0) return -1;
    if (refs[ma.orbit].nbrglue || refs[mb.orbit].nbrglue) return glue_nbr(i, ma, j, mb, fi, fj, dout);
    const auto& vi = ci.cyc[fi]; const auto& vj = cj.cyc[fj];
    int n = (int)vi.size();
    if ((int)vj.size() != n) return -1;
    if (n > 64) return -1;               // no Voronoi face of a crystalline cell comes close
    int code = ci.al.empty() ? -1 : ci.al[fi];
    if (code < 0) return -1;
    V3 Pi[64], Pj[64];
    for (int k = 0; k < n; k++) {
        int q = code >= 1000 ? special_maps[code - 1000][k] : (code >= 128 ? ((code - 128 - k) % n + n) % n : (code + k) % n);
        Pi[k] = Pof(i, ma, vi[k]); Pj[k] = Pof(j, mb, vj[q]);
    }
    V3 d = Dof(i, ma, fi);
    M3 A{}, B{}, Bi;
    for (int k = 0; k < 3; k++) for (int r = 0; r < 3; r++) { A[3 * r + k] = Pi[k][r] - d[r]; B[3 * r + k] = Pj[k][r]; }
    if (!inv3(B, Bi)) return -1;
    int g = nearest_pg(mul(A, Bi));
    double err = 0;
    for (int k = 0; k < n; k++) { V3 y = mvec(PG[g], Pj[k]); for (int r = 0; r < 3; r++) err = std::max(err, std::fabs(y[r] + d[r] - Pi[k][r])); }
    if (err > TOL) return -1;
    if (dout) *dout = d;
    return g;
}

// the gluing, accepted only if the fit passes from both sides: near the tolerance the two fits can disagree, and the
// labelling check (one order) and the gluing output (the other order) must give the same answer
static std::atomic<long long> n_glue{0};
static int glue(int i, const Iso& ma, int j, const Iso& mb, V3* dout) {
    n_glue.fetch_add(1, std::memory_order_relaxed);
    int g = glue_one(i, ma, j, mb, dout);
    if (g < 0 || glue_one(j, mb, i, ma, nullptr) < 0) return -1;
    return g;
}

static V3 minimage(V3 x) {
    if (tri) {                                        // reduce along c, then b, then a (as LAMMPS does)
        double m = std::round(x[2] / Lb[2]); x[0] -= m * txz; x[1] -= m * tyz; x[2] -= m * Lb[2];
        m = std::round(x[1] / Lb[1]); x[0] -= m * txy; x[1] -= m * Lb[1];
        m = std::round(x[0] / Lb[0]); x[0] -= m * Lb[0];
        // in a strongly tilted box the result can be one box vector away from the nearest image (a vector more than
        // half the box height along c is reduced into a long one); the 26 neighboring images are checked, and the
        // result is unchanged wherever the reduction above is already the nearest image
        V3 best = x; double bn = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
        for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int k = -1; k <= 1; k++) {
            if (!i && !j && !k) continue;
            V3 y = {x[0] + i * Lb[0] + j * txy + k * txz, x[1] + j * Lb[1] + k * tyz, x[2] + k * Lb[2]};
            double yn = y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
            if (yn < bn * (1 - 1e-12)) { best = y; bn = yn; }
        }
        return best;
    }
    for (int k = 0; k < 3; k++) if (per[k]) x[k] -= Lb[k] * std::round(x[k] / Lb[k]);
    return x;
}

// mean-square misfit of cell i's actual neighbour vectors vs. the ideal ones of labelling m
static double resid(int i, const Iso& m) {
    const Cell& c = cells[i];
    int nf = (int)c.nbr.size();
    double H[9] = {0}, AA = 0;
    std::vector<V3> A(nf), B(nf);
    for (int f = 0; f < nf; f++) {
        A[f] = Dof(i, m, f);
        V3 x = {X[c.nbr[f]][0] - X[i][0], X[c.nbr[f]][1] - X[i][1], X[c.nbr[f]][2] - X[i][2]};
        B[f] = minimage(x);
    }
    // best orthogonal map by a simple polar-decomposition iteration (Higham), then residual
    M3 Mm{};
    for (int f = 0; f < nf; f++) for (int r = 0; r < 3; r++) for (int s = 0; s < 3; s++) Mm[3 * r + s] += B[f][r] * A[f][s];
    M3 Q = Mm;
    for (int it = 0; it < 30; it++) {
        M3 Qi; if (!inv3(Q, Qi)) break;
        M3 Qn; for (int r = 0; r < 3; r++) for (int s = 0; s < 3; s++) Qn[3 * r + s] = 0.5 * (Q[3 * r + s] + Qi[3 * s + r]);
        Q = Qn;
    }
    double num = 0;
    for (int f = 0; f < nf; f++) { V3 y = mvec(Q, A[f]); for (int r = 0; r < 3; r++) num += B[f][r] * y[r]; AA += A[f][0] * A[f][0] + A[f][1] * A[f][1] + A[f][2] * A[f][2]; }
    double s = num / AA, e = 0;
    for (int f = 0; f < nf; f++) { V3 y = mvec(Q, A[f]); for (int r = 0; r < 3; r++) { double q = B[f][r] - s * y[r]; e += q * q; } }
    (void)H;
    return e / nf;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: vb_engine refs.txt sample.txt\n"); return 1; }
    double T0 = omp_get_wtime(), Tl = T0;
    auto tick = [&](const char* w) { if (getenv("VB_TIMING")) { double t = omp_get_wtime(); fprintf(stderr, "%-18s %.3f s  (glue calls so far %lld)\n", w, t - Tl, (long long)n_glue.load()); Tl = t; } };
    // ---------------- references
    FILE* fr = fopen(argv[1], "r");
    int nref, G;
    if (fscanf(fr, "%d %d %lf", &nref, &G, &TOL) != 3) return 2;
    PG.resize(G);
    for (int g = 0; g < G; g++) for (int k = 0; k < 9; k++) if (fscanf(fr, "%lf", &PG[g][k]) != 1) return 2;
    refs.resize(nref);
    for (int r = 0; r < nref; r++) {
        char tag[8]; int nf, nv;
        if (fscanf(fr, "%7s %d %d %d %d", tag, &refs[r].species, &nf, &nv, &refs[r].mode) != 5) return 2;
        if (nf >= 32767 || nv >= 32767) { fprintf(stderr, "reference cell too large (labelling maps are 16-bit)\n"); return 2; }
        refs[r].nbrglue = refs[r].mode & 4; refs[r].mode &= 3;
        Cell& c = refs[r].cell;
        c.nbr.assign(nf, -1);
        std::vector<std::vector<int>> rc(nf);
        refs[r].flab_sig.resize(nf); refs[r].flab_sp.resize(nf);
        for (int f = 0; f < nf; f++) {
            unsigned long long h; int sp, n;
            if (fscanf(fr, "%llu %d %d", &h, &sp, &n) != 3) return 2;
            refs[r].flab_sig[f] = h; refs[r].flab_sp[f] = sp;
            rc[f].resize(n);
            for (int k = 0; k < n; k++) if (fscanf(fr, "%d", &rc[f][k]) != 1) return 2;
            c.nbr[f] = -(f + 10);         // dummy distinct ids
        }
        refs[r].d.resize(nf); refs[r].P.resize(nv);
        for (int f = 0; f < nf; f++) if (fscanf(fr, "%lf %lf %lf", &refs[r].d[f][0], &refs[r].d[f][1], &refs[r].d[f][2]) != 3) return 2;
        for (int v = 0; v < nv; v++) if (fscanf(fr, "%lf %lf %lf", &refs[r].P[v][0], &refs[r].P[v][1], &refs[r].P[v][2]) != 3) return 2;
        c.cyc.set(rc);
        build(c, -1, true);
        Cell& m = refs[r].mirror;
        m.nbr = c.nbr;
        for (auto& cy : rc) std::reverse(cy.begin(), cy.end());
        m.cyc.set(rc);
        build(m, -1, true);
    }
    fclose(fr);

    tick("read_refs");
    // ---------------- sample and tessellation
    // text: "N Lx Ly Lz px py pz use_sp" then N lines "x y z [species]";
    // binary (file name ending in .bin): int32 N, 3 doubles L, 4 int32 (px py pz use_sp), 3N doubles, [N int32]
    int N, use_sp;
    std::vector<int> spc;
    std::string sname(argv[2]);
    if (sname.size() > 4 && sname.compare(sname.size() - 4, 4, ".bin") == 0) {
        FILE* fs = fopen(argv[2], "rb");
        int32_t n32, flags[4];
        if (fread(&n32, 4, 1, fs) != 1 || fread(Lb, 8, 3, fs) != 3 || fread(flags, 4, 4, fs) != 4) return 3;
        N = n32; for (int k = 0; k < 3; k++) per[k] = flags[k]; use_sp = flags[3];
        if (flags[0] == 2) {                          // triclinic, fully periodic: three tilts follow
            double t[3];
            if (fread(t, 8, 3, fs) != 3) return 3;
            tri = true; txy = t[0]; txz = t[1]; tyz = t[2]; per[0] = per[1] = per[2] = 1;
        }
        X.resize(N); spc.assign(N, 0);
        if (fread(X.data(), 8, 3 * (size_t)N, fs) != 3 * (size_t)N) return 3;
        if (use_sp) { std::vector<int32_t> t(N); if (fread(t.data(), 4, N, fs) != (size_t)N) return 3; for (int i = 0; i < N; i++) spc[i] = t[i]; }
        fclose(fs);
    } else {
        FILE* fs = fopen(argv[2], "r");
        if (fscanf(fs, "%d %lf %lf %lf %d %d %d %d", &N, &Lb[0], &Lb[1], &Lb[2], &per[0], &per[1], &per[2], &use_sp) != 8) return 3;
        X.resize(N); spc.assign(N, 0);
        for (int i = 0; i < N; i++) {
            if (fscanf(fs, "%lf %lf %lf", &X[i][0], &X[i][1], &X[i][2]) != 3) return 3;
            if (use_sp && fscanf(fs, "%d", &spc[i]) != 1) return 3;
        }
        fclose(fs);
    }
    int threads = omp_get_max_threads();
    double dens = N / (Lb[0] * Lb[1] * Lb[2]);
    int nx = std::max(1, int(Lb[0] * pow(dens / 5, 1. / 3))), ny = std::max(1, int(Lb[1] * pow(dens / 5, 1. / 3))),
        nz = std::max(1, int(Lb[2] * pow(dens / 5, 1. / 3)));
    cells.resize(N);
    auto tessellate = [&](auto& con) {
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
                    Cell& c = cells[pid];
                    canonical(nb, fv, c.nbr, c.cyc);
                    c.nbr.shrink_to_fit();
                    c.ok = true;
                }
        }
    };
    if (tri) { container_triclinic con(Lb[0], txy, Lb[1], txz, tyz, Lb[2], nx, ny, nz, 8, threads); tessellate(con); }
    else { container_3d con(0, Lb[0], 0, Lb[1], 0, Lb[2], nx, ny, nz, per[0], per[1], per[2], 8, threads); tessellate(con); }
    tick("tessellate");
    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < N; i++) if (cells[i].ok) build(cells[i], i, false);
    // align every shared face once (vertex k of i's copy -> the vertex of j's copy with the same global identity, as
    // the gluings did), then free the vertex identities, which nothing else needs
    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < N; i++) {
        Cell& c = cells[i];
        if (!c.ok || !c.simple) continue;
        c.al.assign(c.nbr.size(), -1);
        for (size_t f = 0; f < c.nbr.size(); f++) {
            int j = c.nbr[f];
            if (j < 0 || !cells[j].ok) continue;
            bool firstf = true; for (size_t g = 0; g < f; g++) if (c.nbr[g] == j) firstf = false;
            if (!firstf) continue;                                  // gluings use the first face with j
            const Cell& d = cells[j];
            int fj = -1; for (size_t g = 0; g < d.nbr.size(); g++) if (d.nbr[g] == i) { fj = (int)g; break; }
            if (fj < 0) continue;
            const auto& vi = c.cyc[f]; const auto& vj = d.cyc[fj];
            int n = (int)vi.size();
            if ((int)vj.size() != n || n > 64) continue;
            std::vector<int> q(n, -1); bool all = true;
            for (int k = 0; k < n && all; k++) {
                for (int x = 0; x < n; x++) if (d.vkey[vj[x]] == c.vkey[vi[k]]) { q[k] = x; break; }
                if (q[k] < 0) all = false;
            }
            if (!all) continue;
            int o = q[0]; bool fwd = true, bwd = true;
            for (int k = 0; k < n; k++) { if (q[k] != (o + k) % n) fwd = false; if (q[k] != ((o - k) % n + n) % n) bwd = false; }
            if (fwd) c.al[f] = (int16_t)o;
            else if (bwd) c.al[f] = (int16_t)(o + 128);
            else {
                #pragma omp critical
                { c.al[f] = (int16_t)(1000 + (int)special_maps.size()); special_maps.push_back(q); }
            }
        }
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) std::vector<std::array<int, 4>>().swap(cells[i].vkey);
    tick("align");

    tick("build");
    // ---------------- candidate labellings
    isos.assign(N, {});
    std::map<uint64_t, std::vector<int>> by_sig;
    for (int r = 0; r < nref; r++) by_sig[refs[r].cell.sig].push_back(r);
    // all labellings of a cell (mode 2): two isomorphisms onto a reference differ by an automorphism of the reference,
    // and an orientation-preserving automorphism of a polyhedron is fixed by the image of one dart. So once the first
    // labelling phi is found (first dart D1 in table order), the one starting at a later dart D is beta phi, where
    // beta is the automorphism taking D1 to D, if there is one (no isomorphism starts at D otherwise); it is kept if
    // its face labels match, which is the check the search itself makes. auts[r][mir]: the automorphisms of the
    // reference cell (mir 0) and of its mirror image (mir 1)
    std::vector<std::array<std::vector<IsoFull>, 2>> auts(nref);
    autall.assign(nref, {});
    for (int r = 0; r < nref; r++) {
        if (refs[r].mode != 2) continue;
        for (int mir = 0; mir < 2; mir++) {
            const Cell& T = mir ? refs[r].mirror : refs[r].cell;
            for (auto& kv : T.he) { IsoFull a; if (isomorphism(T, T, kv.first, a, nullptr, nullptr, nullptr, false)) auts[r][mir].push_back(a); }
            autall[r].insert(autall[r].end(), auts[r][mir].begin(), auts[r][mir].end());
        }
        if (autall[r].size() >= 32767) { fprintf(stderr, "reference with too many automorphisms\n"); return 2; }
    }
    lpool.assign(N, {});
    auto add_base = [&](int i, const IsoFull& m) {             // the maps of m into cell i's pool: their offset
        int off = (int)lpool[i].size();
        lpool[i].insert(lpool[i].end(), m.vmap.begin(), m.vmap.end()); lpool[i].insert(lpool[i].end(), m.fmap.begin(), m.fmap.end());
        return off;
    };
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < N; i++) {
        Cell& c = cells[i];
        if (!c.ok || !c.simple) continue;
        bool wall = false;
        for (int nbv : c.nbr) if (nbv < 0) wall = true;
        if (wall) continue;
        auto it = by_sig.find(c.sig);
        if (it == by_sig.end()) continue;
        std::vector<uint64_t> fsig(c.nbr.size()); std::vector<int> fsp(c.nbr.size());
        for (size_t f = 0; f < c.nbr.size(); f++) {
            const Cell& nc = cells[c.nbr[f]];
            fsig[f] = nc.ok ? nc.sig : 0;
            fsp[f] = spc[c.nbr[f]];
        }
        build_he(c);
        for (int r : it->second) {
            const Ref& R = refs[r];
            if (use_sp && spc[i] != R.species) continue;
            bool labels = R.mode != 0;
            for (int mir = 0; mir < 2; mir++) {
                const Cell& T = mir ? R.mirror : R.cell;
                bool found = false;
                if (R.mode == 2) {
                    const auto& A = auts[r][mir];
                    int a0 = mir ? (int)auts[r][0].size() : 0;    // A[g] is autall[r][a0 + g]
                    std::vector<int> beta;                    // per dart of T (table order): the automorphism D1 -> it
                    IsoFull phi; int off = 0;
                    for (size_t t = 0; t < T.he.v.size(); t++) {
                        const auto& D = T.he.v[t].first;
                        if (!found) {
                            phi.orbit = r;
                            if (!isomorphism(c, T, D, phi, labels ? &fsig : nullptr, &fsp, &R, use_sp)) continue;
                            off = add_base(i, phi); isos[i].push_back(Iso{(int16_t)r, -1, off}); found = true;
                            beta.assign(T.he.v.size(), -1);
                            for (size_t g = 0; g < A.size(); g++) {
                                auto e = T.he.find({A[g].vmap[D.first], A[g].vmap[D.second]});
                                if (e != T.he.end() && beta[e - T.he.begin()] < 0) beta[e - T.he.begin()] = (int)g;
                            }
                            continue;
                        }
                        if (beta[t] < 0) continue;
                        const IsoFull& b = A[beta[t]];
                        bool okl = true;
                        for (size_t f = 0; f < phi.fmap.size() && okl; f++) {
                            int F = b.fmap[phi.fmap[f]];
                            if (labels && R.flab_sig[F] != 0 && fsig[f] != R.flab_sig[F]) okl = false;   // label 0: any neighbor
                            if (labels && use_sp && fsp[f] != R.flab_sp[F]) okl = false;
                        }
                        if (okl) isos[i].push_back(Iso{(int16_t)r, (int16_t)(a0 + beta[t]), off});
                    }
                } else
                    for (auto& kv : T.he) {
                        IsoFull m; m.orbit = r;
                        if (isomorphism(c, T, kv.first, m, labels ? &fsig : nullptr, &fsp, &R, use_sp)) {
                            isos[i].push_back(Iso{(int16_t)r, -1, add_base(i, m)}); found = true;
                            break;
                        }
                    }
                if (found && R.mode != 2) break;
            }
        }
        free_he(c);
    }

    tick("labellings");
    // ---------------- classes of equivalent labellings of each cell (a labelling's class: the first equivalent one)
    std::vector<std::vector<int>> lcls(N);
    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < N; i++) {
        size_t m = isos[i].size();
        lcls[i].assign(m, -1);
        std::vector<int> reps;
        for (size_t a = 0; a < m; a++) {
            for (int r : reps) if (equiv_lab(i, isos[i][r], isos[i][a])) { lcls[i][a] = r; break; }
            if (lcls[i][a] < 0) { lcls[i][a] = (int)a; reps.push_back((int)a); }
        }
    }
    tick("classes");
    // ---------------- whether two labelling classes of neighbouring cells glue, for every face of a cell with several
    // labellings: support is the same for a whole class (see equiv_lab), so it is computed once per pair of classes,
    // from their first labellings, in parallel. Arc consistency and region growth below ask exactly these questions
    // (with the same labellings) and read the answers here. A table entry holds
    // bit (ka * nc(j) + kb) for class ka of i and kb of j (nc(i) * nc(j) <= 64); sfill marks the entries filled.
    std::vector<int8_t> ncl(N, 0);
    std::vector<std::vector<int8_t>> cidx(N);
    for (int i = 0; i < N; i++) {
        size_t m = isos[i].size();
        if (m < 2) { ncl[i] = (int8_t)m; continue; }
        cidx[i].assign(m, -1); int k = 0;
        for (size_t a = 0; a < m; a++) if (lcls[i][a] == (int)a && k < 127) cidx[i][a] = (int8_t)k++;
        for (size_t a = 0; a < m; a++) cidx[i][a] = cidx[i][lcls[i][a]];
        ncl[i] = (int8_t)k;
    }
    std::vector<int64_t> soff(N + 1, 0);
    for (int i = 0; i < N; i++) soff[i + 1] = soff[i] + (isos[i].size() > 1 ? (int64_t)cells[i].nbr.size() : 0);
    std::vector<uint64_t> sup(soff[N], 0);
    std::vector<char> sfill(soff[N], 0);
    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < N; i++) {
        // a cell with more than 64 classes has no table entries (nc(i) * nc(j) > 64 for every neighbor j), and the
        // class index of a cell with more than 127 classes is -1 beyond the 127th (int8)
        if (isos[i].size() < 2 || ncl[i] > 64) continue;
        std::vector<int> ri(ncl[i]);
        for (size_t a = 0; a < isos[i].size(); a++) if (lcls[i][a] == (int)a) ri[cidx[i][a]] = (int)a;
        const auto& nb = cells[i].nbr;
        for (size_t f = 0; f < nb.size(); f++) {
            int j = nb[f];
            if (j < 0 || isos[j].empty() || ncl[i] * ncl[j] > 64) continue;
            bool firstf = true; for (size_t g = 0; g < f; g++) if (nb[g] == j) firstf = false;
            if (!firstf) continue;                              // gluings use the first face with j
            std::vector<int> rj(ncl[j]);
            for (size_t b = 0; b < isos[j].size(); b++) if (lcls[j][b] == (int)b) rj[ncl[j] > 1 ? cidx[j][b] : 0] = (int)b;
            uint64_t m = 0;
            for (int ka = 0; ka < ncl[i]; ka++) for (int kb = 0; kb < ncl[j]; kb++)
                if (glue(i, isos[i][ri[ka]], j, isos[j][rj[kb]], nullptr) >= 0) m |= 1ULL << (ka * ncl[j] + kb);
            sup[soff[i] + f] = m; sfill[soff[i] + f] = 1;
        }
    }
    auto kof = [&](int i, int a) { return isos[i].size() > 1 ? (int)cidx[i][a] : 0; };
    auto entry = [&](int i, int j) -> int64_t {              // table slot of the first face of i with j, -1 if none
        if (isos[i].size() < 2) return -1;
        const auto& nb = cells[i].nbr;
        for (size_t f = 0; f < nb.size(); f++) if (nb[f] == j) return sfill[soff[i] + f] ? soff[i] + (int64_t)f : -1;
        return -1;
    };
    // glue(i, a, j, b) >= 0, read from the table when a and b are the first labellings of their classes (the
    // labellings it was computed with), computed otherwise
    auto glues = [&](int i, int a, int j, int b) {
        if (lcls[i][a] == a && lcls[j][b] == b) {
            int64_t e = entry(i, j);
            if (e >= 0) return (bool)((sup[e] >> (kof(i, a) * ncl[j] + kof(j, b))) & 1);
            e = entry(j, i);
            if (e >= 0) return (bool)((sup[e] >> (kof(j, b) * ncl[i] + kof(i, a))) & 1);
        }
        return glue(i, isos[i][a], j, isos[j][b], nullptr) >= 0;
    };
    tick("class_gluings");
    if (getenv("VB_HASH")) {
        uint64_t h0 = 0, h1 = 0;
        for (int i = 0; i < N; i++) { h0 = hmix(h0, cells[i].nbr.size()); for (int x : cells[i].nbr) h0 = hmix(h0, x); for (int x : cells[i].cyc.v) h0 = hmix(h0, x);
            h1 = hmix(h1, isos[i].size()); for (auto& m : isos[i]) { h1 = hmix(h1, m.orbit); for (size_t f = 0; f < cells[i].nbr.size(); f++) h1 = hmix(h1, fm(i, m, (int)f)); } }
        fprintf(stderr, "HASH cells %016llx isos %016llx\n", (unsigned long long)h0, (unsigned long long)h1);
    }
    // ---------------- arc consistency, as in the close-packed route: a labelling of a cell is removed if some
    // neighbour has no remaining labelling that agrees with it, until nothing changes. It runs in rounds, each taking
    // all its decisions from the domains at its start, so the result depends on no order (and the cells of a round
    // are independent of one another). Every cell is revised, so of two neighbours that cannot both be right neither
    // is kept. A cell left with none no longer constrains its neighbours in that pass; once the pass ends, such cells
    // are defects and arc consistency starts again without them, until no cell is left with none. Support is the same
    // for a whole class (see equiv_lab), so it is decided once per class, from the first labelling of each (the
    // table's questions)
    // gluings between cells with a single labelling do not depend on the gauge: computed once, in parallel, for arc
    // consistency, region growth, and the output (cg = -2: not cached); cells that arc consistency leaves with a single
    // labelling are added after it
    std::vector<std::vector<int8_t>> cg(N);              // point-group index (< 128), -1 no gluing, -2 not cached
    auto cache_gluings = [&](auto single, auto lab) {
        #pragma omp parallel for schedule(dynamic, 256)
        for (int i = 0; i < N; i++) {
            if (!single(i)) continue;
            size_t nf = cells[i].nbr.size();
            if (cg[i].empty()) cg[i].assign(nf, -2);
            for (size_t f = 0; f < nf; f++) {
                int j = cells[i].nbr[f];
                if (j < 0 || cg[i][f] != -2 || !single(j)) continue;
                cg[i][f] = (int8_t)glue(i, isos[i][lab(i)], j, isos[j][lab(j)], nullptr);
            }
        }
    };
    cache_gluings([&](int i) { return isos[i].size() == 1; }, [&](int) { return 0; });
    std::vector<std::vector<int>> dom(N);
    std::vector<char> cand(N, 0), alive(N, 0);
    auto glue_ok = [&](int i, int ai, int j, int aj) {       // i, j neighbours; ai, aj labelling indices
        if (!cg[i].empty() && dom[j].size() == 1) {
            const auto& nb = cells[i].nbr;
            for (size_t f = 0; f < nb.size(); f++) if (nb[f] == j && cg[i][f] != -2) return cg[i][f] >= 0;
        }
        return glues(i, ai, j, aj);
    };
    for (int i = 0; i < N; i++) alive[i] = !isos[i].empty();
    std::vector<char> mark(N, 0);
    while (true) {
        std::vector<int> act, next;
        for (int i = 0; i < N; i++) {
            cand[i] = alive[i]; dom[i].clear();
            if (alive[i]) for (size_t a = 0; a < isos[i].size(); a++) dom[i].push_back((int)a);
            if (alive[i]) act.push_back(i);
        }
        while (!act.empty()) {
            std::vector<std::vector<int>> nd(act.size());
            std::vector<char> changed(act.size(), 0);
            #pragma omp parallel
            {
                std::vector<std::pair<int, char>> done;      // class of i -> supported by every neighbour
                std::vector<int> keep;
                #pragma omp for schedule(dynamic, 64)
                for (size_t t = 0; t < act.size(); t++) {
                    int i = act[t];
                    const auto& nb = cells[i].nbr;
                    done.clear(); keep.clear();
                    for (int a : dom[i]) {
                        int ca = lcls[i][a]; int hit = -1;
                        for (auto& d : done) if (d.first == ca) { hit = d.second; break; }
                        if (hit < 0) {
                            bool all = true;
                            for (size_t f = 0; f < nb.size() && all; f++) {
                                int j = nb[f];
                                if (j < 0 || !cand[j]) continue;
                                bool firstf = true; for (size_t g = 0; g < f; g++) if (nb[g] == j) firstf = false;
                                if (!firstf) continue;              // gluings use the first face with j
                                bool any = false;
                                for (int b : dom[j]) if (lcls[j][b] == b && glue_ok(i, ca, j, b)) { any = true; break; }
                                all = any;
                            }
                            hit = all; done.push_back({ca, (char)all});
                        }
                        if (hit) keep.push_back(a);
                    }
                    if (keep.size() != dom[i].size()) { changed[t] = 1; nd[t] = keep; }
                }
            }
            next.clear();
            for (size_t t = 0; t < act.size(); t++) if (changed[t]) {
                int i = act[t];
                dom[i].swap(nd[t]);
                if (dom[i].empty()) cand[i] = 0;
            }
            for (size_t t = 0; t < act.size(); t++) if (changed[t])
                for (int k : cells[act[t]].nbr) if (k >= 0 && cand[k] && !mark[k]) { mark[k] = 1; next.push_back(k); }
            for (int k : next) mark[k] = 0;
            act.swap(next);
        }
        bool any = false;
        for (int i = 0; i < N; i++) if (alive[i] && !cand[i]) { alive[i] = 0; any = true; }
        if (!any) break;
    }

    if (getenv("VB_TIMING")) { std::map<size_t, long> h; for (int i = 0; i < N; i++) if (!isos[i].empty()) h[isos[i].size()]++; fprintf(stderr, "labellings per cell:"); for (auto& kv : h) fprintf(stderr, " %zu:%ld", kv.first, kv.second); fprintf(stderr, "\n"); }
    tick("arc_consistency");
    if (getenv("VB_DEBUG")) { size_t multi = 0, emptied = 0; for (int i = 0; i < N; i++) { if (!isos[i].empty() && !cand[i]) emptied++; if (cand[i] && dom[i].size() > 1) multi++; } fprintf(stderr, "after arc consistency: %zu cells with several labellings, %zu emptied\n", multi, emptied); }
    // ---------------- gauge fixing by region growth (with the pseudo-symmetry tie-break)
    std::vector<int> chosen(N, -1), region(N, -1);
    std::vector<char> kept(N, 0);
    std::vector<int> order;
    for (int i = 0; i < N; i++) if (cand[i]) order.push_back(i);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return dom[a].size() < dom[b].size(); });
    std::vector<int> nfixed(N, 0);
    int nreg = 0, npseudo = 0;
    cache_gluings([&](int i) { return cand[i] && dom[i].size() == 1; }, [&](int i) { return dom[i][0]; });
    // two labellings of one cell are equivalent if a crystal symmetry maps one onto the other (the same reading in
    // another frame): choosing between them is a choice of gauge and cannot change any result
    auto equiv = [&](int j, const Iso& a, const Iso& b) { return equiv_lab(j, a, b); };
    int nforced = 0, ndeferred = 0;
    double min_margin = 1e300;           // smallest rms misfit of the runner-up minus that of the best (A), over regions with several readings
    // cells next to the growing region, by (number of fixed neighbours, index), largest first: an indexed binary heap
    // whose keys only increase (each cell at most once; the same order as a heap of (count, cell) entries that
    // skips outdated ones). hpos: a cell's place in the heap, -1 if absent (all -1 between growths)
    std::vector<int> hpos(N, -1), hp;
    // waiting cells by (fixed neighbours, -time it first waited), largest first: the cell a forced decision takes is
    // the first, in order of waiting, of those with the most fixed neighbours (an indexed heap; empty between growths)
    std::vector<int> wpos(N, -1), wp, wseq(N, 0);
    auto grow = [&](int seed, int a0, std::vector<int>& fixed, std::vector<char>& avail, std::vector<int>& nfx) {
        kept[seed] = 1; chosen[seed] = a0; avail[seed] = 0; fixed.push_back(seed);
        auto key = [&](int j) { return ((int64_t)nfx[j] << 32) | (uint32_t)j; };
        auto up = [&](int t) {
            int j = hp[t]; int64_t kj = key(j);
            while (t > 0) { int u = (t - 1) / 2; if (key(hp[u]) >= kj) break; hp[t] = hp[u]; hpos[hp[t]] = t; t = u; }
            hp[t] = j; hpos[j] = t;
        };
        auto bump = [&](int j) {                         // one more fixed neighbour of j
            nfx[j]++;
            if (hpos[j] < 0) { hp.push_back(j); hpos[j] = (int)hp.size() - 1; }
            up(hpos[j]);
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
        int nseq = 0;
        auto wkey = [&](int j) { return ((int64_t)nfx[j] << 32) | (uint32_t)(0x7fffffff - wseq[j]); };
        auto wup = [&](int t) {
            int j = wp[t]; int64_t kj = wkey(j);
            while (t > 0) { int u = (t - 1) / 2; if (wkey(wp[u]) >= kj) break; wp[t] = wp[u]; wpos[wp[t]] = t; t = u; }
            wp[t] = j; wpos[j] = t;
        };
        auto wpop = [&]() {
            int top = wp[0]; wpos[top] = -1;
            int last = wp.back(); wp.pop_back();
            if (!wp.empty()) {
                int t = 0, n = (int)wp.size(); int64_t kl = wkey(last);
                while (true) {
                    int c = 2 * t + 1; if (c >= n) break;
                    if (c + 1 < n && wkey(wp[c + 1]) > wkey(wp[c])) c++;
                    if (wkey(wp[c]) <= kl) break;
                    wp[t] = wp[c]; wpos[wp[t]] = t; t = c;
                }
                wp[t] = last; wpos[last] = t;
            }
            return top;
        };
        auto wait = [&](int j) {                         // j waits (again): it keeps the time it first waited
            if (wpos[j] >= 0) return;
            wseq[j] = nseq++; wp.push_back(j); wpos[j] = (int)wp.size() - 1; wup(wpos[j]);
        };
        auto bump2 = [&](int j) { bump(j); if (wpos[j] >= 0) wup(wpos[j]); };
        for (int j : cells[seed].nbr) if (j >= 0 && avail[j]) bump2(j);
        // the labelling of cell j: among those consistent with its fixed neighbours, the one that fits the actual
        // neighbour vectors best. Geometry only orders equivalent choices (a gauge); when inequivalent ones remain,
        // the cell waits for more fixed neighbours (force: all cells are waiting, and one must be decided)
        std::vector<int> ok;
        std::vector<std::pair<int, char>> done;              // consistency is the same for a whole class
        auto settle = [&](int j, bool force) {
            ok.clear(); done.clear();
            for (int a : dom[j]) {
                int ca = lcls[j][a]; int hit = -1;
                for (auto& d : done) if (d.first == ca) { hit = d.second; break; }
                if (hit < 0) {
                    bool okk = true;
                    const auto& nb = cells[j].nbr;
                    for (size_t f = 0; f < nb.size(); f++) {
                        int i = nb[f];
                        if (i < 0 || !kept[i]) continue;
                        // both single-labelled: the cached gluing, on j's side (a gluing passes from both sides or not
                        // at all, so it is the answer on i's side too)
                        bool g = (!cg[j].empty() && cg[j][f] != -2) ? cg[j][f] >= 0 : glue_ok(i, chosen[i], j, a);
                        if (!g) { okk = false; break; }
                    }
                    hit = okk; done.push_back({ca, (char)okk});
                }
                if (hit) ok.push_back(a);
            }
            if (ok.empty()) { avail[j] = 0; return; }                    // no consistent labelling: a defect
            bool one_class = true;
            for (size_t t = 1; t < ok.size(); t++) if (lcls[j][ok[t]] != lcls[j][ok[0]]) { one_class = false; break; }
            if (!force && !one_class) { wait(j); ndeferred++; return; }
            // equivalent labellings fit equally well (the misfit is taken after the best orthogonal map): take the
            // first; between inequivalent ones (a forced decision), the best fit
            int pick = ok[0]; double best = 1e300;
            if (!one_class) for (int a : ok) { double r = resid(j, isos[j][a]); if (r < best) { best = r; pick = a; } }
            avail[j] = 0; kept[j] = 1; chosen[j] = pick; fixed.push_back(j);
            for (int k : cells[j].nbr) if (k >= 0 && avail[k]) bump2(k);
        };
        while (true) {
            while (!hp.empty()) {
                int j = pop();
                if (!avail[j]) continue;
                settle(j, false);
            }
            // every reachable cell is settled or waiting: decide the waiting cell with the most fixed neighbours
            int w = -1;
            while (!wp.empty()) { int j = wpop(); if (avail[j] && !kept[j]) { w = j; break; } }   // others never return
            if (w < 0) break;
            nforced++;
            settle(w, true);
        }
    };
    std::vector<char> avail(cand.begin(), cand.end());
    for (int seed : order) {
        if (!avail[seed]) continue;
        // the seed's labellings up to the crystal's symmetries: each class is one reading of the region. Several
        // classes (a crystal whose topology is more symmetric than the crystal, e.g. an oblique lattice) are told
        // apart once, by growing the region from each and comparing the fit of all its cells (geometry decides only
        // here, with the evidence of the whole region); the margin to the runner-up is reported
        std::vector<int> cls;
        for (int a : dom[seed]) {
            bool seen = false;
            for (int c : cls) if (equiv(seed, isos[seed][c], isos[seed][a])) { seen = true; break; }
            if (!seen) cls.push_back(a);
        }
        int best = cls[0];
        std::vector<int> fixed;
        if (cls.size() > 1) {
            npseudo++;
            double bt = 1e300, bt2 = 1e300;
            // growth is deterministic from a given state, so the best trial's region is kept as it is (growing it
            // again from the same state would give the same cells and labellings)
            std::vector<int> bfixed, bchosen, bnf; std::vector<char> bav; int bforced = 0, bdeferred = 0;
            for (int c : cls) {
                std::vector<int> fx; std::vector<char> av = avail; std::vector<int> nf2 = nfixed;
                int f0 = nforced, d0 = ndeferred;
                grow(seed, c, fx, av, nf2);
                std::vector<double> r(fx.size());
                #pragma omp parallel for schedule(dynamic, 1024)
                for (size_t t = 0; t < fx.size(); t++) r[t] = resid(fx[t], isos[fx[t]][chosen[fx[t]]]);
                double tot = 0; for (double x : r) tot += x; tot /= fx.size();      // summed in order
                bool better = tot < bt;
                if (better) {
                    bfixed = fx; bav.swap(av); bnf.swap(nf2); bforced = nforced - f0; bdeferred = ndeferred - d0;
                    bchosen.resize(fx.size()); for (size_t t = 0; t < fx.size(); t++) bchosen[t] = chosen[fx[t]];
                }
                for (int i : fx) { kept[i] = 0; chosen[i] = -1; }
                if (better) { bt2 = bt; bt = tot; best = c; } else if (tot < bt2) bt2 = tot;
            }
            if (getenv("VB_DEBUG")) fprintf(stderr, "pseudo seed %d: %zu readings, best fit %.6g, next %.6g\n", seed, cls.size(), bt, bt2);
            min_margin = std::min(min_margin, std::sqrt(bt2) - std::sqrt(bt));
            fixed.swap(bfixed); avail.swap(bav); nfixed.swap(bnf);
            for (size_t t = 0; t < fixed.size(); t++) { kept[fixed[t]] = 1; chosen[fixed[t]] = bchosen[t]; }
            nforced += bforced; ndeferred += bdeferred;   // counted as if grown again
        } else
            grow(seed, best, fixed, avail, nfixed);
        for (int i : fixed) region[i] = nreg;
        nreg++;
    }

    tick("gauge");
    if (getenv("VB_HASH")) { uint64_t h = 0; for (int i = 0; i < N; i++) h = hmix(h, chosen[i]); fprintf(stderr, "HASH chosen %016llx\n", (unsigned long long)h); }
    // from here on only each kept cell's chosen labelling is needed: the others, and the tables of the search, go
    std::vector<std::vector<int>>().swap(dom); std::vector<std::vector<int>>().swap(lcls);
    std::vector<std::vector<int8_t>>().swap(cidx); std::vector<uint64_t>().swap(sup); std::vector<char>().swap(sfill); std::vector<int64_t>().swap(soff);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < N; i++) {
        if (kept[i]) { if (isos[i].size() > 1) { Iso m = isos[i][chosen[i]]; std::vector<Iso>(1, m).swap(isos[i]); chosen[i] = 0; } }
        else std::vector<Iso>().swap(isos[i]);
    }
    // ---------------- gluings and rotation transport
    // gluings (i, j, g, fi), in order of i, then of i's faces: x_i = PG[g] x_j + d, where d is the ideal vector of
    // i's face fi (the first with j) under i's labelling, so it is looked up, not stored
    std::vector<std::array<int, 4>> gl;
    auto dvec = [&](size_t e) -> V3 { int i = gl[e][0]; return Dof(i, isos[i][chosen[i]], gl[e][3]); };
    std::vector<int> rglob(N, -1);
    {
        std::vector<std::vector<std::pair<int, int>>> lj(N);     // per cell i: (j, g) with j > i
        std::vector<std::vector<int8_t>> fis(N);                  // and the face of i each is fitted on
        #pragma omp parallel for schedule(dynamic, 256)
        for (int i = 0; i < N; i++) {
            if (!kept[i]) continue;
            for (size_t f = 0; f < cells[i].nbr.size(); f++) {
                int j = cells[i].nbr[f];
                if (j <= i || !kept[j]) continue;
                int g;
                if (!cg[i].empty() && cg[i][f] != -2) g = cg[i][f];
                else g = glue(i, isos[i][chosen[i]], j, isos[j][chosen[j]], nullptr);
                if (g < 0) continue;
                int fi = 0; while (cells[i].nbr[fi] != j) fi++;           // the face the gluing is fitted on
                lj[i].push_back({j, g}); fis[i].push_back((int8_t)fi);
            }
        }
        std::vector<std::vector<int8_t>>().swap(cg);
        std::vector<int64_t> go(N + 1, 0);
        for (int i = 0; i < N; i++) go[i + 1] = go[i] + (int64_t)lj[i].size();
        gl.resize(go[N]);
        #pragma omp parallel for schedule(dynamic, 1024)
        for (int i = 0; i < N; i++) {
            for (size_t k = 0; k < lj[i].size(); k++) gl[go[i] + k] = {i, lj[i][k].first, lj[i][k].second, fis[i][k]};
            std::vector<std::pair<int, int>>().swap(lj[i]); std::vector<int8_t>().swap(fis[i]);
        }
    }
    tick("gluings");
    // each cell's gluings (other cell, gluing index), in order of the gluing index (CSR)
    std::vector<int64_t> ao(N + 1, 0);
    for (auto& e : gl) { ao[e[0] + 1]++; ao[e[1] + 1]++; }
    for (int i = 0; i < N; i++) ao[i + 1] += ao[i];
    std::vector<std::pair<int, int>> adj(ao[N]);
    {
        std::vector<int64_t> fill(ao.begin(), ao.end() - 1);
        for (size_t e = 0; e < gl.size(); e++) { adj[fill[gl[e][0]]++] = {gl[e][1], (int)e}; adj[fill[gl[e][1]]++] = {gl[e][0], (int)e}; }
    }
    int conflicts = 0;
    std::vector<int> root_of(nreg, -1);
    for (int i = 0; i < N; i++) if (kept[i] && root_of[region[i]] < 0) root_of[region[i]] = i;
    int identity = nearest_pg({1, 0, 0, 0, 1, 0, 0, 0, 1});
    const int NG = (int)PG.size();
    std::vector<int> prod(NG * NG), inv(NG);          // group tables: PG[a] PG[b] and PG[a]^T
    for (int a = 0; a < NG; a++) {
        M3 t; for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) t[3 * p + q] = PG[a][3 * q + p];
        inv[a] = nearest_pg(t);
        for (int b = 0; b < NG; b++) prod[a * NG + b] = nearest_pg(mul(PG[a], PG[b]));
    }
    for (int r = 0; r < nreg; r++) {
        int root = root_of[r];
        if (root < 0) continue;
        rglob[root] = identity;
        std::deque<int> q{root};
        while (!q.empty()) {
            int i = q.front(); q.pop_front();
            for (int64_t a = ao[i]; a < ao[i + 1]; a++) {
                auto [j, e] = adj[a];
                // x_a = R x_b + d with (a,b) = (gl.i, gl.j); transport R_j = R_i * R_ij
                int gij = gl[e][0] == i ? gl[e][2] : inv[gl[e][2]];
                int gj = prod[rglob[i] * NG + gij];
                if (rglob[j] >= 0) { if (rglob[j] != gj) conflicts++; continue; }
                rglob[j] = gj; q.push_back(j);
            }
        }
    }

    std::vector<std::pair<int, int>>().swap(adj); std::vector<int64_t>().swap(ao);
    tick("gluings_transport");
    // ---------------- output (binary if a prefix is given as argv[3], else text)
    int nk = 0; for (int i = 0; i < N; i++) nk += kept[i];
    if (argc > 3) {
        std::string pre(argv[3]);
        std::vector<int32_t> buf;                  // files are written from whole arrays
        auto put = [&](const char* ext) { FILE* fo = fopen((pre + ext).c_str(), "wb"); fwrite(buf.data(), 4, buf.size(), fo); fclose(fo); buf.clear(); };
        buf.reserve((size_t)4 * N);
        for (int i = 0; i < N; i++) { buf.push_back((int32_t)kept[i]); buf.push_back(kept[i] ? isos[i][chosen[i]].orbit : -1); buf.push_back(region[i]); buf.push_back(rglob[i]); }
        put(".cells");
        int32_t off = 0;
        buf.push_back(off);
        for (int i = 0; i < N; i++) { off += (int32_t)cells[i].nbr.size(); buf.push_back(off); }
        for (int i = 0; i < N; i++) for (int j : cells[i].nbr) buf.push_back(j);
        put(".nbrs");
        // face cycles of the cells that are not kept (the defect clusters): the cluster analysis needs the faces of
        // these cells only (neighbor lists of all cells are in .nbrs). Per cell: id, face count (-1 if the cell could
        // not be computed), then per face: neighbor, vertex count, local vertex ids.
        for (int i = 0; i < N; i++) {
            if (kept[i]) continue;
            buf.push_back(i); buf.push_back(cells[i].ok ? (int32_t)cells[i].nbr.size() : -1);
            if (!cells[i].ok) continue;
            for (size_t f = 0; f < cells[i].nbr.size(); f++) {
                buf.push_back(cells[i].nbr[f]); buf.push_back((int32_t)cells[i].cyc[f].size());
                for (int x : cells[i].cyc[f]) buf.push_back(x);
            }
        }
        put(".faces");
        // the Voronoi faces, edges, and vertices at the cells that are not kept, as rows of the cells that meet there
        // (sorted within each row; all rows sorted, without repeats): per face (cell, neighbor), per edge on exactly
        // two faces of a cell (cell, the two neighbors), per vertex on exactly three faces (cell, the three
        // neighbors); and the cells with an edge or vertex that is not (irregular). The cluster analysis reads these
        // instead of grouping the faces itself (shell.py, _Faces.structure_all).
        {
            int nth = omp_get_max_threads();
            std::vector<std::vector<std::array<int32_t, 2>>> P2(nth); std::vector<std::vector<std::array<int32_t, 3>>> P3(nth);
            std::vector<std::vector<std::array<int32_t, 4>>> P4(nth); std::vector<std::vector<int32_t>> irr(nth);
            #pragma omp parallel
            {
                int t = omp_get_thread_num();
                std::vector<int> cnt, first3; std::vector<std::array<int, 3>> e;     // (min, max, neighbor) per edge use
                #pragma omp for schedule(dynamic, 256)
                for (int i = 0; i < N; i++) {
                    if (kept[i] || !cells[i].ok) continue;
                    const Cell& c = cells[i];
                    bool bad = false;
                    int nv = 0; for (int x : c.cyc.v) nv = std::max(nv, x + 1);
                    cnt.assign(nv, 0); first3.assign(3 * (size_t)nv, 0); e.clear();
                    for (size_t f = 0; f < c.nbr.size(); f++) {
                        int nb = c.nbr[f];
                        std::array<int32_t, 2> r2 = {i, nb}; if (r2[0] > r2[1]) std::swap(r2[0], r2[1]);
                        P2[t].push_back(r2);
                        auto cy = c.cyc[f]; int m = (int)cy.size();
                        for (int k = 0; k < m; k++) {
                            int u = cy[k], v = cy[(k + 1) % m];
                            if (cnt[u] < 3) first3[3 * u + cnt[u]] = nb;
                            cnt[u]++;
                            e.push_back({std::min(u, v), std::max(u, v), nb});
                        }
                    }
                    for (int u = 0; u < nv; u++) {
                        if (cnt[u] == 0) continue;
                        if (cnt[u] != 3) { bad = true; continue; }
                        std::array<int32_t, 4> r = {i, first3[3 * u], first3[3 * u + 1], first3[3 * u + 2]};
                        std::sort(r.begin(), r.end()); P4[t].push_back(r);
                    }
                    std::sort(e.begin(), e.end(), [](const std::array<int, 3>& x, const std::array<int, 3>& y) { return x[0] != y[0] ? x[0] < y[0] : x[1] < y[1]; });
                    for (size_t a = 0; a < e.size();) {
                        size_t b = a; while (b < e.size() && e[b][0] == e[a][0] && e[b][1] == e[a][1]) b++;
                        if (b - a != 2) bad = true;
                        else { std::array<int32_t, 3> r = {i, e[a][2], e[a + 1][2]}; std::sort(r.begin(), r.end()); P3[t].push_back(r); }
                        a = b;
                    }
                    if (bad) irr[t].push_back(i);
                }
            }
            auto merge = [&](auto& parts) {
                typename std::decay<decltype(parts[0])>::type all;
                size_t tot = 0; for (auto& x : parts) tot += x.size();
                all.reserve(tot);
                for (auto& x : parts) { all.insert(all.end(), x.begin(), x.end()); std::vector<typename std::decay<decltype(x[0])>::type>().swap(x); }
                __gnu_parallel::sort(all.begin(), all.end()); all.erase(std::unique(all.begin(), all.end()), all.end());
                return all;
            };
            auto A2 = merge(P2); auto A3 = merge(P3); auto A4 = merge(P4);
            std::vector<int32_t> bads; for (auto& x : irr) bads.insert(bads.end(), x.begin(), x.end()); std::sort(bads.begin(), bads.end());
            FILE* fs = fopen((pre + ".struct").c_str(), "wb");
            int64_t hd[4] = {(int64_t)A2.size(), (int64_t)A3.size(), (int64_t)A4.size(), (int64_t)bads.size()};
            fwrite(hd, 8, 4, fs); fwrite(A2.data(), 8, A2.size(), fs); fwrite(A3.data(), 12, A3.size(), fs); fwrite(A4.data(), 16, A4.size(), fs);
            fwrite(bads.data(), 4, bads.size(), fs); fclose(fs);
        }
        std::vector<int32_t>().swap(buf);
        // gluings sorted by (i, j), so that the caller's lookup needs no sort (they are in order of i already: each
        // cell's are sorted by j, ties by gluing index)
        std::vector<size_t> go_(gl.size());
        for (size_t e = 0; e < gl.size(); e++) go_[e] = e;
        std::vector<int64_t> gs(N + 1, 0);                 // gluings of cell i (as the first cell): gs[i] .. gs[i + 1]
        for (auto& e : gl) gs[e[0] + 1]++;
        for (int i = 0; i < N; i++) gs[i + 1] += gs[i];
        #pragma omp parallel for schedule(dynamic, 4096)
        for (int i = 0; i < N; i++)
            std::sort(go_.begin() + gs[i], go_.begin() + gs[i + 1], [&](size_t x, size_t y) { return gl[x][1] != gl[y][1] ? gl[x][1] < gl[y][1] : x < y; });
        {
            FILE* fg = fopen((pre + ".glue").c_str(), "wb");
            const size_t CH = 1 << 20;                     // records per chunk (36 bytes each)
            std::vector<char> gb;
            for (size_t t0 = 0; t0 < gl.size(); t0 += CH) {
                size_t nt = std::min(CH, gl.size() - t0);
                gb.resize(nt * 36);
                #pragma omp parallel for schedule(static)
                for (long long t = 0; t < (long long)nt; t++) {
                    size_t e = go_[t0 + t]; int32_t r3[3] = {gl[e][0], gl[e][1], gl[e][2]};
                    V3 d = dvec(e); memcpy(gb.data() + 36 * t, r3, 12); memcpy(gb.data() + 36 * t + 12, d.data(), 24);
                }
                fwrite(gb.data(), 1, gb.size(), fg);
            }
            fclose(fg);
        }
        std::vector<size_t>().swap(go_); std::vector<int64_t>().swap(gs);
        tick("write_files");
        // the lab orientation of each region: the sum over its gluings of (actual neighbour vector) (ideal vector, in the
        // region's root frame)^T, which the caller turns into the best rotation (computed here, in parallel)
        {
            int nth = omp_get_max_threads();
            std::vector<std::vector<double>> acc(nth, std::vector<double>((size_t)nreg * 9, 0.0));
            #pragma omp parallel for schedule(static)
            for (long long e = 0; e < (long long)gl.size(); e++) {
                int i = gl[e][0], j = gl[e][1];
                if (rglob[i] < 0 || region[i] < 0) continue;
                V3 A = mvec(PG[rglob[i]], dvec(e));
                V3 B = minimage(V3{X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]});
                double* M = acc[omp_get_thread_num()].data() + (size_t)region[i] * 9;
                for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) M[3 * x + y] += B[x] * A[y];
            }
            std::vector<double> tot((size_t)nreg * 9, 0.0);
            for (auto& a2 : acc) for (size_t k = 0; k < tot.size(); k++) tot[k] += a2[k];
            FILE* fq = fopen((pre + ".qfit").c_str(), "wb");
            int32_t nr = nreg; fwrite(&nr, 4, 1, fq); fwrite(tot.data(), 8, tot.size(), fq);
            fclose(fq);
        }
        // the smallest loops: every triangle i < j < k of mutually glued kept cells must close (net rotation the
        // identity, net translation below the tolerance); counted here, where every frame relation is at hand
        long long n_tri = 0, n_open = 0; double worst_open = 0;
        {
            // each cell's relations as (neighbor, gluing index, direction), sorted by neighbor; R and d are read from
            // the gluing on the fly (reverse direction: R^T and -R^T d), so nothing is copied
            std::vector<int64_t> go(N + 1, 0);
            for (auto& e : gl) { go[e[0] + 1]++; go[e[1] + 1]++; }
            for (int i = 0; i < N; i++) go[i + 1] += go[i];
            std::vector<int> gdst(go[N]); std::vector<int32_t> gref(go[N]);   // gluing index + 1, negative: reverse
            {
                std::vector<int64_t> fill(go.begin(), go.end() - 1);
                for (size_t e = 0; e < gl.size(); e++) {
                    int64_t a = fill[gl[e][0]]++; gdst[a] = gl[e][1]; gref[a] = (int32_t)e + 1;
                    int64_t b = fill[gl[e][1]]++; gdst[b] = gl[e][0]; gref[b] = -((int32_t)e + 1);
                }
            }
            #pragma omp parallel for schedule(dynamic, 1024)
            for (int i = 0; i < N; i++) {
                std::vector<std::pair<int, int32_t>> v;
                for (int64_t a = go[i]; a < go[i + 1]; a++) v.push_back({gdst[a], gref[a]});
                std::sort(v.begin(), v.end());
                for (size_t q = 0; q < v.size(); q++) { gdst[go[i] + q] = v[q].first; gref[go[i] + q] = v[q].second; }
            }
            auto rel = [&](int64_t a, M3& Rm, V3& d) {          // x_i = R x_j + d for the relation at slot a
                int32_t r = gref[a]; size_t e = (size_t)(std::abs(r) - 1);
                const M3& R0 = PG[gl[e][2]];
                if (r > 0) { Rm = R0; d = dvec(e); return; }
                for (int x = 0; x < 3; x++) for (int c = 0; c < 3; c++) Rm[3 * x + c] = R0[3 * c + x];
                d = mvec(Rm, dvec(e)); for (int k = 0; k < 3; k++) d[k] = -d[k];
            };
            #pragma omp parallel for schedule(dynamic, 1024) reduction(+:n_tri, n_open) reduction(max:worst_open)
            for (int i = 0; i < N; i++)
                for (int64_t a = go[i]; a < go[i + 1]; a++) {
                    int j = gdst[a]; if (j <= i) continue;
                    M3 Ra; V3 da; rel(a, Ra, da);
                    for (int64_t b = go[j]; b < go[j + 1]; b++) {
                        int k = gdst[b]; if (k <= j) continue;
                        auto lo = gdst.begin() + go[k], hi = gdst.begin() + go[k + 1];
                        auto it = std::lower_bound(lo, hi, i);
                        if (it == hi || *it != i) continue;
                        int64_t c = it - gdst.begin();
                        M3 Rb, Rc; V3 db, dc; rel(b, Rb, db); rel(c, Rc, dc);
                        M3 Rab = mul(Ra, Rb), Rt = mul(Rab, Rc);
                        V3 t1 = mvec(Ra, db), t2 = mvec(Rab, dc);
                        double rot = 0, tr = 0;
                        for (int r = 0; r < 3; r++) for (int s2 = 0; s2 < 3; s2++) rot = std::max(rot, std::fabs(Rt[3 * r + s2] - (r == s2 ? 1.0 : 0.0)));
                        for (int r = 0; r < 3; r++) tr = std::max(tr, std::fabs(da[r] + t1[r] + t2[r]));
                        bool rot_ok = rot < 1e-3, closed = rot_ok && tr < TOL;
                        n_tri++;
                        if (!closed) { n_open++; if (rot_ok) worst_open = std::max(worst_open, tr); }
                    }
                }
        }
        tick("smallest_loops");
        printf("stats kept %d candidates %d regions %d pseudo %d conflicts %d gluings %zu forced %d deferred %d margin %.4g smallest_loops %lld open_loops %lld worst_open %.6g\n", nk, (int)order.size(), nreg, npseudo, conflicts, gl.size(), nforced, ndeferred, npseudo ? min_margin : 0.0, n_tri, n_open, worst_open);
        tick("output");
        fflush(stdout); fflush(stderr);
        _exit(0);                       // all output is written and closed: skip freeing the tessellation
    }
    printf("cells %d\n", N);
    for (int i = 0; i < N; i++) printf("%d %d %d %d %d\n", i, (int)kept[i], kept[i] ? isos[i][chosen[i]].orbit : -1, region[i], rglob[i]);
    printf("nbrs\n");
    for (int i = 0; i < N; i++) {
        printf("%d", (int)cells[i].nbr.size());
        for (int j : cells[i].nbr) printf(" %d", j);
        printf("\n");
    }
    printf("gluings %zu\n", gl.size());
    for (size_t e = 0; e < gl.size(); e++) { V3 d = dvec(e); printf("%d %d %d %.12g %.12g %.12g\n", gl[e][0], gl[e][1], gl[e][2], d[0], d[1], d[2]); }
    printf("stats kept %d candidates %d regions %d pseudo %d conflicts %d\n", nk, (int)order.size(), nreg, npseudo, conflicts);
    return 0;
}
