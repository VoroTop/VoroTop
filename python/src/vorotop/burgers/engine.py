"""Fast path: the combinatorial connection computed by the C++ engine (tools/vb_engine.cc).

FastConnection exposes the interface used by loops.py (pts, L, periodic, kept, cells[i].nbr,
minimage, holonomy, burgers, ref) but holds only compact tables: per cell (kept, orbit, region,
transport rotation) and per glued face pair (R index, d)."""

import os
import re
import subprocess
import tempfile
from collections import deque

import numpy as np

from . import voro

_HERE = os.path.dirname(os.path.abspath(__file__))
VB_ENGINE = os.environ.get("VOROTOP_BURGERS_ENGINE", os.path.join(_HERE, "bin", "vb_engine"))   # built and installed with the package
_FNV, _PRIME, _MASK = 1469598103934665603, 1099511628211, (1 << 64) - 1


def sig_hash(signature):
    """Same hash as vb_engine: FNV-1a over the sorted face sizes, then the face count."""
    h = _FNV
    for x in signature:
        h = ((h ^ x) * _PRIME) & _MASK
    return ((h ^ len(signature)) * _PRIME) & _MASK


def write_reference(ref, fn):
    with open(fn, "w") as f:
        f.write(f"{len(ref.refs)} {len(ref.R)} {ref.tol:.12g}\n")
        for R in ref.R:
            f.write(" ".join(f"{x:.15g}" for x in R.ravel()) + "\n")
        for r in ref.refs:
            mode = 0 if r["rigid"] else (1 if r["rigid_lab"] else 2)
            mode += 4 if r.get("glue_by_neighbors") else 0
            c = r["cell"]
            nv = max(max(cy) for cy in c.cyc) + 1
            f.write(f"ref {r['species']} {len(c.cyc)} {nv} {mode}\n")
            for F, cy in enumerate(c.cyc):
                sig, sp = r["Flab"][F]
                # label 0 matches any neighbor (degenerate references); an empty label of an ordinary reference matches none
                h = sig_hash(sig) if sig else (0 if r.get("glue_by_neighbors") else 1)
                f.write(f"{h} {sp} {len(cy)} " + " ".join(map(str, cy)) + "\n")
            for F in range(len(c.cyc)):
                f.write(" ".join(f"{x:.15g}" for x in r["d"][F]) + "\n")
            for v in range(nv):
                f.write(" ".join(f"{x:.15g}" for x in r["P"][v]) + "\n")


class _NbrCell:
    __slots__ = ("nbr",)

    def __init__(self, nbr):
        self.nbr = nbr


class _Cells:
    """Neighbor lists of all cells, served from the engine's compact arrays on demand (a Python list per atom would
    cost several hundred bytes per atom)."""

    def __init__(self, offs, data):
        self._o, self._d = offs, data

    def __len__(self):
        return len(self._o) - 1

    def __getitem__(self, i):
        return _NbrCell(self._d[self._o[i]:self._o[i + 1]].tolist())

    def __iter__(self):
        return (self[i] for i in range(len(self)))


def _unique_rows(A):
    """Sorted unique rows of an integer array (lexicographic, as sorted tuples), faster than np.unique(axis=0)."""
    if not len(A):
        return A
    o = np.lexsort(A.T[::-1]); A = A[o]
    return A[np.r_[True, np.any(A[1:] != A[:-1], axis=1)]]


class _Faces:
    """Faces of the cells that are not kept, from the engine's own tessellation, as flat arrays: per face its cell,
    neighbor, and vertex cycle (local vertex ids). faces[i] gives [(neighbor, [vertex ids]), ...] for such a cell and
    None for every other cell (the cluster analysis needs no others); structure(members) gives the Voronoi faces,
    edges, and vertices at the member cells (see shell._structure) by sorting and grouping these arrays."""

    def __init__(self, n, fcell, fnb, foff, flen, vflat, missing):
        self._n = n
        self.fcell, self.fnb, self.foff, self.flen, self.vflat = fcell, fnb, foff, flen, vflat
        self.missing = missing
        self.rows = None                   # the engine's rows (see _read_struct), when it wrote them
        order = np.argsort(fcell, kind="stable")
        self._order, self._sorted = order, fcell[order]

    def __len__(self):
        return self._n

    def _faces_of(self, cells):
        cells = np.asarray(cells, np.int64)
        lo = np.searchsorted(self._sorted, cells, "left"); hi = np.searchsorted(self._sorted, cells, "right")
        return np.concatenate([self._order[a:b] for a, b in zip(lo, hi)]) if len(cells) else np.zeros(0, np.int64)

    def __getitem__(self, i):
        fi = self._faces_of([int(i)])
        if not len(fi):
            return None
        return [(int(self.fnb[f]), self.vflat[self.foff[f]:self.foff[f] + self.flen[f]].tolist()) for f in fi]

    def structure_all(self, label):
        """structure() for every cell at once, each row with the cluster it belongs to: label[i] is the cluster of
        atom i (-1 for crystalline atoms); a Voronoi edge or vertex at a defect cell belongs to that cell's cluster,
        since defect cells that share one share a face. Returns (pairs, triples, quads, their clusters, bad cells),
        rows sorted and unique."""
        rows_lab = lambda A: np.max(np.where(A >= 0, label[np.maximum(A, 0)], -1), axis=1)
        if self.rows is not None:
            # the engine's rows are those of every cell that is not kept; they are the rows below whenever every such
            # cell belongs to a cluster (no mask, no wall cells), else they are computed here from the faces
            P, T, Q, irregular = self.rows
            if (label[self.fcell] >= 0).all() and all(label[i] >= 0 for i in self.missing):
                bad = set(int(i) for i in irregular) | set(int(i) for i in self.missing)
                return P, T, Q, rows_lab(P), rows_lab(T), rows_lab(Q), bad
        fi = np.arange(len(self.fcell))
        fi = fi[label[self.fcell] >= 0]
        cell, nb, m, o = self.fcell[fi], self.fnb[fi], self.flen[fi], self.foff[fi]
        bad = set(int(i) for i in self.missing if label[i] >= 0)
        pairs = _unique_rows(np.sort(np.c_[cell, nb], axis=1))
        rep = np.repeat(np.arange(len(fi)), m)
        k = np.arange(len(rep)) - np.repeat(np.cumsum(m) - m, m)
        u = self.vflat[o[rep] + k]; v = self.vflat[o[rep] + (k + 1) % m[rep]]
        c, w = cell[rep], nb[rep]
        VB = np.int64(1 << 12)                        # local vertex ids are small
        def groups(key, size):
            srt = np.argsort(key, kind="stable"); K = key[srt]
            new = np.r_[True, K[1:] != K[:-1]]
            start = np.flatnonzero(new); cnt = np.diff(np.r_[start, len(K)])
            bad.update(int(x) for x in c[srt][start[cnt != size]])
            okg = start[cnt == size]
            nbs = np.stack([w[srt][okg + j] for j in range(size)], 1)
            return _unique_rows(np.sort(np.c_[c[srt][okg], nbs], axis=1))
        if len(c):
            if (np.maximum(u, v) >= VB).any():
                raise ValueError("a cell with too many vertices for structure_all")
            triples = groups((c.astype(np.int64) * VB + np.minimum(u, v)) * VB + np.maximum(u, v), 2)
            quads = groups(c.astype(np.int64) * VB + u, 3)
        else:
            triples, quads = np.zeros((0, 3), np.int64), np.zeros((0, 4), np.int64)
        return pairs, triples, quads, rows_lab(pairs), rows_lab(triples), rows_lab(quads), bad

    def structure(self, members):
        """(pairs (P, 2), triples (T, 3), quads (Q, 4), bad cells), rows sorted and unique, as shell._structure."""
        members = np.asarray(members, np.int64)
        bad = set(int(i) for i in members if i in self.missing)
        fi = self._faces_of(members)
        cell, nb, m, o = self.fcell[fi], self.fnb[fi], self.flen[fi], self.foff[fi]
        pairs = _unique_rows(np.sort(np.c_[cell, nb], axis=1))
        rep = np.repeat(np.arange(len(fi)), m)
        k = np.arange(len(rep)) - np.repeat(np.cumsum(m) - m, m)
        u = self.vflat[o[rep] + k]; v = self.vflat[o[rep] + (k + 1) % m[rep]]
        c, w = cell[rep], nb[rep]

        def groups(keys, size):
            srt = np.lexsort(keys[::-1]); K = np.stack(keys, 1)[srt]
            new = np.r_[True, np.any(K[1:] != K[:-1], axis=1)]
            start = np.flatnonzero(new); cnt = np.diff(np.r_[start, len(K)])
            bad.update(int(x) for x in K[start[cnt != size], 0])
            ok = start[cnt == size]
            nbs = np.stack([w[srt][ok + j] for j in range(size)], 1)
            return _unique_rows(np.sort(np.c_[K[ok, 0], nbs], axis=1))

        triples = groups([c, np.minimum(u, v), np.maximum(u, v)], 2) if len(c) else np.zeros((0, 3), np.int64)
        quads = groups([c, u], 3) if len(c) else np.zeros((0, 4), np.int64)
        return pairs, triples, quads, bad


def _read_struct(path):
    """The engine's rows of the Voronoi faces, edges, and vertices at the cells that are not kept (vb_engine .struct:
    int64 counts of pairs, triples, quads, irregular cells; then the rows as int32, sorted, without repeats), or None."""
    if not os.path.exists(path):
        return None
    raw = np.fromfile(path, np.int8)
    n2, n3, n4, nb = np.frombuffer(raw[:32].tobytes(), np.int64)
    a = np.frombuffer(raw[32:].tobytes(), np.int32)
    P = a[:2 * n2].reshape(-1, 2).astype(np.int64); a = a[2 * n2:]
    T = a[:3 * n3].reshape(-1, 3).astype(np.int64); a = a[3 * n3:]
    Q = a[:4 * n4].reshape(-1, 4).astype(np.int64); a = a[4 * n4:]
    return P, T, Q, a[:nb].astype(np.int64)


def _read_faces(path, n):
    raw = np.fromfile(path, dtype=np.int32)
    fcell, fnb, foff, flen, missing = [], [], [], [], set()
    k, L = 0, len(raw)
    rl = raw.tolist()
    while k < L:
        i, nf = rl[k], rl[k + 1]; k += 2
        if nf < 0:
            missing.add(i); continue
        for _ in range(nf):
            fcell.append(i); fnb.append(rl[k]); flen.append(rl[k + 1]); foff.append(k + 2); k += 2 + rl[k + 1]
    return _Faces(n, np.array(fcell, np.int64), np.array(fnb, np.int64), np.array(foff, np.int64),
                  np.array(flen, np.int64), raw.astype(np.int64), missing)


class FastConnection:
    def __init__(self, pts, L, periodic, ref, species=None, threads=None, orientation=None):
        self.pts, self.L, self.periodic, self.ref = np.asarray(pts, float), np.asarray(L, float), tuple(periodic), ref
        # a 3x3 L is a fully periodic triclinic box: rows are the lattice vectors a, b, c in restricted-triclinic form
        # (a = (lx,0,0), b = (xy,ly,0), c = (xz,yz,lz)), as in LAMMPS
        self.tri = self.L.shape == (3, 3)
        if self.tri:
            H = self.L
            if abs(H[0, 1]) + abs(H[0, 2]) + abs(H[1, 2]) > 1e-9 * np.abs(H).max():
                raise ValueError("triclinic box must be restricted-triclinic: a = (lx,0,0), b = (xy,ly,0)")
            if not all(periodic):
                raise ValueError("a triclinic box must be periodic in all three directions")
        self.orientation = np.eye(3) if orientation is None else np.asarray(orientation, float)
        if not np.all(np.isfinite(self.pts)):
            raise ValueError("non-finite particle coordinates")
        n = len(self.pts)
        with tempfile.TemporaryDirectory() as tmp:
            rf, sf = os.path.join(tmp, "ref.txt"), os.path.join(tmp, "sample.bin")
            write_reference(ref, rf)
            with open(sf, "wb") as f:          # binary sample (see tools/vb_engine.cc)
                np.array([n], np.int32).tofile(f)
                if self.tri:                    # lengths, flag 2 = triclinic, then the tilts xy, xz, yz
                    np.array([H[0, 0], H[1, 1], H[2, 2]], np.float64).tofile(f)
                    np.array([2, 1, 1, int(species is not None)], np.int32).tofile(f)
                    np.array([H[1, 0], H[2, 0], H[2, 1]], np.float64).tofile(f)
                else:
                    np.asarray(L, np.float64).tofile(f)
                    np.array([int(periodic[0]), int(periodic[1]), int(periodic[2]), int(species is not None)], np.int32).tofile(f)
                np.ascontiguousarray(self.pts, np.float64).tofile(f)
                if species is not None:
                    np.asarray(species, np.int32).tofile(f)
            env = dict(os.environ, OMP_NUM_THREADS=str(threads or voro.default_threads()))
            pre = os.path.join(tmp, "out")
            run = subprocess.run([VB_ENGINE, rf, sf, pre], capture_output=True, text=True, env=env, check=True)
            self.stats, self.engine_log = run.stdout.strip(), run.stderr
            self._load(pre, n)

    def _load(self, pre, n):
        """Read the engine's output files (also written by vb_closepacked) and set up lookups and region frames."""
        R = self.ref.R
        tab = np.fromfile(pre + ".cells", dtype=np.int32).reshape(n, 4)
        nb = np.fromfile(pre + ".nbrs", dtype=np.int32)
        rec = np.dtype([("ijg", np.int32, 3), ("d", np.float64, 3)])
        g = np.fromfile(pre + ".glue", dtype=rec)
        qfit = np.fromfile(pre + ".qfit", dtype=np.uint8) if os.path.exists(pre + ".qfit") else None
        self.defect_faces = _read_faces(pre + ".faces", n)
        self.defect_faces.rows = _read_struct(pre + ".struct")
        self.kept = tab[:, 0].astype(bool)
        self.orbit, self.region, rg = tab[:, 1], tab[:, 2], tab[:, 3]
        offs, data = nb[:n + 1], nb[n + 1:]
        self._nb_offs, self._nb_data = offs.astype(np.int64), data.astype(np.int64)
        self.cells = _Cells(self._nb_offs, self._nb_data)
        I, J, Gi = g["ijg"][:, 0], g["ijg"][:, 1], g["ijg"][:, 2]
        D = g["d"]
        self._I, self._J, self._Gi, self._D = I, J, Gi, D
        # gluing lookup by sorted keys i*n + j (dicts over all faces are slow to build in Python)
        key = I.astype(np.int64) * n + J
        if len(key) < 2 or np.all(key[1:] > key[:-1]):    # the engine writes them sorted
            self._korder = np.arange(len(key)); self._keys = key
        else:
            self._korder = np.argsort(key); self._keys = key[self._korder]
        self._n = n
        self._Rg = rg
        # lab orientation of each region's root frame, from the gluing vectors (vectorised)
        # (sums of outer products per region, accumulated in chunks: no array over all gluings is formed)
        self.Qlab, self.G = {}, {}
        idx = np.flatnonzero(rg[I] >= 0)
        rids = np.unique(self.region[I[idx]])
        Macc = np.zeros((len(rids), 3, 3))
        if qfit is not None:                          # the sums, computed by the engine
            nr = int(np.frombuffer(qfit[:4].tobytes(), np.int32)[0])
            tot = np.frombuffer(qfit[4:4 + 72 * nr].tobytes(), np.float64).reshape(nr, 3, 3)
            Macc = tot[rids]
            idx = idx[:0]
        for c0 in range(0, len(idx), 1 << 20):
            c = idx[c0:c0 + (1 << 20)]
            A = np.einsum("nab,nb->na", R[rg[I[c]]], D[c])
            B = self.minimage(self.pts[J[c]] - self.pts[I[c]])
            k = np.searchsorted(rids, self.region[I[c]])
            for x in range(3):
                for y in range(3):
                    Macc[:, x, y] += np.bincount(k, weights=B[:, x] * A[:, y], minlength=len(rids))
        for kr, rid in enumerate(rids):
            U, S, Vt = np.linalg.svd(Macc[kr])
            Q = U @ Vt
            self.Qlab[int(rid)] = Q
            errs = [np.linalg.norm(Q @ Gm.T - self.orientation) for Gm in R]
            self.G[int(rid)] = R[int(np.argmin(errs))]

    def _find(self, i, j):
        q = i * self._n + j
        p = int(np.searchsorted(self._keys, q))
        return int(self._korder[p]) if p < len(self._keys) and self._keys[p] == q else None

    @property
    def _row(self):
        if not hasattr(self, "_rowdict"):
            self._rowdict = dict(zip(zip(self._I.tolist(), self._J.tolist()), range(len(self._I))))
        return self._rowdict

    @property
    def Rglob(self):
        """Rotation of each kept cell's frame into its region's root frame."""
        if not hasattr(self, "_Rglob"):
            rg = self._Rg
            self._Rglob = {int(i): self.ref.R[rg[i]] for i in np.flatnonzero(self.kept & (rg >= 0))}
        return self._Rglob

    @property
    def region_of(self):
        if not hasattr(self, "_region_of"):
            self._region_of = {i: int(self.region[i]) for i in self.Rglob}
        return self._region_of

    @property
    def _g(self):
        """All gluings as a dict (built on demand; both directions)."""
        if not hasattr(self, "_gdict"):
            R = self.ref.R
            out = {}
            for (i, j), k in self._row.items():
                Rm, d = R[self._Gi[k]], self._D[k]
                out[(i, j)] = (Rm, d)
                out[(j, i)] = (Rm.T, -Rm.T @ d)
            self._gdict = out
        return self._gdict

    def minimage(self, x):
        x = np.array(x, float)
        if self.tri:                                 # reduce along c, then b, then a (as LAMMPS and the engine do)
            H = self.L
            for k in (2, 1, 0):
                x -= np.round(x[..., k] / H[k, k])[..., None] * H[k]
            # in a strongly tilted box the result can be one box vector away from the nearest image; the 26
            # neighboring images are checked (as in the engine), and the result is unchanged where it already is one
            T = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if i or j or k], float) @ H
            n0 = np.einsum("...k,...k->...", x, x)
            # a vector no longer than half of lam, a lower bound on every nonzero box vector, is already the nearest
            # image (|x + v| >= |v| - |x| >= lam / 2 >= |x|): only longer ones need the 26 images. lam = the smallest
            # singular value of H, since |n H| >= sigma_min |n| >= sigma_min for every nonzero integer n
            if not hasattr(self, "_lam2"):
                self._lam2 = float(np.linalg.svd(H, compute_uv=False).min()) ** 2
            far = n0 > 0.25 * self._lam2 * (1 - 1e-9)
            if x.ndim == 1 or not far.all():
                if x.ndim == 1:
                    if not far:
                        return x
                else:
                    sub = self.minimage_far(x[far], n0[far], T)
                    x[far] = sub
                    return x
            for v in T:
                y = x + v
                n1 = np.einsum("...k,...k->...", y, y)
                m = n1 < n0 * (1 - 1e-12)
                x = np.where(m[..., None], y, x); n0 = np.where(m, n1, n0)
            return x
        for k in range(3):
            if self.periodic[k]:
                x[..., k] -= self.L[k] * np.round(x[..., k] / self.L[k])
        return x

    @staticmethod
    def minimage_far(x, n0, T):
        for v in T:
            y = x + v
            n1 = np.einsum("...k,...k->...", y, y)
            m = n1 < n0 * (1 - 1e-12)
            x = np.where(m[..., None], y, x); n0 = np.where(m, n1, n0)
        return x

    def g_many(self, I, J):
        """g for many pairs at once: (found, R (k, 3, 3), d (k, 3)), with the identity and zero where not found."""
        I, J = np.asarray(I, np.int64), np.asarray(J, np.int64)
        R = np.asarray(self.ref.R)
        Rout = np.tile(np.eye(3), (len(I), 1, 1)); dout = np.zeros((len(I), 3)); found = np.zeros(len(I), bool)
        for fwd in (True, False):
            q = (I * self._n + J) if fwd else (J * self._n + I)
            p = np.minimum(np.searchsorted(self._keys, q), len(self._keys) - 1)
            hit = (self._keys[p] == q) & ~found
            k = self._korder[p[hit]]
            Rm, d = R[self._Gi[k]], self._D[k]
            if fwd:
                Rout[hit], dout[hit] = Rm, d
            else:
                Rout[hit] = np.transpose(Rm, (0, 2, 1)); dout[hit] = -np.einsum("kba,kb->ka", Rm, d)
            found |= hit
        return found, Rout, dout

    def g(self, i, j):
        k = self._find(i, j)
        if k is not None:
            return self.ref.R[self._Gi[k]], self._D[k]
        k = self._find(j, i)
        if k is not None:
            Rm = self.ref.R[self._Gi[k]]
            return Rm.T, -Rm.T @ self._D[k]
        return None

    def holonomy(self, loop):
        loop = list(loop)
        if loop[0] != loop[-1]:
            loop.append(loop[0])
        R, t = np.eye(3), np.zeros(3)
        for a, b in zip(loop[:-1], loop[1:]):
            gg = self.g(a, b)
            if gg is None:
                return None
            t = t + R @ gg[1]
            R = R @ gg[0]
        return R, t

    def main_region(self):
        """The crystalline region with the most perfect cells (the bulk crystal)."""
        if not hasattr(self, "_main_region"):
            reg = np.asarray(self.region)
            ids, cnt = np.unique(reg[reg >= 0], return_counts=True)
            self._main_region = int(ids[np.argmax(cnt)]) if len(ids) else -1
        return self._main_region

    def burgers(self, loop):
        # a closed loop can be read from any of its cells; read it from a cell of the bulk crystal when it passes
        # through one, so that the vector is expressed in the bulk crystal's setting rather than that of an
        # isolated pocket of perfect cells (each region has its own, arbitrary setting)
        cyc = list(loop[:-1]) if len(loop) > 1 and loop[0] == loop[-1] else list(loop)
        if len(cyc) > 2:
            main = self.main_region()
            k = next((i for i, c in enumerate(cyc) if self.region[c] == main and self.kept[c]), 0)
            loop = cyc[k:] + cyc[:k]
        h = self.holonomy(loop)
        start = loop[0]
        if h is None or not self.kept[start] or self._Rg[start] < 0:
            return None
        R, t = h
        b_root = -(self.ref.R[self._Rg[start]] @ t)
        rid = int(self.region[start])
        b_cr = self.G[rid] @ b_root
        return dict(b_crystal=b_cr, lattice=self.ref.lattice_coords(b_cr), b_lab=self.Qlab[rid] @ b_root,
                    R_identity=bool(np.abs(R - np.eye(3)).max() < 1e-3), region=rid)

    def gluing_stats(self):
        """(glued kept-kept face pairs, kept-kept face pairs NOT glued by a lattice symmetry)."""
        ok = bad = 0
        for i in np.flatnonzero(self.kept):
            for j in self.cells[i].nbr:
                if j > i and self.kept[j]:
                    if (i, j) in self._row or (j, i) in self._row:
                        ok += 1
                    else:
                        bad += 1
        return ok, bad

    def loop_closure(self):
        """The smallest loops: every triangle of mutually glued kept cells. Returns (number of triangles, number that
        do not close, largest net translation of those whose rotation closes). Counted by the engine (all frame
        relations are at hand there); a loop closes when its net rotation is the identity and its net translation is
        below the tolerance of the frame relations."""
        st = dict(re.findall(r"(\w+) (-?[\d.]+(?:e[-+]?\d+)?)", self.stats))
        return int(st["smallest_loops"]), int(st["open_loops"]), float(st["worst_open"])

    def _loop_closure_numpy(self, chunk=1 << 22):
        """The same count in numpy (for checking the engine's)."""
        n = self._n
        R = self.ref.R
        I, J = self._I.astype(np.int64), self._J.astype(np.int64)
        Ra, Da = R[self._Gi], self._D
        RaT = np.transpose(Ra, (0, 2, 1))
        # directed gluings in both directions, sorted by key i*n + j
        src = np.concatenate([I, J]); dst = np.concatenate([J, I])
        Rd = np.concatenate([Ra, RaT]); Dd = np.concatenate([Da, -np.einsum("nab,nb->na", RaT, Da)])
        key = src * n + dst
        o = np.argsort(key); key, src, dst, Rd, Dd = key[o], src[o], dst[o], Rd[o], Dd[o]
        start = np.searchsorted(src, np.arange(n + 1))
        up = np.flatnonzero(src < dst)                          # edges i -> j with i < j
        deg = start[dst[up] + 1] - start[dst[up]]
        n_tri = n_open = 0
        worst = 0.0
        for c0 in range(0, len(up), max(1, chunk // 16)):
            e = up[c0:c0 + max(1, chunk // 16)]
            d = deg[c0:c0 + len(e)]
            rep = np.repeat(e, d)                               # edge i -> j, once per neighbor k of j
            off = np.arange(int(d.sum())) - np.repeat(np.cumsum(d) - d, d)
            jk = np.repeat(start[dst[e]], d) + off              # edge j -> k
            k = dst[jk]
            i = src[rep]
            m = k > dst[rep]                                    # each triangle once: i < j < k
            rep, jk, i, k = rep[m], jk[m], i[m], k[m]
            q = k * n + i                                       # edge k -> i, if glued
            p = np.searchsorted(key, q)
            ok = (p < len(key)) & (key[np.minimum(p, len(key) - 1)] == q)
            rep, jk, ki = rep[ok], jk[ok], p[ok]
            Rt = np.einsum("nab,nbc->nac", np.einsum("nab,nbc->nac", Rd[rep], Rd[jk]), Rd[ki])
            t = Dd[rep] + np.einsum("nab,nb->na", Rd[rep], Dd[jk]) \
                + np.einsum("nab,nb->na", np.einsum("nab,nbc->nac", Rd[rep], Rd[jk]), Dd[ki])
            rot_ok = np.abs(Rt - np.eye(3)).max(axis=(1, 2)) < 1e-3
            tr = np.abs(t).max(axis=1)
            closed = rot_ok & (tr < self.ref.tol)
            n_tri += len(closed)
            n_open += int((~closed).sum())
            if (rot_ok & ~closed).any():
                worst = max(worst, float(tr[rot_ok & ~closed].max()))
        return n_tri, n_open, worst

    def flatness(self):
        """Holonomy around every triangle of mutually glued kept cells (= around Voronoi edges)."""
        nb = {}
        for (i, j) in self._g:
            nb.setdefault(i, set()).add(j)
        n = bad = 0
        for i, s in nb.items():
            for j in s:
                if j <= i:
                    continue
                for k in s & nb.get(j, set()):
                    if k <= j:
                        continue
                    h = self.holonomy([i, j, k])
                    n += 1
                    if not (np.abs(h[0] - np.eye(3)).max() < 1e-3 and np.abs(h[1]).max() < self.ref.tol):
                        bad += 1
        return n, bad, 0
