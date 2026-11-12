"""Reference crystal: point group, Wyckoff orbits, and one reference Voronoi cell per orbit.

Reference.  From the ideal crystal (lattice + basis + species) we get, via spglib, the point group
and the Wyckoff orbits, and from a large ideal cluster one reference Voronoi cell per orbit:
every face carries the ideal neighbour vector, every vertex its ideal position (Cartesian,
relative to the site). For each orbit we also count
    n_site : site-symmetry operations (point-group elements fixing the labelled neighbour star)
    n_comb : combinatorial automorphisms of the cell, proper plus mirror (with or without labels)
If n_comb == n_site, every combinatorial match of a sample cell onto the reference is a genuine
crystal symmetry, so ONE match suffices ("rigid" orbit). Otherwise the extra combinatorial
symmetries must be removed using the neighbours (arc consistency).

The labelling of a sample, the gluings across shared faces (P_i(x) = R P_j(x) + d for all vertices x of the face,
d = i's ideal vector to j) and the holonomies of loops are computed by the C++ engine (engine.FastConnection).
"""

import itertools

import numpy as np
import spglib

from . import voro
from .cell import Cell, all_isomorphisms


# ====================================================================== reference

def relative_gap(A, frac, stop=None):
    """Minimum relative gap of a crystal's Voronoi vertices: at each vertex (Delaunay circumcentre, radius r) the
    nearest atom other than its four is (1 + gap) r away. Vertices are taken with circumcentres in the unit cell, from
    the atoms of the cell grown by a margin m on every side; flat simplices (on the cluster's faces) are skipped. A
    vertex whose ball of radius 2r could leave the cluster (2r > m) doubles m and starts again, so every vertex used
    and its nearest other atom are those of the crystal. With stop, return as soon as one vertex has a gap below it."""
    from scipy.spatial import Delaunay
    A = np.asarray(A, float); frac = np.mod(np.asarray(frac, float), 1.0)
    vol = abs(np.linalg.det(A))
    d = float(np.cbrt(vol / len(frac)))
    h = np.array([vol / np.linalg.norm(np.cross(A[(i + 1) % 3], A[(i + 2) % 3])) for i in range(3)])  # cell heights
    m = 6 * d
    while True:
        mf = m / h                                                       # the margin in fractional coordinates
        lo, hi = np.floor(-mf).astype(int) - 1, np.ceil(1 + mf).astype(int) + 1
        sh = np.array(list(itertools.product(*[range(l, u + 1) for l, u in zip(lo, hi)])))
        F = (sh[:, None, :] + frac[None]).reshape(-1, 3)
        X = F[np.all((F > -mf) & (F < 1 + mf), axis=1)] @ A
        best, grown = np.inf, False
        for t in Delaunay(X).simplices:
            P = X[t]
            M = 2 * (P[1:] - P[0])
            if abs(np.linalg.det(M)) < 1e-6 * d ** 3:
                continue
            c = np.linalg.solve(M, (P[1:] ** 2 - P[0] ** 2).sum(1))
            fc = np.linalg.solve(A.T, c)
            if np.any(fc < -1e-9) or np.any(fc >= 1 - 1e-9):
                continue
            r = np.linalg.norm(P[0] - c)
            if 2 * r > m:
                m, grown = 2 * m, True
                break
            dist = np.linalg.norm(X - c, axis=1); dist[t] = np.inf
            best = min(best, (dist.min() - r) / r)
            if stop is not None and best < stop:
                return float(best)
        if not grown:
            return float(best)
    return float(best)


class Reference:
    def __init__(self, lattice, frac, numbers, pad=3, tol=None, symprec=1e-4, name=None):
        self.name = name
        self.A = np.asarray(lattice, float)
        frac = np.asarray(frac, float) % 1.0
        numbers = np.asarray(numbers)
        ds = spglib.get_symmetry_dataset((self.A, frac, numbers), symprec=symprec)
        self.spacegroup = ds.international
        At, Ait = self.A.T, np.linalg.inv(self.A.T)
        rots = {}
        for r in ds.rotations:
            Rc = At @ r @ Ait
            Rr = np.rint(Rc)
            rots[tuple((Rr if np.allclose(Rc, Rr, atol=1e-8) else np.round(Rc, 10)).ravel())] = 1
        self.R = np.array([np.array(k, float).reshape(3, 3) for k in rots])
        # centering translations (fractional, with the zero vector): the crystal's translations are the integer
        # vectors plus one of these, and so is the net translation of every closed loop whose net rotation is the identity
        cen = {(0.0, 0.0, 0.0)}
        for r, t in zip(ds.rotations, ds.translations):
            if np.array_equal(r, np.eye(3, dtype=int)):
                cen.add(tuple(np.round(np.asarray(t, float) % 1.0, 8) % 1.0))
        self.centerings = np.array(sorted(cen))
        orbits = ds.equivalent_atoms
        self.frac, self.numbers, self.orbit_of_site = frac, numbers, orbits
        self.scale = float(np.cbrt(abs(np.linalg.det(self.A)) / len(frac)))   # mean spacing d
        # how far the positions are from exact symmetry: the largest distance from the image of an atom under a
        # symmetry operation to the nearest atom
        err = 0.0
        for r, t in zip(ds.rotations, ds.translations):
            dif = (frac @ r.T + t)[:, None, :] - frac[None, :, :]
            dif -= np.rint(dif)
            err = max(err, float(np.linalg.norm(dif @ self.A, axis=2).min(axis=1).max()))
        self.symmetry_error = err
        # frame relations compare ideal positions only, so the tolerance need only absorb that error: rounding alone
        # (1e-6 d) for positions that are exact (a structure snapped onto its space group), else 0.2% of d
        if tol is None:
            tol = 1e-6 if err < 1e-9 * self.scale else 2e-3
        self.tol = tol * self.scale

        # large ideal cluster; reference cells are those of the unit cell at the origin
        d = self.scale
        centre = 0.5 * self.A.sum(0)
        Rc = 0.5 * np.linalg.norm(self.A.sum(0)) + np.linalg.norm(self.A, axis=1).max() + pad * 2.5 * d
        corners = np.array(list(itertools.product([-1, 1], repeat=3))) * Rc + centre
        fr = corners @ np.linalg.inv(self.A)
        lo, hi = np.floor(fr.min(0)).astype(int) - 1, np.ceil(fr.max(0)).astype(int) + 1
        shifts = np.array(list(itertools.product(*[range(l, h + 1) for l, h in zip(lo, hi)])))
        F = (shifts[:, None, :] + frac[None]).reshape(-1, 3)
        site = np.tile(np.arange(len(frac)), len(shifts))
        X = F @ self.A
        m = np.linalg.norm(X - centre, axis=1) < Rc
        F, site, X = F[m], site[m], X[m]
        off = X.min(0) - 1.0
        cells = voro.faces(X - off, (X - off).max(0) + 1.0, (0, 0, 0), threads=4)
        sig = {}

        def signature(j):
            if j not in sig:
                sig[j] = Cell(j, cells[j]).signature() if cells[j] is not None else None
            return sig[j]

        self.refs = []
        for s in sorted(set(orbits)):
            idx = np.flatnonzero((site == s) & np.all(np.isclose(F - frac[s], 0), axis=1))[0]
            c = Cell(idx, cells[idx])
            x0 = X[idx]
            dvec = {f: X[nb] - x0 for f, nb in enumerate(c.nbr)}
            P = {}
            for v, key in c.vkey.items():
                q = np.array([X[j] - x0 for j in key if j != idx])
                P[v] = np.linalg.solve(2 * q[:3], (q[:3] ** 2).sum(1))
            Flab = [(signature(nb), int(numbers[site[nb]])) for nb in c.nbr]
            mirror = c.mirror()
            n_comb = len(all_isomorphisms(c, c)) + len(all_isomorphisms(c, mirror))
            n_comb_lab = len(all_isomorphisms(c, c, Flab, Flab)) + len(all_isomorphisms(c, mirror, Flab, Flab))
            # site symmetry: point-group operations mapping the labelled neighbour star to itself
            D = np.array([dvec[f] for f in range(len(c.nbr))])
            lab = [Flab[f][1] for f in range(len(c.nbr))]
            n_site = 0
            for Rm in self.R:
                img = D @ Rm.T
                dist = np.linalg.norm(img[:, None, :] - D[None], axis=2)
                j = dist.argmin(1)
                if dist[np.arange(len(D)), j].max() < 1e-4 * d and all(lab[a] == lab[b] for a, b in enumerate(j)):
                    n_site += 1
            self.refs.append(dict(site=int(s), species=int(numbers[s]), cell=c, mirror=mirror, d=dvec, P=P,
                                  signature=c.signature(), simple=c.simple, Flab=Flab,
                                  n_comb=n_comb, n_comb_lab=n_comb_lab, n_site=n_site))
        sigs = [r["signature"] for r in self.refs]
        for r in self.refs:
            shared = sigs.count(r["signature"]) > 1
            r["rigid"] = (r["n_comb"] == r["n_site"]) and not shared
            r["rigid_lab"] = r["n_comb_lab"] == r["n_site"]

    def min_gap(self):
        """Minimum relative gap of the ideal crystal's Voronoi vertices (relative_gap)."""
        return relative_gap(self.A, self.frac)

    @property
    def all_simple(self):
        return all(r["simple"] for r in self.refs)

    def summary(self):
        out = [f"{self.name or ''} space group {self.spacegroup}, {len(self.R)} point-group operations, "
               f"{len(self.refs)} orbits, d = {self.scale:.3f}"]
        for r in self.refs:
            kind = "rigid" if r["rigid"] else ("rigid with labels" if r["rigid_lab"] else "needs neighbours")
            out.append(f"  site {r['site']} (Z={r['species']}): {len(r['cell'].cyc)} faces, simple={r['simple']}, "
                       f"|Aut| comb {r['n_comb']} (labelled {r['n_comb_lab']}), site symmetry {r['n_site']}: {kind}")
        return "\n".join(out)

    def lattice_coords(self, t):
        return np.linalg.solve(self.A.T, t)


def from_cif(path, **kw):
    from pymatgen.core import Structure
    s = Structure.from_file(path)
    return Reference(s.lattice.matrix, s.frac_coords, [site.specie.Z for site in s], name=path.split("/")[-1], **kw)


def bcc(a=1.0):
    return Reference(np.eye(3) * a, [[0, 0, 0], [.5, .5, .5]], [1, 1], name="bcc")


def a15(a=1.0, species=(50, 41)):
    nb = [[.25, 0, .5], [.75, 0, .5], [.5, .25, 0], [.5, .75, 0], [0, .5, .25], [0, .5, .75]]
    return Reference(np.eye(3) * a, [[0, 0, 0], [.5, .5, .5]] + nb, [species[0]] * 2 + [species[1]] * 6, name="A15")


# ====================================================================== connection
