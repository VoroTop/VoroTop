"""Close-packed crystals (fcc, hcp, and their stacking faults) through the C++ tool vb_closepacked, behind the
FastConnection interface.

The perfect fcc and hcp cells are degenerate: some of their Voronoi vertices are shared by more than four cells, and
any perturbation splits them, so a thermal cell has extra faces and its topology alone need not say which of its faces
join first neighbours. A perturbation small enough to leave the generic vertices intact changes the topology only at
the degenerate ones, so the cells it can produce are enumerated in full from the resolutions of those vertices
(VoroTypes): for each topology (Weinberg vector), the sets of its faces that join the twelve first neighbours, over all
its isomorphisms with the enumerated cells. Shipped tables: data/cp_fcc.json.gz (2,815 topologies),
data/cp_hcp_1.6237.json.gz (hcp at the axial ratio of Mg, 9,556 topologies), and data/cp_hcp_ideal.json.gz.

The tool (src/cpp/vb_closepacked.cc, where the steps are described) reads each cell's first neighbours from the table
with neighbour agreement, labels every decided cell whose twelve neighbours form an ideal shell (the templates below),
checks that neighbouring labelled cells are related by crystal symmetries, and writes the same files as vb_engine, so
everything after it (clusters, boundary surfaces, Burgers vectors, branch labels) is shared with the stable route."""

import gzip
import itertools
import json
import os
import subprocess
import tempfile

import numpy as np
import spglib

from . import voro
from .engine import FastConnection, _Cells

_HERE = os.path.dirname(os.path.abspath(__file__))
VB_CLOSEPACKED = os.environ.get("VOROTOP_BURGERS_CLOSEPACKED", os.path.join(_HERE, "bin", "vb_closepacked"))
DATA = os.path.join(_HERE, "data")
TABLES = {"fcc": "cp_fcc.json.gz", "hcp": "cp_hcp_1.6237.json.gz", "hcp_ideal": "cp_hcp_ideal.json.gz"}


def _templates(named):
    """{kind: (V, adjacency)} from shells of 12 vectors; two shell members are adjacent if they are first neighbours of
    each other (their distance is a bond length of the shell)."""
    out = {}
    for kind, V in named.items():
        V = np.asarray(V, float)
        bl = np.unique(np.round(np.linalg.norm(V, axis=1), 6))
        adj = {i: set() for i in range(12)}
        for i, j in itertools.combinations(range(12), 2):
            if np.min(np.abs(bl - np.linalg.norm(V[i] - V[j]))) < 0.02 * bl.min():
                adj[i].add(j); adj[j].add(i)
        out[kind] = (V, adj)
    return out


class FCCReference:
    """Ideal fcc crystal in units of the conventional cell a: ideal vectors are integers in units of a/2, so lattice
    coordinates (in units of a) are v / 2. Templates: the cuboctahedral shell, and the hcp-like shell of
    stacking-fault atoms in both stacking senses (the 6 first neighbours in the (111) plane, the 3 above it, and their
    mirror images through the plane; the other sense is its rotation by 180 degrees about [111], which is not a cubic
    symmetry)."""

    def __init__(self, a):
        self.a = a
        self.scale = (a ** 3 / 4) ** (1 / 3)
        self.A = a * np.eye(3)                     # the cubic cell; lattice coordinates are in its units
        self.centerings = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
        self.R = np.array([np.diag(s) @ np.eye(3, dtype=int)[list(p)] for p in itertools.permutations(range(3))
                           for s in itertools.product([1, -1], repeat=3)])
        shell = np.array([v for v in itertools.product([-1, 0, 1], repeat=3) if sum(map(abs, v)) == 2], float)
        n = np.ones(3) / np.sqrt(3)
        top = shell[shell @ n > 1e-9]
        hcp = np.vstack([shell[np.abs(shell @ n) < 1e-9], top, top - 2 * np.outer(top @ n, n)])
        self.templates = _templates({"fcc": shell, "hcp": hcp, "hcp2": -hcp})

    def lattice_coords(self, t):
        return np.asarray(t, float) / 2.0


class HCPReference:
    """Ideal hcp crystal (lattice a, c). Templates: the hcp shells of the two sites, and the fcc-like shells of atoms
    at basal stacking faults (the in-plane hexagon, the upper triad, and its inversion through the atom)."""

    def __init__(self, a, c):
        self.a, self.c = a, c
        self.A = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])
        frac = np.array([[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]])
        X = (np.array(list(itertools.product(range(-2, 3), repeat=3)))[:, None, :] + frac[None]).reshape(-1, 3) @ self.A
        x0 = frac[0] @ self.A
        d = np.linalg.norm(X - x0, axis=1)
        shell = X[(d > 1e-6) & (d <= np.sort(d)[12] + 1e-6)] - x0
        assert len(shell) == 12
        up = shell[shell[:, 2] > 1e-6]
        fccsh = np.vstack([shell[np.abs(shell[:, 2]) < 1e-6], up, -up])
        self.templates = _templates({"hcp": shell, "hcp2": -shell, "fcc": fccsh, "fcc2": -fccsh})
        ds = spglib.get_symmetry_dataset((self.A, frac, [1, 1]), symprec=1e-5)
        At, Ait = self.A.T, np.linalg.inv(self.A.T)
        self.R = np.unique(np.round([At @ r @ Ait for r in ds.rotations], 10).reshape(-1, 9), axis=0).reshape(-1, 3, 3)
        self.scale = float(np.cbrt(abs(np.linalg.det(self.A)) / 2))
        self.centerings = np.zeros((1, 3))         # a primitive cell

    def lattice_coords(self, t):
        return np.linalg.solve(self.A.T, np.asarray(t, float))


def _write_table(names, path):
    """The topology tables (shipped names or paths), merged, as int32: entries; per entry its faces and readings."""
    out = []
    n = 0
    for name in names:
        p = os.path.join(DATA, TABLES[name]) if name in TABLES else name
        with gzip.open(p, "rt") as fh:
            for _, ref, sets in json.load(fh):
                out.append(len(ref))
                for f in ref:
                    out.append(len(f)); out.extend(f)
                out.append(len(sets))
                for s in sets:
                    out.append(len(s)); out.extend(s)
                n += 1
    np.array([n] + out, np.int32).tofile(path)


def _write_params(path, R, templates, order, ambiguous, perfect):
    with open(path, "w") as f:
        f.write(f"{len(R)}\n")
        for M in R:
            f.write(" ".join(f"{x:.15g}" for x in np.asarray(M, float).ravel()) + "\n")
        f.write(f"{len(order)}\n")
        for k in order:
            V, adj = templates[k]
            f.write(f"{k} {int(k in ambiguous)} {int(k in perfect)}\n")
            for v in np.asarray(V, float):
                f.write(" ".join(f"{x:.15g}" for x in v) + "\n")
            for a in range(12):
                f.write(" ".join(str(int(b in adj[a])) for b in range(12)) + "\n")


class CPConnection(FastConnection):
    """kind: "fcc" (cubic fcc; partials: stacking-fault atoms with hcp shells may be labelled) or "hcp" (a, c;
    faults: fcc-like shells of basal faults are kept). Orthogonal, fully periodic boxes."""

    def __init__(self, pts, L, kind, a, c=None, threads=None, orientation=None, partials=True, faults=True, tables=None):
        self.pts, self.L, self.periodic = np.asarray(pts, float), np.asarray(L, float), (True, True, True)
        self.tri = False
        self.orientation = np.eye(3) if orientation is None else np.asarray(orientation, float)
        if kind == "fcc":
            self.ref = FCCReference(a)
            order = ["fcc", "hcp", "hcp2"] if partials else ["fcc"]
            ambiguous, perfect = {"hcp", "hcp2"}, set(order)
            tables = tables or ("fcc",)
        else:
            self.ref = HCPReference(a, c)
            order = ["hcp", "hcp2"] + (["fcc", "fcc2"] if faults else [])
            ambiguous, perfect = set(order), set(order)
            tables = tables or ("hcp",)
        self.kind_names = order
        n = len(self.pts)
        with tempfile.TemporaryDirectory() as tmp:
            pf, tf, sf, pre = (os.path.join(tmp, x) for x in ("params.txt", "table.bin", "sample.bin", "out"))
            _write_params(pf, self.ref.R, self.ref.templates, order, ambiguous, perfect)
            _write_table(tables, tf)
            with open(sf, "wb") as f:
                np.array([n], np.int32).tofile(f)
                np.asarray(self.L, np.float64).tofile(f)
                np.array([1, 1, 1], np.int32).tofile(f)
                np.ascontiguousarray(self.pts % self.L, np.float64).tofile(f)
            env = dict(os.environ, OMP_NUM_THREADS=str(threads or voro.default_threads()))
            run = subprocess.run([VB_CLOSEPACKED, pf, tf, sf, pre], capture_output=True, text=True, env=env, check=True)
            self.stats, self.engine_log = run.stdout.strip(), run.stderr
            self.status = np.fromfile(pre + ".cp", dtype=np.int8)
            vb = np.fromfile(pre + ".vnbrs", dtype=np.int32)
            self._load(pre, n)
        # face-sharing (Voronoi) neighbours of every atom, which define the clusters and their shells; cells[i].nbr are
        # the labelled first neighbours (which leave out pairs across a split octahedral hole)
        self.vcells = _Cells(vb[:n + 1].astype(np.int64), vb[n + 1:].astype(np.int64))
        self.certified = self.kept.copy()           # atoms whose first neighbours are not decided join the clusters

    def burgers(self, loop):
        b = super().burgers(loop)
        if b is not None and isinstance(self.ref, FCCReference):
            b["b_lab"] = b["b_lab"] * (self.ref.a / 2)  # fcc ideal vectors are in units of a/2
        return b
