"""The analysis as described in the paper, for crystals with topologically stable Voronoi cells: label every cell by
its topology against the ideal crystal, glue neighboring crystalline cells by exact frame relations, group the other atoms
into defect clusters, and read each cluster's Burgers vectors from its boundary surface.

    crystal = vorotop.burgers.api.crystal_from_cif("theta.cif")          # snapped onto its space group
    result = vorotop.burgers.api.analyze(points, box, crystal, species)
    for c in result["clusters"]: print(c["n_atoms"], c["vectors"])

Everything is read from the atoms and the periodic box; there is no window and nothing to tune."""
import re
import warnings

import numpy as np
import spglib

from . import crystal as _crystal, shell, voro
from .engine import FastConnection

# tolerances (A) tried, smallest first, when snapping a structure onto its space group
SNAP_TOLS = (1e-5, 1e-4, 1e-3, 3e-3, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2)
SMALL_GAP = 0.2            # below this relative Voronoi gap, thermal noise and strain blur the topology (paper, SI)


class Crystal:
    """The ideal crystal: its reference cells, and how the input was snapped onto its space group."""

    def __init__(self, lattice, frac, numbers, name=None, spacegroup=None, snap=True):
        lattice, frac, numbers = np.asarray(lattice, float), np.asarray(frac, float) % 1.0, np.asarray(numbers, int)
        self.snap_note = "not snapped"
        self.cif_lattice, self.cif_frac, self.cif_numbers = lattice, frac, numbers   # the input cell, before snapping
        if snap:
            lattice, frac, numbers, self.snap_note = _snap(lattice, frac, numbers, spacegroup)
        self.kind, self.cp = _close_packed(lattice, frac, numbers)
        if self.kind != "stable":                     # fcc or hcp: degenerate cells, read with the shipped tables
            self.reference = None
            self.lattice, self.frac, self.numbers = lattice, frac, numbers
            self.gap = 0.0
            self.r_max = self.cp["r_max"]
            ds = spglib.get_symmetry_dataset((lattice, frac, numbers), symprec=1e-3)
            self.spacegroup = ds.international if ds is not None else "?"
            self._to_conventional = np.eye(3)
            self._cell_symmetry(ds)
            return
        self.reference = _crystal.Reference(lattice, frac, numbers, name=name)
        self.lattice, self.frac, self.numbers = self.reference.A, self.reference.frac, self.reference.numbers
        self.gap = float(self.reference.min_gap())
        if self.gap <= 0:                              # degenerate cells other than fcc and hcp (e.g. omega-Ti): every
            from . import degenerate                   # perturbed topology a reference cell of its own (VoroTypes)
            self.kind = "degenerate"
            self.reference = degenerate.reference(self.lattice, self.frac, self.numbers, name=name)
        d_of = lambda r: np.array(list(r["d"].values()) if isinstance(r["d"], dict) else r["d"], float)
        self.r_max = max(float(np.linalg.norm(d_of(r), axis=1).max()) for r in self.reference.refs)   # longest neighbor
        ds = spglib.get_symmetry_dataset((self.lattice, self.frac, self.numbers), symprec=1e-3)
        self.spacegroup = ds.international if ds is not None else "?"
        self._to_conventional = np.asarray(ds.transformation_matrix) if ds is not None else np.eye(3)
        self._cell_symmetry(ds)

    def _cell_symmetry(self, ds):
        """The point group of the input cell (rotations in its lattice coordinates) and the map from conventional-cell
        to input-cell coordinates, from spglib at the smallest tolerance that finds the snapped space group."""
        self._cell_rotations, self._from_conventional, self._conv_rotations = None, None, None
        if ds is None:
            return
        std = spglib.get_symmetry_dataset((ds.std_lattice, ds.std_positions, ds.std_types), symprec=1e-3)
        if std is not None:
            self._conv_rotations = np.unique(np.asarray(std.rotations), axis=0)
        cell = (self.cif_lattice, self.cif_frac, self.cif_numbers)
        for tol in SNAP_TOLS:
            d0 = spglib.get_symmetry_dataset(cell, symprec=tol)
            if d0 is not None and d0.number == ds.number:
                self._cell_rotations = np.unique(np.asarray(d0.rotations), axis=0)
                self._from_conventional = np.linalg.inv(np.asarray(d0.transformation_matrix))
                return

    def conventional(self, n):
        """Lattice coordinates n (of self.lattice) in the conventional cell of the space group, in spglib's setting:
        the same as the input cell's up to a symmetry of the crystal."""
        return self._to_conventional @ np.asarray(n, float)

    def family(self, conventional):
        """The family of a vector given in the conventional cell, e.g. "1/2<110>": its symmetry image with the fewest
        negative and nonzero entries, largest first."""
        x = np.asarray(conventional, float)
        if self._conv_rotations is not None:
            images = np.round(np.einsum("kij,j->ki", self._conv_rotations, x), 6) + 0.0
            x = min(images, key=lambda y: ((y < 0).sum(), (y != 0).sum(), tuple(-y)))
        from fractions import Fraction
        f = [Fraction(float(v)).limit_denominator(24) for v in x]
        den = int(np.lcm.reduce([q.denominator for q in f]))
        num = "".join(str(int(q * den)) if q >= 0 else "-" + str(int(-q * den)) for q in f)
        sep = "," if any(abs(int(q * den)) > 9 for q in f) else ""
        if sep:
            num = ",".join(str(int(q * den)) for q in f)
        return (f"1/{den}" if den > 1 else "") + f"<{num}>"

    def axes(self, orient):
        """Rows: the unit vectors, in the input cell's Cartesian frame, of the sample's x, y and z axes. orient: three
        directions in lattice indices of the input cell (as LAMMPS orient x, y, z), e.g. [[1,1,-2],[1,1,1],[1,-1,0]]."""
        E = np.asarray(orient, float) @ self.cif_lattice
        E /= np.linalg.norm(E, axis=1)[:, None]
        if np.abs(E @ E.T - np.eye(3)).max() > 1e-3:
            raise ValueError("the orientation directions must be mutually perpendicular")
        if np.linalg.det(E) < 0:
            raise ValueError("the orientation directions must form a right-handed set")
        return E

    def in_cell(self, conventional, b_lab, E):
        """Exact coordinates, in the input cell, of a Burgers vector given in the conventional cell (up to a symmetry
        of the crystal): the symmetry image closest to the laboratory vector b_lab of a sample with axes E (see axes).
        None if no image is close, which means the orientation given does not describe the sample."""
        measured = np.linalg.solve(self.cif_lattice.T, E.T @ np.asarray(b_lab, float))
        if self._cell_rotations is None:
            return None
        x0 = self._from_conventional @ np.asarray(conventional, float)
        images = np.einsum("kij,j->ki", self._cell_rotations, x0)
        err = np.linalg.norm((images - measured) @ self.cif_lattice, axis=1)
        k = int(np.argmin(err))
        return images[k] if err[k] < 0.05 * np.linalg.norm(x0 @ self.cif_lattice) else None


HCP_TABLES = {"hcp": 1.6237, "hcp_ideal": float(np.sqrt(8 / 3))}     # shipped hcp tables and their c/a


def _close_packed(lattice, frac, numbers):
    """("fcc", params) or ("hcp", params) for a one-species fcc or hcp crystal, else ("stable", None)."""
    if len(set(np.asarray(numbers).tolist())) != 1:
        return "stable", None
    ds = spglib.get_symmetry_dataset((lattice, frac, numbers), symprec=1e-3)
    if ds is None:
        return "stable", None
    L = np.asarray(ds.std_lattice)
    if ds.number == 225 and len(ds.std_positions) == 4:
        a = float(np.linalg.norm(L[0]))
        return "fcc", dict(a=a, r_max=a / np.sqrt(2) * 1.0001, table="fcc")
    if ds.number == 194 and len(ds.std_positions) == 2:
        a, c = float(np.linalg.norm(L[0])), float(np.linalg.norm(L[2]))
        table = min(HCP_TABLES, key=lambda k: abs(HCP_TABLES[k] - c / a))
        return "hcp", dict(a=a, c=c, r_max=max(a, np.sqrt(a * a / 3 + c * c / 4)), table=table,
                           table_ratio=HCP_TABLES[table])
    return "stable", None


def _snap(lattice, frac, numbers, spacegroup):
    """Exact positions of the space group (spglib standardize_cell with idealization, primitive cell), at the smallest
    tolerance at which spglib finds it. With no space group given, the one found at 1e-3 A is used."""
    cell = (lattice, frac, numbers)
    if spacegroup is None:
        ds = spglib.get_symmetry_dataset(cell, symprec=1e-3)
        if ds is None:
            return lattice, frac, numbers, "no space group found; used as given"
        spacegroup = ds.number
    for tol in SNAP_TOLS:
        ds = spglib.get_symmetry_dataset(cell, symprec=tol)
        if ds is not None and ds.number == int(spacegroup):
            L, F, N = spglib.standardize_cell(cell, to_primitive=True, no_idealize=False, symprec=tol)
            return np.asarray(L), np.asarray(F), np.asarray(N), f"snapped onto space group {spacegroup} at {tol:g} A"
    return lattice, frac, numbers, f"space group {spacegroup} not found up to {SNAP_TOLS[-1]} A; used as given"


def crystal_from_cif(path, spacegroup=None, snap=True):
    """The ideal crystal from a CIF file (needs pymatgen). spacegroup: the declared space-group number, if known (for
    example from the Materials Project); otherwise the one spglib finds."""
    from pymatgen.core import Structure
    s = Structure.from_file(path)
    return Crystal(s.lattice.matrix, s.frac_coords, [site.specie.Z for site in s], name=str(path),
                   spacegroup=spacegroup, snap=snap)


def _heights(box):
    """Distances between opposite faces of the periodic box."""
    H = np.diag(box) if np.ndim(box) == 1 else np.asarray(box, float)
    V = abs(np.linalg.det(H))
    return np.array([V / np.linalg.norm(np.cross(H[(k + 1) % 3], H[(k + 2) % 3])) for k in range(3)])


def analyze(points, box, crystal, species=None, threads=None, orient=None):
    """Burgers content of every defect cluster of a fully periodic sample.

    points: (N, 3) positions; box: edge lengths (orthogonal box with its corner at the origin) or a 3x3 restricted-
    triclinic matrix (rows a = (lx,0,0), b = (xy,ly,0), c = (xz,yz,lz)); species: atomic numbers (needed when the
    crystal has more than one species); orient: the crystal directions, in lattice indices of the input cell, along the
    sample's x, y and z axes (as LAMMPS orient), if known. Every vector is given in the laboratory frame (exact, up to
    sign), in the reference lattice and in the conventional cell (both up to a symmetry of the crystal: without orient,
    the analysis cannot tell which of the symmetry-equivalent orientations the sample has), and, when orient is
    given, in the input cell ("cell", exact up to sign). Returns a dict with the clusters, per-atom arrays and
    warnings. Per-atom "branch": the Burgers vector (lattice coordinates, up to sign) of the branch each atom of a
    dislocation cluster lies on (shell.branch_labels), NaN elsewhere; "branch_ties": atoms equally near several
    branches, with all their vectors."""
    points = np.asarray(points, float)
    if orient is not None:
        return _with_cell_coords(analyze(points, box, crystal, species, threads), crystal, crystal.axes(orient))
    if crystal.kind in ("fcc", "hcp"):
        return _analyze_close_packed(points, box, crystal, threads)
    ref = crystal.reference
    multi = len(set(ref.numbers.tolist())) > 1
    if multi and species is None:
        raise ValueError("the crystal has several species: pass species (atomic numbers) for every atom")
    tri = np.ndim(box) == 2
    H = np.asarray(box, float)
    if tri:
        frac = points @ np.linalg.inv(H)
        points = (frac % 1.0) @ H
    else:
        points = points % H
    warn = []
    if crystal.kind == "stable" and crystal.gap < SMALL_GAP:
        warn.append(f"small Voronoi gap g = {crystal.gap:.3g} (< {SMALL_GAP}): thermal noise and strain may destroy the "
                    "ideal cell topology; expect wide cores, merged or missing lines (never wrong vectors)")
    h = _heights(H)
    if (h < 3 * crystal.r_max).any():
        warn.append(f"box height {h.min():.2f} A is less than three times the longest neighbor distance "
                    f"({crystal.r_max:.2f} A), so a smallest loop can wrap around the box: enlarge the sample")
    F = FastConnection(points, H, (1, 1, 1), ref, species=np.asarray(species, int) if multi else None,
                       threads=threads, orientation=np.eye(3))
    stats = {k: float(v) for k, v in re.findall(r"(\w+) (-?[\d.]+(?:e-?\d+)?)", F.stats)}
    if stats.get("forced", 0) > 0:
        warn.append(f"{int(stats['forced'])} cells were labeled by geometric fit where topology could not decide")
    nt, no, worst = F.loop_closure()
    stats["smallest_loops"], stats["open_loops"] = nt, no
    if no:
        warn.append(f"{no} of {nt} smallest loops do not close (largest net translation of those whose rotation "
                    f"closes: {worst:.3g} A)")
    if stats.get("pseudo", 0) > 0 and stats.get("margin", 1.0) < 0.01:
        warn.append(f"pseudo-symmetric crystal with a near tie between labelings (misfit margin {stats['margin']:.3g} A)")
    R = shell.analyse(F, faces=None, mask=None)             # faces of the defect cells, from the engine
    branch, branch_sets = shell.branch_labels(F, R, F.defect_faces)
    clusters, cluster_of = [], np.full(len(points), -1)
    for k, c in enumerate(R):
        cells = np.asarray(c["cells"], int)
        cluster_of[cells] = k
        vecs = []
        for pc in c["pieces"]:
            for v in pc.get("vectors", []):
                rid = pc["region"]
                lab = F.Qlab[rid] @ np.linalg.solve(F.G[rid], np.asarray(v, float) @ ref.A)
                vecs.append(dict(lattice=np.round(np.asarray(v, float), 6).tolist(),
                                 conventional=np.round(crystal.conventional(v), 6).tolist(),
                                 family=crystal.family(crystal.conventional(v)),
                                 lab=np.round(lab, 6).tolist(), length=round(float(np.linalg.norm(lab)), 6),
                                 region=int(rid)))
        clusters.append(dict(id=k, n_atoms=int(len(cells)), genus=[pc["genus"] for pc in c["pieces"]], vectors=vecs,
                             atoms=cells))
    for w in warn:
        warnings.warn(w)
    return dict(clusters=clusters, crystalline=np.asarray(F.kept, bool), cluster=cluster_of, branch=branch,
                branch_ties={a: [np.asarray(v).tolist() for v in l] for a, l in branch_sets.items() if len(l) > 1},
                crystalline_fraction=float(np.mean(F.kept)), regions=int(stats.get("regions", 0)), engine=stats,
                warnings=warn, crystal=dict(spacegroup=crystal.spacegroup, snap=crystal.snap_note,
                                            gap=crystal.gap))


def _with_cell_coords(res, crystal, E):
    """Adds each vector's exact coordinates in the input cell, the symmetry image singled out by the sample axes E."""
    missed = 0
    for c in res["clusters"]:
        for v in c["vectors"]:
            x = crystal.in_cell(v["conventional"], v["lab"], E)
            v["cell"] = None if x is None else np.round(x, 6).tolist()
            missed += x is None
    if missed:
        w = f"{missed} vectors match no symmetry image in the orientation given: check orient"
        res["warnings"].append(w)
        warnings.warn(w)
    return res


def _analyze_close_packed(points, box, crystal, threads=None):
    """fcc and hcp metals: first neighbors from the table of perturbed cell topologies and neighbor agreement, then the
    same boundary-surface analysis. Orthogonal boxes only."""
    if np.ndim(box) == 2:
        H = np.asarray(box, float)
        if np.abs(H[np.tril_indices(3, -1)]).max() > 1e-9:
            raise ValueError("the close-packed analysis currently needs an orthogonal box")
        box = np.diag(H)
    L = np.asarray(box, float)
    points = np.asarray(points, float) % L
    cp, warn = crystal.cp, []
    h = _heights(L)
    if (h < 3 * crystal.r_max).any():
        warn.append(f"box height {h.min():.2f} A is less than three times the longest neighbor distance "
                    f"({crystal.r_max:.2f} A), so a smallest loop can wrap around the box: enlarge the sample")
    from .closepacked import CPConnection
    if crystal.kind == "fcc":
        F = CPConnection(points, L, "fcc", cp["a"], threads=threads)
        A = cp["a"] * np.eye(3)                       # vectors are in units of the cubic cell
    else:
        if abs(cp["c"] / cp["a"] - cp["table_ratio"]) > 0.005 * cp["table_ratio"]:
            warn.append(f"c/a = {cp['c'] / cp['a']:.4f}; the nearest shipped table is for {cp['table_ratio']:.4f}: "
                        "build a table for this ratio (VoroTypes) for exact first neighbors")
        F = CPConnection(points, L, "hcp", cp["a"], cp["c"], threads=threads, tables=(cp["table"],))
        A = F.ref.A
    R = shell.analyse(F, faces=None, mask=None)             # faces of the defect cells, from the tool
    several = sum(pc.get("regions", 0) > 1 for c in R for pc in c["pieces"])
    if several:
        warn.append(f"{several} boundary pieces face more than one crystalline region; their curves were not combined")
    branch, branch_sets = shell.branch_labels(F, R, F.defect_faces)
    clusters, cluster_of = [], np.full(len(points), -1)
    for k, c in enumerate(R):
        cells = np.asarray(c["cells"], int)
        cluster_of[cells] = k
        vecs = []
        for pc in c["pieces"]:
            for v in pc.get("vectors", []):
                rid = pc["region"]
                lab = F.Qlab[rid] @ np.linalg.solve(F.G[rid], np.asarray(v, float) @ A)
                vecs.append(dict(lattice=np.round(np.asarray(v, float), 6).tolist(),
                                 conventional=np.round(np.asarray(v, float), 6).tolist(),
                                 family=crystal.family(v),
                                 lab=np.round(lab, 6).tolist(), length=round(float(np.linalg.norm(lab)), 6),
                                 region=int(rid)))
        clusters.append(dict(id=k, n_atoms=int(len(cells)), genus=[pc["genus"] for pc in c["pieces"]], vectors=vecs,
                             atoms=cells))
    for w in warn:
        warnings.warn(w)
    return dict(clusters=clusters, crystalline=np.asarray(F.kept, bool), cluster=cluster_of, branch=branch,
                branch_ties={a: [np.asarray(v).tolist() for v in l] for a, l in branch_sets.items() if len(l) > 1},
                crystalline_fraction=float(np.mean(F.kept)), regions=len(getattr(F, "G", {})), engine={},
                warnings=warn, crystal=dict(spacegroup=crystal.spacegroup, snap=crystal.snap_note, gap=0.0,
                                            route=crystal.kind))
