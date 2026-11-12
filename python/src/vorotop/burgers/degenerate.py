"""Crystals with degenerate Voronoi cells, by enumeration: every cell topology that a small perturbation can produce
becomes a reference cell of its own.

A perturbation small enough to leave the generic Voronoi vertices of the ideal crystal intact changes a cell only at its
degenerate vertices, each of which resolves in finitely many ways; VoroTypes enumerates all compatible combinations up
to the site symmetry. Each resolved cell keeps, for every face, the atom of the ideal crystal across it, so it carries
an exact labeling: the ideal vector to that atom. Its vertices are the circumcenters of its Delaunay tetrahedra in the
ideal crystal (the vertices created by one degenerate vertex coincide there). Registered as reference cells of the
crystal, these are matched, glued and transported by the engine exactly as the cells of a stable crystal; a cell
whose topology has several labelings is settled by arc consistency on the gluings. Neighbor labels are wildcards (the
neighbors of a perturbed cell have perturbed topologies of their own).

Used for omega-Ti, whose degenerate vertices never add or remove a face. Needs VoroTypes (voronoi_enumerate) on the
path ($VOROTYPES, if VoroTypes is not installed)."""
import os
import sys
from itertools import combinations

import numpy as np

from .cell import Cell, all_isomorphisms
from .crystal import Reference


def _need_vorotypes():
    """The enumeration of perturbed cells comes from VoroTypes (its voronoi_enumerate package), an optional
    dependency needed only for crystals with degenerate cells other than fcc and hcp; $VOROTYPES may point to a
    source checkout."""
    if os.environ.get("VOROTYPES"):
        sys.path.insert(0, os.environ["VOROTYPES"])
    try:
        import voronoi_enumerate  # noqa: F401
    except ImportError as e:
        raise ImportError("crystals with degenerate cells need VoroTypes (voronoi_enumerate); install it or set "
                          "$VOROTYPES to its source directory") from e


def _enumerate(lattice, frac, site):
    _need_vorotypes()
    from voronoi_enumerate.cell import orient_faces, star_to_faces
    from voronoi_enumerate.combine import _make_prune_fn, _precompute_compatibility
    from voronoi_enumerate.crystal import Crystal
    from voronoi_enumerate.resolution import enumerate_resolutions
    from voronoi_enumerate.symmetry import compute_site_symmetry
    from voronoi_enumerate.voronoi import analyze_voronoi
    cryst = Crystal(lattice, frac, ["A"] * len(frac))
    vertices, c0, coords, _ = analyze_voronoi(cryst, atom_index=site, n_images=3, near_gap_threshold=None)
    generic = frozenset(tuple(sorted(v.atom_indices)) for v in vertices if not v.is_degenerate)
    degen = [v for v in vertices if v.is_degenerate]
    vo = []
    for v in degen:
        order = [c0] + [a for a in v.atom_indices if a != c0]
        res, _ = enumerate_resolutions(coords[order], central=0)
        vo.append([(frozenset(tuple(sorted(order[j] for j in t)) for t in r.star), False) for r in res])
    sym = compute_site_symmetry(cryst, c0, coords, degen, [[t for t, u in o if not u] for o in vo], verbose=False)
    adj, compat, hc = _precompute_compatibility(generic, vo, c0)
    prune = _make_prune_fn(adj, compat) if hc else None
    # each tetrahedron's vertex: a generic vertex is its own; a tetrahedron of a resolution belongs to the degenerate
    # vertex that was resolved (two degenerate vertices can share four atoms)
    gpos = {tuple(sorted(v.atom_indices)): np.asarray(v.position, float) for v in vertices if not v.is_degenerate}
    dpos = [np.asarray(v.position, float) for v in degen]
    out = []
    for combo, _ in sym.orderly_generate([len(o) for o in vo], prune=prune):
        pos = dict(gpos)
        for vi, ri in enumerate(combo):
            for t in vo[vi][ri][0]:
                pos[t] = dpos[vi]
        star = sorted(pos)
        if not _closed(star, c0):
            continue
        faces, fn = star_to_faces(star, central=c0)
        if faces:
            out.append((star, orient_faces(faces), fn, [pos[t] for t in star]))
    return out, c0, coords, len(degen)


def _closed(star, c0):
    """Is the star a closed triangulated ball around c0: every triangle through c0 in exactly two tetrahedra, and no two
    tetrahedra the same? (Combinations of resolutions that are only pairwise compatible can fail this.)"""
    from collections import Counter
    if len(set(star)) != len(star):
        return False
    cnt = Counter(tuple(sorted(f)) for t in star for f in combinations(t, 3) if c0 in f)
    return all(n == 2 for n in cnt.values())


def _same_up_to(D1, D2, R, tol):
    """Is the labeling D2 (face -> vector) a point-group image of D1?"""
    return any(np.abs(D1 @ Rm.T - D2).max() < tol for Rm in R)


def reference(lattice, frac, numbers, **kw):
    """A crystal.Reference whose reference cells are the resolved cells of each site: one per topology and distinct
    labeling (up to the point group). Combinations of resolutions that do not close into a cell are skipped."""
    _need_vorotypes()
    from voronoi_enumerate.weinberg import weinberg_vector
    ref = Reference(lattice, frac, numbers, **kw)
    A = np.asarray(lattice, float)
    tol = 1e-4 * ref.scale
    types = {}                                   # Weinberg vector -> (reference cell, [(site, D, P)])
    for s in sorted(set(ref.orbit_of_site)):
        cells, c0, coords, ndeg = _enumerate(A, np.asarray(frac, float), int(s))
        for star, faces, fn, tpos in cells:
            w = weinberg_vector(faces)
            if w is None:
                continue
            c = Cell(c0, list(zip(fn, faces)))
            if not c.simple:
                continue
            d = {f: coords[nb] - coords[c0] for f, nb in enumerate(c.nbr)}
            P = {v: tpos[v] - coords[c0] for v in c.vkey}
            w = tuple(w)
            if w not in types:
                types[w] = (c, [(int(s), np.array([d[f] for f in range(len(c.nbr))]), P)])
                continue
            c_ref, labs = types[w]
            for tgt in (c_ref, c_ref.mirror()):
                for vmap, fmap in all_isomorphisms(c, tgt):
                    D = np.zeros((len(c.nbr), 3))
                    for f, F in fmap.items():
                        D[F] = d[f]
                    if not any(si == s and _same_up_to(D0, D, ref.R, tol) for si, D0, _ in labs):
                        labs.append((int(s), D, {V: P[v] for v, V in vmap.items()}))
    refs = []
    for w, (c, labs) in types.items():
        mirror = c.mirror()
        n_comb = len(all_isomorphisms(c, c)) + len(all_isomorphisms(c, mirror))
        for s, D, P in labs:
            n_site = sum(_same_up_to(D, D[perm], [Rm], tol) for Rm in ref.R
                         for perm in [np.linalg.norm((D @ Rm.T)[:, None] - D[None], axis=2).argmin(1)])
            refs.append(dict(site=s, species=int(numbers[s]), cell=c, mirror=mirror,
                             d={f: D[f] for f in range(len(D))}, P=P, signature=c.signature(), simple=True,
                             Flab=[(None, int(numbers[s]))] * len(D), n_comb=n_comb, n_comb_lab=n_comb, n_site=n_site,
                             rigid=(n_comb == n_site and len(labs) == 1), rigid_lab=False, weinberg=w,
                             glue_by_neighbors=True))
    ref.refs = refs
    return ref


def omega(a, c, Z=22):
    """omega phase (P6/mmm): 1a (0,0,0) and 2d (1/3,2/3,1/2), (2/3,1/3,1/2)."""
    A = np.array([[a, 0, 0], [-0.5 * a, np.sqrt(3) / 2 * a, 0], [0, 0, c]])
    return reference(A, [[0, 0, 0], [1 / 3, 2 / 3, 0.5], [2 / 3, 1 / 3, 0.5]], [Z] * 3)
