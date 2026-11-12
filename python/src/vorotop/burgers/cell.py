"""Combinatorial Voronoi cells: half-edge structure, global vertex keys, isomorphisms."""

from collections import defaultdict


class Cell:
    """One Voronoi cell as a planar map.

    nbr[f]  : the particle across face f (negative: container wall)
    cyc[f]  : the face's vertex cycle (local vertex ids), oriented consistently over the cell
    vkey[v] : the global identity of vertex v, the sorted tuple of the cells meeting there; this
              is what lets two cells sharing a face agree on its vertices combinatorially
    simple  : every vertex lies on exactly 3 faces (4 cells meet there)
    """

    def __init__(self, i, faces):
        self.i = i
        self.nbr = [nb for nb, _ in faces]
        self.cyc = [list(c) for _, c in faces]
        inc = defaultdict(list)
        for f, c in enumerate(self.cyc):
            for v in c:
                inc[v].append(f)
        self.vkey = {v: tuple(sorted({i, *[self.nbr[f] for f in fs]})) for v, fs in inc.items()}
        self.simple = all(len(fs) == 3 for fs in inc.values())
        self.he = {}
        for f, c in enumerate(self.cyc):
            m = len(c)
            for k in range(m):
                self.he[(c[k], c[(k + 1) % m])] = f
        self._sig = None

    def next(self, u, v):
        c = self.cyc[self.he[(u, v)]]
        return v, c[(c.index(v) + 1) % len(c)]

    def signature(self):
        """Cheap invariant (sorted face sizes) used to shortlist candidate references."""
        if self._sig is None:
            self._sig = tuple(sorted(len(c) for c in self.cyc))
        return self._sig

    def mirror(self):
        """The mirror-image cell: every face cycle reversed."""
        return Cell(self.i, [(nb, c[::-1]) for nb, c in zip(self.nbr, self.cyc)])


def isomorphism(cell, ref, start=0, flab=None, Flab=None):
    """Orientation-preserving isomorphism cell -> ref that sends cell's first half-edge to ref's
    half-edge number `start` (iterating `start` over all half-edges gives all isomorphisms).
    Optional face labels flab[f] / Flab[F] must be preserved. Returns (vmap, fmap) or None."""
    h0 = next(iter(cell.he))
    H0 = list(ref.he)[start]
    vmap, fmap = {}, {}
    stack = [(h0, H0)]
    while stack:
        (u, v), (U, V) = stack.pop()
        if (u, v) not in cell.he or (U, V) not in ref.he:
            return None
        for a, A in ((u, U), (v, V)):
            if vmap.setdefault(a, A) != A:
                return None
        f, F = cell.he[(u, v)], ref.he[(U, V)]
        if f in fmap:
            if fmap[f] != F:
                return None
            continue
        if len(cell.cyc[f]) != len(ref.cyc[F]):
            return None
        if flab is not None and flab[f] != Flab[F]:
            return None
        fmap[f] = F
        h, Hh = (u, v), (U, V)
        for _ in range(len(cell.cyc[f])):
            stack.append(((h[1], h[0]), (Hh[1], Hh[0])))
            h, Hh = cell.next(*h), ref.next(*Hh)
            if vmap.setdefault(h[0], Hh[0]) != Hh[0]:
                return None
    if len(fmap) != len(cell.cyc) or len(fmap) != len(ref.cyc):
        return None
    return vmap, fmap


def all_isomorphisms(cell, ref, flab=None, Flab=None, first_only=False):
    out = []
    if cell.signature() != ref.signature():
        return out
    for s in range(len(ref.he)):
        m = isomorphism(cell, ref, s, flab, Flab)
        if m is not None:
            out.append(m)
            if first_only:
                break
    return out
