"""Defect clusters, their boundary surfaces, and the Burgers content of each cluster from a basis of loops on
its boundary (theory note: "one layer suffices", "canonical loops").

A defect cluster D is a connected component (through Voronoi faces) of the atoms that are not kept. Its boundary
Sigma is the closed surface made of the Voronoi faces between D and the other atoms. Everything about Sigma is
read from groups of atoms, with no geometry:
    a face of Sigma   <->  a mixed pair      {d, x}   (d in D, x not in D)
    an edge of Sigma  <->  a mixed triple    (a Voronoi edge with atoms both in D and not)
    a vertex of Sigma <->  a mixed quadruple (a Voronoi vertex likewise)
so chi = #quadruples - #triples + #pairs = 2 - 2g on each connected piece. A tree-cotree decomposition of the dual
cell structure (faces as nodes, edges as links, vertices as polygons) gives 2g independent closed curves on each
piece. Each becomes a closed walk through the first layer (the atoms x, consecutive ones sharing a face), whose
holonomy is the Burgers content of that curve. Spheres (g = 0) carry none. Integer combinations of the basis
curves that do not wrap the periodic box give the Burgers vectors of the cluster's dislocations.

Requires generic Voronoi geometry (three cells per edge, four per vertex), as in any thermal snapshot, and a
connection whose kept cells glue across faces (FastConnection or crystal.Connection: crystals with stable cells).
"""
import os
from collections import deque, defaultdict
from fractions import Fraction

import numpy as np

from . import voro


def _structure(faces, members):
    """Pairs, triples (Voronoi edges) and quadruples (Voronoi vertices) at the cells `members`, as sorted tuples,
    from each cell's faces: an edge of a cell lies on two of its faces, a vertex on three. Also returns the cells
    whose structure is missing or not generic."""
    pairs, triples, quads, bad = set(), set(), set(), set()
    for i in members:
        fs = faces[i]
        if fs is None:
            bad.add(i); continue
        at_vertex, at_edge = defaultdict(list), defaultdict(list)
        for nb, cyc in fs:
            pairs.add((min(i, nb), max(i, nb)))
            for m in range(len(cyc)):
                u, v = cyc[m], cyc[(m + 1) % len(cyc)]
                at_vertex[u].append(nb)
                at_edge[(min(u, v), max(u, v))].append(nb)
        for nbs in at_edge.values():
            if len(nbs) == 2:
                triples.add(tuple(sorted((i, *nbs))))
            else:
                bad.add(i)
        for nbs in at_vertex.values():
            if len(nbs) == 3:
                quads.add(tuple(sorted((i, *nbs))))
            else:
                bad.add(i)
    return pairs, triples, quads, bad


def _good(conn):
    """Atoms whose cells are kept and, for degenerate crystals analysed with certification, whose first
    neighbours are uniquely determined by Voronoi topology."""
    return np.asarray(conn.kept, bool) & np.asarray(getattr(conn, "certified", np.ones(len(conn.pts), bool)), bool)


def clusters(conn, mask=None, faces=None):
    """Connected components (through Voronoi faces) of the atoms that are not kept, within `mask`. Face adjacency
    is taken from `faces` (the raw tessellation) when given, else from the connection. Wall cells (cells with a face
    on a container wall, in a non-periodic direction) are never defects: they belong to the boundary of the sample,
    and including them would merge every defect near a free surface with the whole surface layer."""
    n = len(conn.pts)
    bad = ~_good(conn)
    if hasattr(faces, "structure"):               # the engine's faces: use its face-sharing neighbor lists
        faces = None
    vc = getattr(conn, "vcells", conn.cells)       # face-sharing neighbors (= cells for the stable engine)
    nbr = lambda a: vc[a].nbr if faces is None or faces[a] is None else [nb for nb, _ in faces[a]]
    wall = np.array([any(b < 0 for b in nbr(a)) for a in range(n)]) if not all(conn.periodic) else np.zeros(n, bool)
    bad &= ~wall
    if mask is not None:
        bad &= np.asarray(mask, bool)
    seen = np.zeros(n, bool)
    out = []
    for s in np.flatnonzero(bad):
        if seen[s]:
            continue
        comp, dq = [], deque([s]); seen[s] = True
        while dq:
            a = dq.popleft(); comp.append(a)
            for b in nbr(a):
                if b >= 0 and bad[b] and not seen[b]:
                    seen[b] = True; dq.append(b)
        out.append(np.array(sorted(comp)))
    return out


def _integer_kernel(W):
    """Integer basis of {c : c @ W = 0} for an integer matrix W (m x 3), by exact Gaussian elimination."""
    W = [[Fraction(int(v)) for v in row] for row in np.asarray(W).T]       # 3 x m: solve W^T c = 0
    m = len(W[0]) if W else 0
    piv, r = [], 0
    for col in range(m):
        p = next((i for i in range(r, len(W)) if W[i][col] != 0), None)
        if p is None:
            continue
        W[r], W[p] = W[p], W[r]
        W[r] = [v / W[r][col] for v in W[r]]
        for i in range(len(W)):
            if i != r and W[i][col] != 0:
                W[i] = [a - W[i][col] * b for a, b in zip(W[i], W[r])]
        piv.append(col); r += 1
    free = [c for c in range(m) if c not in piv]
    basis = []
    for f in free:
        c = [Fraction(0)] * m; c[f] = Fraction(1)
        for i, pc in enumerate(piv):
            c[pc] = -W[i][f]
        den = np.lcm.reduce([x.denominator for x in c])
        basis.append(np.array([int(x * den) for x in c]))
    return basis


def _span(vs, den=12):
    """Generators of the lattice spanned by the vectors vs (rational, denominators dividing `den`), in Hermite normal
    form (unique for each lattice): integer row reduction (Euclid, column by column) of den * vs, pivots positive, and
    the entries above each pivot reduced into [0, pivot)."""
    M = [np.rint(np.asarray(v, float) * den).astype(np.int64) for v in vs]
    M = [r for r in M if r.any()]
    piv_rows, piv_cols = [], []
    for col in range(3):
        while True:
            nz = [k for k, r in enumerate(M) if r[col] != 0]
            if len(nz) <= 1:
                break
            p = min(nz, key=lambda k: abs(M[k][col]))
            for k in nz:
                if k != p:
                    M[k] = M[k] - (M[k][col] // M[p][col]) * M[p]
        nz = [k for k, r in enumerate(M) if r[col] != 0]
        if nz:
            piv = M.pop(nz[0])
            piv_rows.append(piv if piv[col] > 0 else -piv); piv_cols.append(col)
        M = [r for r in M if r.any()]
    for i in range(len(piv_rows)):
        for j in range(i):
            c = piv_cols[i]
            piv_rows[j] = piv_rows[j] - (piv_rows[j][c] // piv_rows[i][c]) * piv_rows[i]
    return [r / den for r in piv_rows]


def _nonwrapping_span(W, B, den=12):
    """Generators (Hermite normal form) of the lattice of net translations of the integer combinations of curves whose
    total winding is zero: rows (W | den * B) reduced by integer row operations (Euclid) on the winding columns; the
    rows left with zero winding generate exactly the combinations with zero winding (the integer kernel, not a
    rational basis of it scaled to integers, which can generate a sublattice)."""
    M = np.c_[np.asarray(W, np.int64), np.rint(np.asarray(B, float) * den).astype(np.int64)]
    rows = [r for r in M]
    for col in range(M.shape[1] - 3):
        while True:
            nz = [k for k, r in enumerate(rows) if r[col] != 0]
            if len(nz) <= 1:
                break
            p = min(nz, key=lambda k: abs(rows[k][col]))
            for k in nz:
                if k != p:
                    rows[k] = rows[k] - (rows[k][col] // rows[p][col]) * rows[p]
        nz = [k for k, r in enumerate(rows) if r[col] != 0]
        if nz:
            rows.pop(nz[0])
    return _span([r[-3:] / den for r in rows if r[-3:].any()], den)


def to_region(conn, b, rid_from, rid_to):
    """A Burgers vector in lattice coordinates of crystalline region rid_from, expressed in those of rid_to. Each
    region has its own crystal setting (its global frame); the vector in the lab frame is the same in both."""
    if rid_from == rid_to:
        return np.asarray(b, float)
    b_root = np.linalg.solve(conn.G[rid_from], np.asarray(b, float) @ conn.ref.A)
    b_lab = conn.Qlab[rid_from] @ b_root
    return conn.ref.lattice_coords(conn.G[rid_to] @ np.linalg.solve(conn.Qlab[rid_to], b_lab))


def in_span(b, vs, den=12):
    """Whether b is an integer combination of the vectors vs (Frank's rule for a loop around several branches):
    reduce den * b by the echelon generators of their lattice; b is in the span iff nothing remains."""
    r = np.rint(np.asarray(b, float) * den).astype(np.int64)
    for e in _span(vs, den):
        e = np.rint(e * den).astype(np.int64)
        col = int(np.flatnonzero(e)[0])
        if r[col] % e[col]:
            return False
        r = r - (r[col] // e[col]) * e
    return not r.any()


def _bonded_walk(conn, walk, inD):
    """The walk, with every step between atoms that share a Voronoi face but have no gluing (a split vertex of a
    degenerate cell for which the connection provides none) replaced by a two-step detour through a common first
    neighbour outside the cluster. (fcc and hcp connections glue the corners of a split octahedron directly.) All such detours must have the same holonomy (they differ
    by first-neighbour triangles, which close); otherwise None. For crystals with stable cells every face is a
    first-neighbour face and the walk is returned unchanged."""
    out = []
    good = _good(conn)
    W = walk + walk[:1]
    for a, b in zip(W[:-1], W[1:]):
        out.append(a)
        if b in conn.cells[a].nbr or conn.g(a, b) is not None:     # a first-neighbour step, or opposite corners of
            continue                                                # an octahedral hole (gluing from their shells)
        ys = [y for y in set(conn.cells[a].nbr) & set(conn.cells[b].nbr) if y >= 0 and good[y] and not inD[y]]
        if not ys:
            return None
        hs = []
        for y in ys:
            g1, g2 = conn.g(a, y), conn.g(y, b)
            if g1 is None or g2 is None:
                continue
            R = g1[0] @ g2[0]; t = g1[1] + g1[0] @ g2[1]
            hs.append((y, R, t))
        if not hs:
            return None
        if any(np.abs(R - hs[0][1]).max() > 1e-6 or np.abs(t - hs[0][2]).max() > 1e-6 for _, R, t in hs[1:]):
            return None
        out.append(hs[0][0])
    return out


def _regions(conn, good):
    """Connected regions of good atoms (the connection's gauge regions), as an array (-1 elsewhere)."""
    reg = getattr(conn, "region", None)
    if isinstance(reg, np.ndarray):
        out = np.where(good, reg, -1)
    else:
        out = -np.ones(len(conn.pts), np.int64)
        for i, r in conn.region_of.items():
            if good[i]:
                out[i] = r
    return out


def analyse(conn, faces=None, mask=None, min_size=1, walks=False):
    """For every defect cluster: its boundary pieces, each with Euler characteristic, genus, the basis curves
    (as first-layer walks) with their Burgers content and box wrapping, and the Burgers vectors of the
    combinations that do not wrap the box. Returns a list of dicts, one per cluster. With the engine's faces and a
    FastConnection, curves are evaluated from the integer voltages of the shell (_pieces_fast); walks=True follows each
    curve's walk instead (the reference computation, and the walks themselves, for diagnostics)."""
    if faces is None:          # the engine's own faces of the cells that are not kept, else a tessellation of all
        faces = getattr(conn, "defect_faces", None) or voro.faces(conn.pts, conn.L, conn.periodic)
    L = np.asarray(conn.L, float)
    per = np.asarray(conn.periodic, bool)
    good0 = _good(conn)
    regs = _regions(conn, good0)
    out = []
    Ds = [D for D in clusters(conn, mask, faces) if len(D) >= min_size]
    fast = not walks and hasattr(faces, "structure_all") and hasattr(conn, "g_many")
    G = shell_graph_many(conn, Ds) if fast else None
    if fast:                                     # the Voronoi faces, edges, and vertices at every cluster, at once
        label = np.full(len(conn.pts), -1, np.int64)
        for ci, D in enumerate(Ds):
            label[D] = ci
        allP, allT, allQ, lP, lT, lQ, allbad = faces.structure_all(label)
        def by_cluster(A, lab_):
            o_ = np.argsort(lab_, kind="stable"); A, lab_ = A[o_], lab_[o_]
            st = np.searchsorted(lab_, np.arange(len(Ds))); en = np.searchsorted(lab_, np.arange(len(Ds)), "right")
            return A, st, en
        allP, sP, eP = by_cluster(allP, lP); allT, sT, eT = by_cluster(allT, lT); allQ, sQ, eQ = by_cluster(allQ, lQ)
    inD = np.zeros(len(conn.pts), bool); prevD = None
    done = _spheres_batch(conn, Ds, allP, allT, allQ, sP, eP, sT, eT, sQ, eQ, allbad, good0, regs) if fast else {}
    for ci, D in enumerate(Ds):
        if prevD is not None:
            inD[prevD] = False
        inD[D] = True; prevD = D
        if ci in done:
            out.append(done[ci]); continue
        if fast:
            P3, T3, Q3 = allP[sP[ci]:eP[ci]], allT[sT[ci]:eT[ci]], allQ[sQ[ci]:eQ[ci]]
            badc = {int(x) for x in D if int(x) in allbad}
            def mixed_arr(A):
                inn = (A >= 0) & inD[np.maximum(A, 0)]
                return A[inn.any(1) & ~inn.all(1)]
            r = _pieces_fast(conn, D, inD, mixed_arr(P3), mixed_arr(T3), mixed_arr(Q3), badc, good0, regs, G, ci)
            if r is not None:
                out.append(r); continue
        isin = lambda a: a >= 0 and inD[a]
        if hasattr(faces, "structure"):            # flat arrays: keep only the mixed rows (sorted, as below)
            P3, T3, Q3, badc = faces.structure(D)
            def mixed_rows(A):
                inn = (A >= 0) & inD[np.maximum(A, 0)]
                return [tuple(r) for r in A[inn.any(1) & ~inn.all(1)].tolist()]
            mp, mt, mq = mixed_rows(P3), mixed_rows(T3), mixed_rows(Q3)
        else:
            pairs, triples, quads, badc = _structure(faces, D)
            mixed = lambda s: any(isin(a) for a in s) and not all(isin(a) for a in s)
            mp = sorted(p for p in pairs if mixed(p))
            mt = sorted(t for t in triples if mixed(t))
            mq = sorted(q for q in quads if mixed(q))
        # every crystalline region carries its own frame (no reference bulk is needed): each connected piece of the
        # boundary faces a single region, since the outer atoms of two adjacent faces are equal or share a face
        good = good0
        why = dict(badc=bool(badc), wall_pair=sum(min(p) < 0 for p in mp), bad_cell=sum(not good[a] for p in mp for a in p if not isin(a)), triple=0, quad=0)
        open_ = bool(badc) or why["wall_pair"] > 0 or why["bad_cell"] > 0
        pid = {p: k for k, p in enumerate(mp)}
        links, tid = [], {}                       # dual links: mixed triple -> its two mixed pairs
        for t in mt:
            ps = [pid[q] for q in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2])) if q in pid]
            if len(ps) != 2:
                why["triple"] += 1; open_ = True; continue
            tid[t] = len(links); links.append((ps[0], ps[1]))
        adj = defaultdict(list)
        for k, (a, b) in enumerate(links):
            adj[a].append((b, k)); adj[b].append((a, k))
        qlinks = []                               # each quadruple as a polygon of dual links
        for q in mq:
            ks = [tid[t] for t in ((q[0], q[1], q[2]), (q[0], q[1], q[3]), (q[0], q[2], q[3]), (q[1], q[2], q[3]))
                  if t in tid]
            qlinks.append(ks)
        lq = defaultdict(list)                    # quadruples on each link (two, at the ends of the Voronoi edge)
        for qi, ks in enumerate(qlinks):
            for k in ks:
                lq[k].append(qi)
        why["quad"] = sum(len(v) != 2 for v in lq.values())
        if why["quad"]:
            open_ = True
        outer = lambda p: mp[p][1] if isin(mp[p][0]) else mp[p][0]      # the first-layer atom of a face of Sigma
        seen = np.zeros(len(mp), bool)
        res = dict(cells=D, open=open_, open_why=why, pieces=[])
        for s in range(len(mp)):
            if seen[s]:
                continue
            # spanning tree of the dual graph on this piece
            par = {s: (None, None)}; seen[s] = True; dq = deque([s]); nodes = []
            while dq:
                a = dq.popleft(); nodes.append(a)
                for b, k in adj[a]:
                    if not seen[b]:
                        seen[b] = True; par[b] = (a, k); dq.append(b)
            ks = sorted(set(k for a in nodes for _, k in adj[a]))
            qs = sorted(set(qi for k in ks for qi in lq[k]))
            chi = len(qs) - len(ks) + len(nodes)
            # loop basis on the valid part of the piece: faces whose outer atom is kept (atoms outside the analysed
            # region, e.g. beyond a cylinder's radius, carry no gluings). Tree-cotree on that surface with holes:
            # a spanning forest of its faces, then of its vertices through the remaining links, with all vertices
            # touching an excluded face merged into one (so loops around the holes are kept, except their sum)
            okn = {a for a in nodes if good[outer(a)]}
            if chi == 2 and len(okn) == len(nodes):
                # a sphere with every face valid: every loop on it contracts, so it has no curves and no content
                res["pieces"].append(dict(n_faces=len(nodes), chi=chi, genus=0.0, n_curves=0, n_invalid=0, curves=[],
                                          region=None, regions=0, n_other_region=0, burgers=[], vectors=[]))
                continue
            vks = [k for k in ks if links[k][0] in okn and links[k][1] in okn]
            par = {}
            for s0 in okn:
                if s0 in par:
                    continue
                par[s0] = (None, None); dq = deque([s0])
                while dq:
                    a = dq.popleft()
                    for b, k in adj[a]:
                        if b in okn and b not in par:
                            par[b] = (a, k); dq.append(b)
            tree = set(k for (_, k) in par.values() if k is not None)
            vkset = set(vks)
            fullq = {q for q in qs if all(k in vkset for k in qlinks[q])}
            qnode = lambda q: q if q in fullq else -1                   # -1: the merged excluded part
            qadj = defaultdict(list)
            for k in vks:
                if k in tree:
                    continue
                ends = [qnode(q) for q in lq[k]]
                if len(ends) == 2:
                    qadj[ends[0]].append((ends[1], k)); qadj[ends[1]].append((ends[0], k))
            cot, qseen = set(), set()
            for q0 in list(qadj):
                if q0 in qseen:
                    continue
                qseen.add(q0); dq = deque([q0])
                while dq:
                    q = dq.popleft()
                    for q2, k in qadj[q]:
                        if q2 not in qseen:
                            qseen.add(q2); cot.add(k); dq.append(q2)
            left = [k for k in vks if k not in tree and k not in cot]

            def to_root(a):
                path = [a]
                while par[path[-1]][0] is not None:
                    path.append(par[path[-1]][0])
                return path
            curves = []
            for k in left:
                a, b = links[k]
                pa, pb = to_root(a), to_root(b)
                while len(pa) > 1 and len(pb) > 1 and pa[-2] == pb[-2]:
                    pa.pop(); pb.pop()
                cyc = pa + pb[::-1][1:]                       # a ... lca ... b, closed by the link b-a
                walk = []
                for p in cyc:
                    x = outer(p)
                    if not walk or walk[-1] != x:
                        walk.append(x)
                while len(walk) > 1 and walk[-1] == walk[0]:
                    walk.pop()
                stp = conn.minimage(np.diff(conn.pts[walk + walk[:1]], axis=0)).sum(0)
                if L.shape == (3, 3):                         # triclinic: wraps are the fractional coordinates
                    wrap = np.round(stp @ np.linalg.inv(L)).astype(int)
                else:
                    wrap = np.where(per, np.round(stp / L), 0).astype(int)
                if len(walk) <= 2:                            # A -> B -> A: trivial holonomy
                    b_, rid_ = np.zeros(3), True
                elif not np.all(good[walk]):                  # through atoms outside the analysed region
                    b_, rid_ = None, None
                else:
                    bw = _bonded_walk(conn, walk, inD)
                    bb = None if bw is None else conn.burgers(bw)
                    b_, rid_ = (None, None) if bb is None else (np.round(bb["lattice"], 6), bb["R_identity"])
                curves.append(dict(walk=walk, wrap=wrap, b=b_, R_identity=rid_, region=int(regs[walk[0]])))
            valid = [c for c in curves if c["b"] is not None and c["R_identity"]]
            piece = dict(n_faces=len(nodes), chi=chi, genus=(2 - chi) / 2, n_curves=len(curves),
                         n_invalid=len(curves) - len(valid), curves=curves)
            # holonomies are expressed in the frame of the crystalline region the piece faces (one region per piece;
            # regions_facing checks this), so the content is read in that region's frame
            piece["regions_facing"] = len(set(int(regs[outer(a)]) for a in okn) - {-1})
            byreg = defaultdict(list)
            for c in valid:
                byreg[c["region"]].append(c)
            # a piece that faces several regions (possible only for close-packed crystals; none did in our samples)
            # is reported, not combined in one of them
            main = next(iter(byreg)) if len(byreg) == 1 else None
            piece["region"] = main                  # the frame of the content: that crystalline region's setting
            piece["regions"] = len(byreg)
            piece["n_other_region"] = sum(len(v) for r, v in byreg.items() if r != main)
            use = byreg.get(main, [])
            if use:
                W = np.array([c["wrap"] for c in use])
                # each curve's net translation is a translation of the crystal, exact up to the rounding of the input
                # structure, which accumulates along long walks (about 1e-2 over hundreds of steps): snapped to it
                B = _snap(conn, np.array([c["b"] for c in use]))
                piece["burgers"] = _nonwrapping_span(W, B)
                piece["vectors"] = []                  # set below: the classes of the shortest circuits through
                piece["shell_atoms"] = np.array(sorted({outer(a) for a in nodes}), np.int64)   # its shell atoms
            else:
                piece["burgers"], piece["vectors"] = [], []
            res["pieces"].append(piece)
        out.append(res)
    # the vectors of each piece: the classes of the shortest circuits through its shell atoms (shell_classes), which
    # depend on no choice of basis; the curves above give the lattice they span
    carry_i = [k for k, res in enumerate(out) if any(len(pc.get("burgers", [])) for pc in res["pieces"])]
    carry = [out[k] for k in carry_i]
    if G is not None:
        scs = shell_classes_many(conn, Ds, graph=G, which=carry_i)
        scs = [scs[k] for k in carry_i]
    else:
        scs = shell_classes_many(conn, [res["cells"] for res in carry])
    for res, sc in zip(carry, scs):
        res["shell"] = sc
        if sc is None:
            continue
        for pc in res["pieces"]:
            if not len(pc.get("burgers", [])):
                continue
            at = pc.pop("shell_atoms")
            k = np.minimum(np.searchsorted(sc["layer"], at), len(sc["layer"]) - 1)
            mk = np.bitwise_or.reduce(sc["mask"][k[sc["layer"][k] == at]]) if len(at) else 0
            pc["vectors"] = [sc["classes"][b] for b in range(len(sc["classes"])) if int(mk) >> b & 1]
    return out


def _spheres_batch(conn, Ds, allP, allT, allQ, sP, eP, sT, eT, sQ, eQ, allbad, good, regs):
    """The clusters whose boundary is a set of spheres facing crystal, found for all clusters at once: their results
    (as the per-cluster computation gives them) by cluster index. Euler characteristics add over pieces, so each piece's
    is counted from the mixed pairs, triples, and quadruples on it. A cluster with any other piece or an open boundary
    is left to the per-cluster computation."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    n = len(conn.pts); nc = len(Ds)
    label = np.full(n, -1, np.int64)
    for ci, D in enumerate(Ds):
        label[D] = ci
    def clus_of(starts, ends, N):
        c = np.full(N, -1, np.int64)
        for ci in range(nc):
            c[starts[ci]:ends[ci]] = ci
        return c
    cP, cT, cQ = clus_of(sP, eP, len(allP)), clus_of(sT, eT, len(allT)), clus_of(sQ, eQ, len(allQ))
    def mixed(A, cA):
        inn = (A >= 0) & (label[np.maximum(A, 0)] == cA[:, None])
        return inn.any(1) & ~inn.all(1)
    mP, mT, mQ = mixed(allP, cP), mixed(allT, cT), mixed(allQ, cQ)
    P, cp = allP[mP], cP[mP]; T, ct = allT[mT], cT[mT]; Q, cq = allQ[mQ], cQ[mQ]
    ids = np.unique(np.r_[P.ravel(), T.ravel(), Q.ravel()])
    B = np.int64(len(ids) + 1)
    cpt = lambda A: np.searchsorted(ids, A).astype(np.int64)
    Pc, Tc, Qc = cpt(P), cpt(T), cpt(Q)
    pkey = Pc[:, 0] * B + Pc[:, 1]; po = np.argsort(pkey); pks = pkey[po]
    def findp(a, b):
        k = a * B + b; i = np.minimum(np.searchsorted(pks, k), max(len(pks) - 1, 0))
        return (po[i] if len(pks) else i), (pks[i] == k if len(pks) else np.zeros(len(k), bool))
    # links (triples with exactly two mixed pairs) and their quadruple counts
    idx = np.zeros((len(T), 3), np.int64); pres = np.zeros((len(T), 3), bool)
    for j, (x, y) in enumerate([(0, 1), (0, 2), (1, 2)]):
        idx[:, j], pres[:, j] = findp(Tc[:, x], Tc[:, y])
    cnt = pres.sum(1)
    tri_bad = np.bincount(ct[cnt != 2], minlength=nc)
    valid = cnt == 2
    j1 = np.argmax(pres, 1); j2 = 2 - np.argmax(pres[:, ::-1], 1)
    la, lb = idx[np.arange(len(T)), j1][valid], idx[np.arange(len(T)), j2][valid]
    # triple keys: (index of its first two atoms as a pair-prefix) * B + third
    pre = Tc[valid, 0] * B + Tc[valid, 1]
    upre = np.unique(pre)
    tkey = np.searchsorted(upre, pre).astype(np.int64) * B + Tc[valid, 2]
    to = np.argsort(tkey); tks = tkey[to]
    nl = int(valid.sum())
    qcount = np.zeros(nl, np.int64)
    qfirst = np.full(len(Q), -1, np.int64)                       # a link of each quadruple (its piece)
    for (x, y, z) in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]:
        pk = Qc[:, x] * B + Qc[:, y]
        pi = np.minimum(np.searchsorted(upre, pk), max(len(upre) - 1, 0))
        okp = (upre[pi] == pk) if len(upre) else np.zeros(len(Q), bool)
        k = pi.astype(np.int64) * B + Qc[:, z]
        ti = np.minimum(np.searchsorted(tks, k), max(len(tks) - 1, 0))
        f = okp & ((tks[ti] == k) if len(tks) else np.zeros(len(Q), bool))
        li = to[ti[f]] if len(tks) else ti[f]
        np.add.at(qcount, li, 1)
        qfirst[f & (qfirst < 0)] = li[(qfirst[f] < 0)]
    lc = cp[la] if nl else np.zeros(0, np.int64)
    quad_bad = np.bincount(lc[(qcount != 2) & (qcount > 0)], minlength=nc) if nl else np.zeros(nc, np.int64)
    # pieces: components of the dual graph (nodes: mixed pairs; edges: links)
    npn = len(P)
    dual = csr_matrix((np.ones(2 * nl), (np.r_[la, lb], np.r_[lb, la])), shape=(npn, npn))
    ncp, piece = connected_components(dual, directed=False)
    nodes_pp = np.bincount(piece, minlength=ncp)
    links_pp = np.bincount(piece[la], minlength=ncp) if nl else np.zeros(ncp, np.int64)
    hasq = qfirst >= 0
    quads_pp = np.bincount(piece[la[qfirst[hasq]]], minlength=ncp) if hasq.any() else np.zeros(ncp, np.int64)
    chi = quads_pp - links_pp + nodes_pp
    outer = np.where(label[np.maximum(P[:, 0], 0)] == cp, P[:, 1], P[:, 0])
    gi = np.where(outer < 0, n - 1, outer)
    okn = good[gi]
    bad_face = np.bincount(piece[~okn], minlength=ncp) if (~okn).any() else np.zeros(ncp, np.int64)
    sphere = (chi == 2) & (bad_face == 0)
    pc_clus = np.full(ncp, -1, np.int64); pc_clus[piece] = cp
    nonsphere = np.bincount(pc_clus[~sphere], minlength=nc) if (~sphere).any() else np.zeros(nc, np.int64)
    # open reasons (as in the per-cluster computation)
    wall = np.bincount(cp[P.min(1) < 0], minlength=nc) if len(P) else np.zeros(nc, np.int64)
    outside = (P >= 0) & (label[np.maximum(P, 0)] != cp[:, None]) | (P < 0)
    badcell = np.bincount(np.repeat(cp, 2)[(outside & ~good[np.where(P < 0, n - 1, P)]).ravel()], minlength=nc)
    badset = np.zeros(nc, bool)
    for x in allbad:
        if label[x] >= 0:
            badset[label[x]] = True
    out = {}
    order = np.argsort(pc_clus, kind="stable")
    firstnode = np.full(ncp, npn); np.minimum.at(firstnode, piece, np.arange(npn))
    for ci in range(nc):
        if nonsphere[ci] or wall[ci] or badcell[ci] or tri_bad[ci] or quad_bad[ci] or badset[ci]:
            continue
        out[ci] = None
    if not out:
        return {}
    by = defaultdict(list)
    for pc_ in np.argsort(firstnode, kind="stable"):
        c_ = pc_clus[pc_]
        if c_ in out:
            by[c_].append(pc_)
    for ci in list(out):
        why = dict(badc=False, wall_pair=0, bad_cell=0, triple=0, quad=0)
        pieces = [dict(n_faces=int(nodes_pp[pc_]), chi=int(chi[pc_]), genus=0.0, n_curves=0, n_invalid=0, curves=[],
                       region=None, regions=0, n_other_region=0, burgers=[], vectors=[]) for pc_ in by[ci]]
        out[ci] = dict(cells=Ds[ci], open=False, open_why=why, pieces=pieces)
    return out


def _pieces_fast(conn, D, inD, P, T, Q, badc, good, regs, G, c):
    """The pieces of one cluster's boundary, as the walk-based loop of analyse computes them, from the integer voltages
    of its shell (shell_graph_many): a curve's net translation and winding are the voltage differences of its two ends
    along the tree of the dual graph plus the voltage of its closing step, so no walk is followed. A curve that takes a
    step without a gluing or across a flagged edge falls back to following its walk (_bonded_walk, conn.burgers)."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components, breadth_first_order
    n = len(conn.pts)
    gi = lambda a: np.where(a < 0, n - 1, a)               # (good[-1] in the walk-based loop: the last atom)
    isin = lambda a: (a >= 0) & inD[np.maximum(a, 0)]
    why = dict(badc=bool(badc), wall_pair=int((P.min(axis=1) < 0).sum()) if len(P) else 0, triple=0, quad=0)
    outside = ~isin(P)
    why["bad_cell"] = int((outside & ~good[gi(P)]).sum())
    open_ = bool(badc) or why["wall_pair"] > 0 or why["bad_cell"] > 0
    ids = np.unique(np.r_[P.ravel(), T.ravel(), Q.ravel()])
    Bk = np.int64(len(ids) + 1)
    if Bk ** 3 > np.iinfo(np.int64).max // 4:
        return None                                        # (too many atoms for the keys: the caller walks)
    cp = lambda A: np.searchsorted(ids, A).astype(np.int64)
    Pc, Tc, Qc = cp(P), cp(T), cp(Q)
    pkey = Pc[:, 0] * Bk + Pc[:, 1] if len(P) else np.zeros(0, np.int64)
    def find(sorted_keys, k):
        if len(sorted_keys) == 0:
            return np.zeros(len(k), np.int64), np.zeros(len(k), bool)
        p_ = np.minimum(np.searchsorted(sorted_keys, k), len(sorted_keys) - 1)
        return p_, sorted_keys[p_] == k
    # dual links: each mixed triple joins its two mixed pairs
    sub = [(0, 1), (0, 2), (1, 2)]
    pp = np.zeros((len(T), 3), np.int64); pr = np.zeros((len(T), 3), bool)
    for j, (x, y) in enumerate(sub):
        pp[:, j], pr[:, j] = find(pkey, Tc[:, x] * Bk + Tc[:, y]) if len(T) else (pp[:, j], pr[:, j])
    cnt = pr.sum(axis=1)
    why["triple"] = int((cnt != 2).sum())
    open_ = open_ or why["triple"] > 0
    vt = np.flatnonzero(cnt == 2)
    prv = pr[vt]; ppv = pp[vt]
    j1 = np.argmax(prv, axis=1); j2 = 2 - np.argmax(prv[:, ::-1], axis=1)       # first and last present sub-pair
    la = ppv[np.arange(len(vt)), j1].astype(np.int64); lb = ppv[np.arange(len(vt)), j2].astype(np.int64)
    nl = len(vt)
    tkey = (Tc[vt, 0] * Bk + Tc[vt, 1]) * Bk + Tc[vt, 2]           # sorted (rows sorted)
    # quadruples: their links (valid sub-triples), and the quadruples on each link
    qsub = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    ql_link, ql_quad = [], []
    for (x, y, z) in qsub:
        if len(Q):
            k_, f_ = find(tkey, (Qc[:, x] * Bk + Qc[:, y]) * Bk + Qc[:, z])
            ql_link.append(np.where(f_, k_, -1)); ql_quad.append(np.arange(len(Q)))
    QL = np.stack(ql_link, 1) if ql_link else np.zeros((0, 4), np.int64)       # (nq, 4) link ids, -1: absent
    pairs_lq = np.concatenate([np.c_[QL[QL[:, j] >= 0, j], np.flatnonzero(QL[:, j] >= 0)] for j in range(QL.shape[1])]) \
        if len(QL) else np.zeros((0, 2), np.int64)
    nq_on = np.bincount(pairs_lq[:, 0], minlength=nl) if len(pairs_lq) else np.zeros(nl, np.int64)
    if len(pairs_lq):                                   # quadruples of each link, by link (in order of quadruple)
        o_ = np.lexsort((pairs_lq[:, 1], pairs_lq[:, 0])); pairs_lq = pairs_lq[o_]
    lq_start = np.searchsorted(pairs_lq[:, 0], np.arange(nl)) if len(pairs_lq) else np.zeros(nl, np.int64)
    lq_end = np.searchsorted(pairs_lq[:, 0], np.arange(nl), "right") if len(pairs_lq) else np.zeros(nl, np.int64)
    why["quad"] = int(((nq_on != 2) & (nq_on > 0)).sum())
    open_ = open_ or why["quad"] > 0
    outer = np.where(isin(P[:, 0]), P[:, 1], P[:, 0]) if len(P) else np.zeros(0, np.int64)
    res = dict(cells=D, open=open_, open_why=why, pieces=[])
    npn = len(P)
    if npn == 0:
        return res
    dual = csr_matrix((np.ones(2 * nl), (np.r_[la, lb], np.r_[lb, la])), shape=(npn, npn))
    ncp, lab = connected_components(dual, directed=False)
    first = np.full(ncp, npn); np.minimum.at(first, lab, np.arange(npn))
    rank = np.argsort(np.argsort(first)); lab = rank[lab]               # pieces in order of their first face
    link_piece = lab[la] if nl else np.zeros(0, np.int64)
    okn_all = good[gi(outer)]
    # step voltages outer(a) -> outer(b) of the shell (0 for the same atom)
    def steps(xa, xb):
        same = xa == xb
        e, f = _edge(G, c, xa, xb)
        f = f & G["ok"][e] if len(e) else f
        lat = np.where((same | ~f)[:, None], 0.0, G["lat"][e]) if len(e) else np.zeros((0, 3))
        img = np.where((same | ~f)[:, None], 0, G["img"][e]) if len(e) else np.zeros((0, 3), np.int64)
        return lat, img, same | f
    for pc in range(ncp):
        nodes = np.flatnonzero(lab == pc)
        ks = np.flatnonzero(link_piece == pc)
        qs = np.unique(np.concatenate([pairs_lq[lq_start[k]:lq_end[k], 1] for k in ks])) if len(pairs_lq) and len(ks) else np.zeros(0, np.int64)
        chi = len(qs) - len(ks) + len(nodes)
        okn = okn_all[nodes]
        if chi == 2 and okn.all():
            res["pieces"].append(dict(n_faces=len(nodes), chi=chi, genus=0.0, n_curves=0, n_invalid=0, curves=[],
                                      region=None, regions=0, n_other_region=0, burgers=[], vectors=[]))
            continue
        okmask = np.zeros(npn, bool); okmask[nodes[okn]] = True
        vks = ks[okmask[la[ks]] & okmask[lb[ks]]]
        # spanning forest of the valid faces (dual tree), from a virtual root joined to each component's first face
        g = csr_matrix((np.ones(2 * len(vks)), (np.r_[la[vks], lb[vks]], np.r_[lb[vks], la[vks]])), shape=(npn, npn))
        okn_nodes = nodes[okn]
        nc2, lab2 = connected_components(g, directed=False)
        _, fi_ = np.unique(lab2[okn_nodes], return_index=True)             # each component's first face, in order
        roots = okn_nodes[np.sort(fi_)].astype(np.int64)
        g2 = csr_matrix((np.ones(2 * len(vks) + len(roots)), (np.r_[la[vks], lb[vks], np.full(len(roots), npn)],
                                                              np.r_[lb[vks], la[vks], roots])), shape=(npn + 1, npn + 1))
        order, pred = breadth_first_order(g2, npn, directed=True, return_predecessors=True)
        order = order[1:]
        # the tree link of each non-root face: a link joining it to its parent (the smallest such)
        lkey = np.r_[la[vks] * npn + lb[vks], lb[vks] * npn + la[vks]]
        lid = np.r_[vks, vks]
        o = np.lexsort((lid, lkey)); lkey, lid = lkey[o], lid[o]
        nonroot = order[pred[order] != npn]
        tk = np.searchsorted(lkey, pred[nonroot].astype(np.int64) * npn + nonroot)
        intree = np.zeros(nl, bool); intree[lid[tk]] = True
        # cotree: a spanning forest of the vertices through the remaining valid links (vertices touching an excluded
        # face merged into one)
        vkset = np.zeros(nl, bool); vkset[vks] = True
        if len(qs):
            Lq = QL[qs]
            fullq = qs[np.all((Lq < 0) | vkset[np.maximum(Lq, 0)], axis=1)]
        else:
            fullq = np.zeros(0, np.int64)
        merged = len(fullq)
        qmap = np.full(len(Q) + 1, merged, np.int64); qmap[fullq] = np.arange(len(fullq))
        cand = vks[~intree[vks]]
        cand = cand[(lq_end[cand] - lq_start[cand]) == 2]
        ce = np.c_[qmap[pairs_lq[lq_start[cand], 1]], qmap[pairs_lq[lq_start[cand] + 1, 1]]] if len(cand) else np.zeros((0, 2), np.int64)
        cl = cand
        incot = np.zeros(nl, bool)
        if len(ce):
            nqn = merged + 1
            qg = csr_matrix((np.ones(2 * len(ce)), (np.r_[ce[:, 0], ce[:, 1]], np.r_[ce[:, 1], ce[:, 0]])), shape=(nqn, nqn))
            nc3, lab3 = connected_components(qg, directed=False)
            ck = np.r_[ce[:, 0] * nqn + ce[:, 1], ce[:, 1] * nqn + ce[:, 0]]; cid_ = np.r_[cl, cl]
            o = np.lexsort((cid_, ck)); ck, cid_ = ck[o], cid_[o]
            touched = np.unique(ce.ravel())
            # one search from a virtual vertex joined to a root of every component (the merged vertex where it is in it)
            lt = lab3[touched]
            _, fi3 = np.unique(lt, return_index=True); fi3 = np.sort(fi3)       # components in order of appearance
            r3 = touched[fi3].astype(np.int64)                                   # their first vertex,
            if len(touched) and touched[-1] == merged:                           # or the merged one
                r3[lt[fi3] == lab3[merged]] = merged
            qg2 = csr_matrix((np.ones(2 * len(ce) + len(r3)), (np.r_[ce[:, 0], ce[:, 1], np.full(len(r3), nqn)],
                                                               np.r_[ce[:, 1], ce[:, 0], r3])), shape=(nqn + 1, nqn + 1))
            ordq, prq = breadth_first_order(qg2, nqn, directed=True, return_predecessors=True)
            vq = ordq[1:]; vq = vq[prq[vq] != nqn]
            incot[cid_[np.searchsorted(ck, prq[vq].astype(np.int64) * nqn + vq)]] = True
        left = vks[~intree[vks] & ~incot[vks]]
        # potentials along the dual tree: the summed step voltages from each root, and whether a step was unglued
        Phi = np.zeros((npn, 3)); Pim = np.zeros((npn, 3), np.int64); taint = np.zeros(npn, bool)
        # depths in the tree: breadth-first distances from the virtual root (order is by depth, so each level is a
        # contiguous run of it)
        from scipy.sparse.csgraph import shortest_path
        dd = shortest_path(g2, directed=True, unweighted=True, indices=npn)[order].astype(np.int64) if len(order) else np.zeros(0, np.int64)
        lev_start = np.searchsorted(dd, np.arange(int(dd.max()) + 2 if len(dd) else 0))
        for lev in range(2, int(dd.max()) + 1 if len(dd) else 0):
            vs = order[lev_start[lev]:lev_start[lev + 1]]; ps = pred[vs]
            sl, si, sok = steps(outer[ps], outer[vs])
            Phi[vs] = Phi[ps] + sl; Pim[vs] = Pim[ps] + si; taint[vs] = taint[ps] | ~sok
        curves = []
        if len(left):
            L_ = np.asarray(left, np.int64); A_, B_ = la[L_], lb[L_]
            sl, si, sok = steps(outer[B_], outer[A_])
            bval = Phi[B_] - Phi[A_] + sl; wval = Pim[B_] - Pim[A_] + si
            bad = taint[A_] | taint[B_] | ~sok
            for j, k in enumerate(left):
                a = int(A_[j])
                if not bad[j]:
                    curves.append(dict(wrap=wval[j].astype(int), b=np.round(bval[j], 6), R_identity=True,
                                       region=int(regs[outer[a]])))
                    continue
                # fall back: follow the walk of this curve (a ... lca ... b, closed by the link)
                b_ = int(B_[j])
                def to_root(x):
                    path = [x]
                    while pred[path[-1]] != npn and pred[path[-1]] >= 0:
                        path.append(int(pred[path[-1]]))
                    return path
                pa, pb = to_root(a), to_root(b_)
                while len(pa) > 1 and len(pb) > 1 and pa[-2] == pb[-2]:
                    pa.pop(); pb.pop()
                cyc = pa + pb[::-1][1:]
                walk = []
                for x in cyc:
                    xo = int(outer[x])
                    if not walk or walk[-1] != xo:
                        walk.append(xo)
                while len(walk) > 1 and walk[-1] == walk[0]:
                    walk.pop()
                if len(walk) <= 2:
                    b2, rid2 = np.zeros(3), True
                elif not np.all(good[walk]):
                    b2, rid2 = None, None
                else:
                    bw = _bonded_walk(conn, walk, inD)
                    bb = None if bw is None else conn.burgers(bw)
                    b2, rid2 = (None, None) if bb is None else (np.round(bb["lattice"], 6), bb["R_identity"])
                curves.append(dict(walk=walk, wrap=wval[j].astype(int), b=b2, R_identity=rid2, region=int(regs[walk[0]])))
        valid = [cv for cv in curves if cv["b"] is not None and cv["R_identity"]]
        piece = dict(n_faces=len(nodes), chi=chi, genus=(2 - chi) / 2, n_curves=len(curves),
                     n_invalid=len(curves) - len(valid), curves=curves)
        piece["regions_facing"] = len(set(int(x) for x in regs[outer[nodes[okn]]]) - {-1})
        byreg = defaultdict(list)
        for cv in valid:
            byreg[cv["region"]].append(cv)
        main = next(iter(byreg)) if len(byreg) == 1 else None    # several regions: reported, not combined
        piece["region"] = main
        piece["regions"] = len(byreg)
        piece["n_other_region"] = sum(len(v) for r, v in byreg.items() if r != main)
        use = byreg.get(main, [])
        if use:
            W = np.array([cv["wrap"] for cv in use]); Bm = _snap(conn, np.array([cv["b"] for cv in use]))
            piece["burgers"] = _nonwrapping_span(W, Bm)
            piece["vectors"] = []
            piece["shell_atoms"] = np.unique(outer[nodes]).astype(np.int64)
        else:
            piece["burgers"], piece["vectors"] = [], []
        res["pieces"].append(piece)
    return res


def _nbrs(conn, a):
    """Face neighbors of atom a from the connection (all atoms), without walls."""
    c = getattr(conn, "vcells", conn.cells)[a]
    return [] if c is None else [w for w in c.nbr if w >= 0]


_MERIDIANS = os.environ.get("VOROTOP_BURGERS_MERIDIANS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin",
                                                                      "vb_meridians"))   # built with the package


def _rank(rows):
    """Exact rank of a small integer matrix (list of rows)."""
    M = [[Fraction(int(v)) for v in r] for r in rows]
    rank, ncol = 0, len(M[0]) if M else 0
    for col in range(ncol):
        p = next((i for i in range(rank, len(M)) if M[i][col] != 0), None)
        if p is None:
            continue
        M[rank], M[p] = M[p], M[rank]
        for i in range(len(M)):
            if i != rank and M[i][col] != 0:
                f = M[i][col] / M[rank][col]
                M[i] = [a - f * b for a, b in zip(M[i], M[rank])]
        rank += 1
    return rank


def _snap(conn, X):
    """Each row of X (lattice coordinates) moved to the nearest translation of the crystal: integer coordinates plus
    one of its centering vectors. The net translation of a closed loop whose net rotation is the identity is such a
    translation, exactly; the input structure's rounding accumulates along long walks and is removed here."""
    X = np.atleast_2d(np.asarray(X, float))
    C = np.asarray(getattr(conn.ref, "centerings", np.zeros((1, 3))), float)
    A = np.asarray(getattr(conn.ref, "A", np.eye(3)), float)
    best, bd = None, None
    for c in C:
        cand = np.rint(X - c) + c
        d = np.linalg.norm((cand - X) @ A, axis=1)
        if best is None:
            best, bd = cand, d
        else:
            m = d < bd - 1e-12
            best[m], bd[m] = cand[m], d[m]
    return best


def _csr(conn):
    """Face-sharing neighbors of every atom as (offsets, data) (negative entries: walls)."""
    vc = getattr(conn, "vcells", None)
    if vc is not None and hasattr(vc, "_o"):
        return vc._o, vc._d
    if hasattr(conn, "_nb_offs"):
        return conn._nb_offs, conn._nb_data
    if not hasattr(conn, "_csr_cache"):
        lst = [[] if c is None else list(c.nbr) for c in conn.cells]
        offs = np.r_[0, np.cumsum([len(l) for l in lst])].astype(np.int64)
        conn._csr_cache = (offs, np.array([x for l in lst for x in l], np.int64))
    return conn._csr_cache


def _gather(offs, data, atoms):
    """(position in atoms, neighbor) for every neighbor of the given atoms (vectorized)."""
    atoms = np.asarray(atoms, np.int64)
    lens = offs[atoms + 1] - offs[atoms]
    start = np.repeat(offs[atoms], lens)
    k = np.arange(int(lens.sum())) - np.repeat(np.cumsum(lens) - lens, lens)
    return np.repeat(np.arange(len(atoms)), lens), data[start + k]


def shell_classes(conn, D, extra=8, spacing=3, budget=0):
    """shell_classes_many for one cluster."""
    return shell_classes_many(conn, [D], extra, spacing, budget)[0]


def shell_graph_many(conn, Ds):
    """The shells of the clusters Ds as one integer voltage graph, each shell with its own copy of its atoms (an atom
    can be in the shells of two clusters; walks never pass from one shell to another).

    Frames are fixed once along a spanning forest of the shells, each tree rooted where vectors are read (a cell with
    a crystal setting, the bulk region first). Each face-sharing pair u -> v of shell atoms with a gluing then carries a
    voltage: the net translation t_u + R_u d - t_v in lattice coordinates, which is zero on tree edges and otherwise the
    translation of the crystal that its cycle closes (snapped to it), and the box image it crosses. The net translation
    and winding of any closed walk in a shell are the sums of its voltages. An edge is flagged where its transported
    frame disagrees with the tree's (a rotation along some loop) or its tree has no crystal setting.

    Returns dict(atom, cid (cluster of each vertex), base (first vertex of each cluster), layers, U, V (directed
    edges, both directions), lat (n_edges, 3), img (n_edges, 3, int), ok (bool), comp (component of each vertex),
    region (of each vertex's tree root, -1 if none))."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components, breadth_first_order
    n = len(conn.pts)
    good = _good(conn)
    offs, data = _csr(conn)
    inD = np.zeros(n, bool)
    layers = []
    for D in Ds:
        D = np.asarray(D, np.int64); inD[D] = True
        _, nb = _gather(offs, data, D)
        nb = nb[nb >= 0]
        layers.append(np.unique(nb[~inD[nb] & good[nb]]))
        inD[D] = False
    sizes = np.array([len(l) for l in layers], np.int64); base = np.r_[0, np.cumsum(sizes)]
    m = int(base[-1])
    atom = np.concatenate(layers) if m else np.zeros(0, np.int64)
    cid = np.repeat(np.arange(len(Ds)), sizes)
    G = dict(atom=atom, cid=cid, base=base, layers=layers, n=n, m=m)
    if m == 0:
        G.update(U=np.zeros(0, np.int64), V=np.zeros(0, np.int64), lat=np.zeros((0, 3)), img=np.zeros((0, 3), np.int64),
                 ok=np.zeros(0, bool), comp=np.zeros(0, np.int64), region=np.zeros(0, np.int64), ncomp=0)
        return G
    key = cid.astype(np.int64) * n + atom
    ui, vg = _gather(offs, data, atom)
    kv = cid[ui].astype(np.int64) * n + vg
    pos = np.minimum(np.searchsorted(key, kv), m - 1)
    keep = (vg >= 0) & (key[pos] == kv)
    U, V = ui[keep], pos[keep]
    if hasattr(conn, "g_many"):
        hit, gR, gt = conn.g_many(atom[U], atom[V])
    else:
        Gl = [conn.g(atom[u], atom[v]) for u, v in zip(U, V)]
        hit = np.array([g is not None for g in Gl], bool)
        gR = np.array([g[0] if g is not None else np.eye(3) for g in Gl]).reshape(-1, 3, 3)
        gt = np.array([g[1] if g is not None else np.zeros(3) for g in Gl]).reshape(-1, 3)
    U, V, gR, gt = U[hit], V[hit], gR[hit], gt[hit]
    disp = conn.minimage(conn.pts[atom[V]] - conn.pts[atom[U]])
    main = conn.main_region()
    if hasattr(conn, "_Rg") and isinstance(getattr(conn, "region", None), np.ndarray):
        inR = np.asarray(conn.kept, bool)[atom] & (np.asarray(conn._Rg)[atom] >= 0)       # (as conn.Rglob)
        rof = np.where(inR, conn.region[atom], -1).astype(np.int64)
        readable = inR & np.isin(rof, list(conn.G))
    else:
        rof = np.array([conn.region_of.get(int(x), -1) for x in atom], np.int64)
        readable = np.array([int(x) in conn.Rglob for x in atom], bool) & np.isin(rof, list(conn.G))
    pref = np.lexsort((np.arange(m), rof != main, ~readable))
    graph = csr_matrix((np.ones(len(U)), (U, V)), shape=(m, m))
    ncomp, comp = connected_components(graph, directed=False)
    root_of = np.full(ncomp, -1)
    for x in pref:
        if root_of[comp[x]] < 0:
            root_of[comp[x]] = x
    g2 = csr_matrix((np.ones(len(U) + ncomp), (np.r_[U, np.full(ncomp, m)], np.r_[V, root_of])), shape=(m + 1, m + 1))
    order, pred = breadth_first_order(g2, m, directed=True, return_predecessors=True)
    order = order[1:]
    depth = np.zeros(m + 1, np.int64)
    for v in order:
        depth[v] = depth[pred[v]] + 1
    ekey = U.astype(np.int64) * m + V
    eord = np.argsort(ekey); esort = ekey[eord]
    Rf = np.tile(np.eye(3), (m, 1, 1)); tf = np.zeros((m, 3)); Pf = conn.pts[atom].astype(float).copy()
    d = depth[order]
    for lev in range(2, int(d.max()) + 1 if len(d) else 0):             # depth 1: the roots
        vs = order[d == lev]; ps = pred[vs]
        k = eord[np.searchsorted(esort, ps.astype(np.int64) * m + vs)]
        Rf[vs] = np.einsum("kab,kbc->kac", Rf[ps], gR[k]); tf[vs] = tf[ps] + np.einsum("kab,kb->ka", Rf[ps], gt[k])
        Pf[vs] = Pf[ps] + disp[k]
    ok = np.abs(np.einsum("kab,kbc->kac", Rf[U], gR) - Rf[V]).max(axis=(1, 2)) < 1e-3
    Mc = np.zeros((ncomp, 3, 3)); hasc = np.zeros(ncomp, bool); regc = np.full(ncomp, -1, np.int64)
    for c, r in enumerate(root_of):
        if readable[r]:                     # root-frame translation -> lattice coordinates, as in conn.burgers
            x = int(atom[r])
            Rx = conn.ref.R[conn._Rg[x]] if hasattr(conn, "_Rg") else conn.Rglob[x]     # (arrays: no dictionaries)
            rx = int(rof[r])
            Mc[c] = np.array([conn.ref.lattice_coords(conn.G[rx] @ -(Rx @ e)) for e in np.eye(3)]).T
            hasc[c] = True; regc[c] = rx
    Mv, hasM, region = Mc[comp], hasc[comp], regc[comp]
    ok &= hasM[U]
    t = tf[U] + np.einsum("kab,kb->ka", Rf[U], gt) - tf[V]
    lat = np.einsum("kab,kb->ka", Mv[U], t)
    nz = np.abs(lat).max(axis=1) > 1e-6     # tree edges carry zero; the others close a cycle: a crystal translation
    lat[nz] = _snap(conn, lat[nz]); lat[~nz] = 0
    Hm = np.asarray(conn.L, float) if np.ndim(conn.L) == 2 else np.diag(np.asarray(conn.L, float))
    img = np.round((Pf[U] + disp - Pf[V]) @ np.linalg.inv(Hm)).astype(np.int64)
    G.update(U=U, V=V, lat=lat, img=img, ok=ok, comp=comp, region=region, ncomp=ncomp, esort=esort, eord=eord)
    return G


def _edge(G, c, x, y):
    """Indices of the shell edges x -> y (atoms) in cluster c's shell, and whether each exists (vectorized)."""
    m, n = G["m"], G["n"]
    if "key" not in G:
        G["key"] = G["cid"].astype(np.int64) * n + G["atom"]
    key = G["key"]
    kx, ky = np.int64(c) * n + np.asarray(x, np.int64), np.int64(c) * n + np.asarray(y, np.int64)
    px = np.minimum(np.searchsorted(key, kx), max(m - 1, 0)); py = np.minimum(np.searchsorted(key, ky), max(m - 1, 0))
    have = (key[px] == kx) & (key[py] == ky)
    ek = px.astype(np.int64) * m + py
    p = np.minimum(np.searchsorted(G["esort"], ek), max(len(G["esort"]) - 1, 0))
    found = have & (G["esort"][p] == ek) if len(G["esort"]) else np.zeros(len(ek), bool)
    return G["eord"][p] if len(G["esort"]) else np.zeros(len(ek), np.int64), found


def shell_classes_many(conn, Ds, extra=8, spacing=3, budget=0, graph=None, which=None):
    """The circuit classes of a cluster's shell: for every shell atom x, the vectors (up to sign) of the shortest
    closed walks through x within the shell whose net translation is not zero and which do not wind around the box.

    With frames fixed once along a spanning forest of the shell, each face-sharing pair carries an integer voltage
    (the translation of the crystal it closes, and the box image) and the net translation of a closed walk is the sum
    of its voltages; the shortest walk of each class through x is then a shortest path in the covering graph, found
    exactly by vb_meridians (src/cpp) (which also certifies far from nodes: lengths change by at most 2 per step). A piece
    of the shell is searched only if some integer combination of its cycles has zero box image and nonzero translation
    (an exact rank test); every atom of it then has such a walk, so every search ends. budget > 0 caps the states of
    one search, leaving the atom unresolved (with a warning).

    All clusters are searched in one graph, each shell with its own copy of its atoms (an atom can be in the shells of
    two clusters, and walks never pass from one shell to another), so the result is that of separate searches.
    Returns per cluster dict(layer: shell atoms (sorted), mask: per shell atom a bitmask over classes, classes: the
    vectors (lattice coordinates, sign fixed), sorted), or None for a cluster without shell atoms."""
    import subprocess, time, warnings
    G = graph if graph is not None else shell_graph_many(conn, Ds)
    out = [None] * len(Ds)
    which = list(range(len(Ds))) if which is None else list(which)
    if G["m"] == 0 or not which:
        return out
    sel_c = np.zeros(len(Ds), bool); sel_c[which] = True
    vsel = sel_c[G["cid"]]
    vid = -np.ones(G["m"], np.int64); vid[vsel] = np.arange(int(vsel.sum()))   # vertices of the searched clusters
    m = int(vsel.sum())
    if m == 0:
        return out
    e = G["ok"] & vsel[G["U"]]
    U, V = vid[G["U"][e]], vid[G["V"][e]]
    volt = np.c_[np.rint(12 * G["lat"][e]), G["img"][e]].astype(np.int32)   # 12: a common denominator of centerings
    if not G["ok"][vsel[G["U"]]].all():
        warnings.warn(f"shell_classes: {int((~G['ok'][vsel[G['U']]]).sum())} shell edges dropped (frame rotation along a "
                      "loop, or a shell piece with no crystal setting); classes near them may not be exact")
    comp = G["comp"][vsel]
    ncomp = G["ncomp"]
    active = np.zeros(m, np.int8)            # pieces with a nonzero non-wrapping class (exact rank test)
    sel = (U < V) & np.any(volt != 0, axis=1)
    cU = comp[U[sel]]; vs_ = volt[sel]
    o = np.argsort(cU, kind="stable"); cU, vs_ = cU[o], vs_[o]
    actc = np.zeros(ncomp, np.int8)
    for c in np.unique(cU):
        a, b = np.searchsorted(cU, c), np.searchsorted(cU, c, "right")
        W = np.unique(vs_[a:b], axis=0)
        if _rank(W.T.tolist()) > _rank(W[:, 3:].T.tolist()):
            actc[c] = 1
    active = actc[comp]
    o = np.argsort(U, kind="stable"); U, V, volt = U[o], V[o], volt[o]
    offs2 = np.r_[0, np.cumsum(np.bincount(U, minlength=m))].astype(np.int32)
    buf = (np.array([m, len(V), 0, extra, spacing, budget], np.int32).tobytes() + offs2.tobytes()
           + V.astype(np.int32).tobytes() + volt.tobytes() + active.tobytes())
    t0 = time.time()
    raw = np.frombuffer(subprocess.run([_MERIDIANS], input=buf, capture_output=True, check=True).stdout, np.int32)
    if time.time() - t0 > 60:
        warnings.warn(f"shell_classes: {time.time() - t0:.0f} s for shells of {m} atoms. The search is slow where many "
                      "shell atoms are far from every branch and nearly as far from two of them, as on a large planar "
                      "fault attached to a dislocation network; the result is still exact")
    # least classes per shell atom (certified atoms carry their one class, with length -1)
    cnt = raw[:m]; rec = raw[m:len(raw) - 5].reshape(-1, 7)        # (the last 5: statistics of the search)
    if (cnt < 0).any():
        warnings.warn(f"shell_classes: {int((cnt < 0).sum())} shell atoms unresolved within the budget of {budget} states")
    owner = np.repeat(np.arange(m), np.maximum(cnt, 0))
    ln = rec[:, 6].astype(np.int64)
    big = np.iinfo(np.int64).max
    mn = np.full(m, big); np.minimum.at(mn, owner, np.where(ln >= 0, ln, big))
    keep = (ln < 0) | (ln == mn[owner])
    owner, rec = owner[keep], rec[keep, :3]
    gl_of = np.flatnonzero(vsel)                                    # searched vertex -> graph vertex
    for c in which:
        lo, hi = G["base"][c], G["base"][c + 1]
        if hi == lo:
            continue
        sel = (gl_of[owner] >= lo) & (gl_of[owner] < hi)
        keys, inv = (np.unique(rec[sel], axis=0, return_inverse=True) if sel.any() else (np.zeros((0, 3), int), np.zeros(0, int)))
        mask = np.zeros(hi - lo, np.int64)
        np.bitwise_or.at(mask, gl_of[owner[sel]] - lo, np.left_shift(np.int64(1), np.asarray(inv).ravel().astype(np.int64)))
        out[c] = dict(layer=G["layers"][c], mask=mask, classes=[np.array(kk, float) / 12.0 for kk in keys])
    return out


def branch_labels(conn, results, faces=None, extra=8, spacing=3, budget=0):
    """The branch of each atom of a cluster that carries Burgers vectors, exact and independent of any processing
    order. Shell atom x: the classes of the shortest circuits through x (shell_classes; computed by analyse and kept in
    the result). Core atom y: the classes whose shell atoms are nearest to y, counted in steps between face-sharing
    atoms through the cluster. Ties are kept.

    Returns (vec, classes): vec (n, 3) lattice coordinates, the first of the atom's classes in lexicographic order
    (NaN outside such clusters, or where no class was found); classes {atom: [vectors]} for every cluster atom (one
    vector, or several where they tie). A cluster with a single class gives it to every atom."""
    n = len(conn.pts)
    vec = np.full((n, 3), np.nan)
    allcls = {}
    offs, data = _csr(conn)
    for res in results:
        if not any(len(pc.get("vectors", [])) for pc in res["pieces"]):
            continue
        D = np.asarray(res["cells"])
        sc = res.get("shell") or shell_classes(conn, D, extra, spacing, budget)
        if sc is None or not sc["classes"]:
            continue
        cls = sc["classes"]
        if len(cls) == 1:
            vec[D] = cls[0]
            for a in D:
                allcls[int(a)] = [cls[0]]
            continue
        inD = np.zeros(n, bool); inD[D] = True
        # core atoms: nearest shell classes, through the cluster, ties kept (one breadth-first search carrying masks)
        cmask = {}
        frontier = {int(x): int(mk) for x, mk in zip(sc["layer"], sc["mask"]) if mk}
        seen = set(frontier)
        while frontier:
            nxt = {}
            for u, mk in frontier.items():
                for w in data[offs[u]:offs[u + 1]].tolist():
                    if w >= 0 and inD[w] and w not in seen:
                        nxt[w] = nxt.get(w, 0) | mk
            seen.update(nxt); cmask.update(nxt); frontier = nxt
        for a in D:
            mk = cmask.get(int(a), 0)
            lst = [cls[b] for b in range(len(cls)) if mk >> b & 1]
            allcls[int(a)] = lst
            if lst:
                vec[a] = lst[0]
    return vec, allcls
