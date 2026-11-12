"""Loops of perfect cells around a line, and automatic tracing of dislocation lines.

Geometry is used here only to decide WHERE to put a loop (around which axis, how far out); the
Burgers vector itself comes from the combinatorial holonomy (engine.FastConnection.burgers).
"""

from collections import deque

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import breadth_first_order


def _tree_coords(conn, x):
    """Coordinates for the spatial tree: wrapped into the box along periodic directions."""
    x = np.array(x, float)
    for k in range(3):
        if conn.periodic[k]:
            x[:, k] %= conn.L[k]
    return x


def _tree(conn):
    """cKDTree of the particle positions (periodic where the box is), cached on conn."""
    if getattr(conn, "_tree_cache", None) is None:
        from scipy.spatial import cKDTree
        X = _tree_coords(conn, conn.pts)
        lo = X.min(0)
        X = X - np.where(conn.periodic, 0.0, lo)          # non-periodic axes: shift to >= 0
        big = 4 * (X.max(0) + 1.0)
        box = np.where(conn.periodic, conn.L, big)          # non-periodic: a box large enough never to wrap
        conn._tree_shift = np.where(conn.periodic, 0.0, lo)
        conn._tree_cache = cKDTree(np.minimum(X, box * (1 - 1e-12)), boxsize=box)
    return conn._tree_cache


def _tri_tree(conn, r):
    """Triclinic box: a plain cKDTree of the points wrapped into the cell plus their periodic images within r of it
    (cached on conn, rebuilt for a larger r), and the particle index of every tree point."""
    c = getattr(conn, "_tri_cache", None)
    if c is None or c[0] < r:
        from scipy.spatial import cKDTree
        H = conn.L
        s = conn.pts @ np.linalg.inv(H)
        s -= np.floor(s)
        # fractional margin per axis: r over the distance between the two faces of the cell normal to that axis
        m = r * np.linalg.norm(np.linalg.inv(H), axis=0)
        S, I = [s], [np.arange(len(s))]
        for sh in np.stack(np.meshgrid(*[(-1, 0, 1)] * 3, indexing="ij"), -1).reshape(-1, 3):
            if not sh.any():
                continue
            t = s + sh
            ok = np.all((t > -m) & (t < 1 + m), axis=1)
            S.append(t[ok]); I.append(np.nonzero(ok)[0])
        S, I = np.concatenate(S), np.concatenate(I)
        c = conn._tri_cache = (r, cKDTree(S @ H), I, H)
    return c


def _near(conn, p, r):
    """Indices of the particles within r of the point p (minimum image), sorted."""
    if getattr(conn, "tri", False):
        _, tree, I, H = _tri_tree(conn, max(r, 30.0))
        s = np.asarray(p, float) @ np.linalg.inv(H)
        q = (s - np.floor(s)) @ H
        return np.unique(I[tree.query_ball_point(q, r)]).astype(np.int64)
    tree = _tree(conn)
    qp = _tree_coords(conn, np.asarray(p, float)[None])[0] - conn._tree_shift
    return np.array(sorted(tree.query_ball_point(np.mod(qp, tree.boxsize[:3]), r)), dtype=np.int64)


def _csr(conn):
    """Face adjacency of the cells as CSR arrays (offsets, neighbour ids; -1 = wall), cached on conn."""
    if getattr(conn, "_csr_cache", None) is None:
        if hasattr(conn, "_nb_offs"):
            conn._csr_cache = (conn._nb_offs, conn._nb_data)
        else:
            deg = np.array([0 if c is None else len(c.nbr) for c in conn.cells], np.int64)
            offs = np.concatenate([[0], np.cumsum(deg)])
            data = np.array([j for c in conn.cells if c is not None for j in c.nbr], np.int64)
            conn._csr_cache = (offs, data)
    return conn._csr_cache


def _perp_basis(t):
    t = np.asarray(t, float) / np.linalg.norm(t)
    a = np.array([1.0, 0, 0]) if abs(t[0]) < 0.9 else np.array([0, 1.0, 0])
    u = np.cross(t, a); u /= np.linalg.norm(u)
    v = np.cross(t, u)
    return t, u, v


def loop_around_axis(conn, p, t, r_min, r_max, dz, allowed=None, outline=None):
    """Shortest closed walk of kept cells that winds once, counter-clockwise about +t, around the
    axis through p along t, staying in the annulus r_min < r < r_max and the slab |s| < dz.
    Found by BFS in the double cover of the kept-cell adjacency graph (the sheet flips when a step
    crosses the half-plane theta = pi). With `outline` (points of a defect cluster), only cells outside the
    cluster's outline in the slab are allowed: per angular sector, farther from the axis than the cluster's
    farthest point in that sector or its neighbours, less 0.5 spacings. Returns a list of cells, or None."""
    t, u, v = _perp_basis(t)
    # only cells within reach of the annulus and slab (spatial tree, cached on conn)
    near = _near(conn, p, np.hypot(r_max, dz) + 1e-9)
    rel = conn.minimage(conn.pts[near] - np.asarray(p, float))
    s_ = rel @ t
    perp = rel - np.outer(s_, t)
    rn = np.linalg.norm(perp, axis=1)
    sel = conn.kept[near] & (rn > r_min) & (rn < r_max) & (np.abs(s_) < dz)
    if allowed is not None:
        sel &= np.asarray(allowed)[near]
    if outline is not None:
        relo = conn.minimage(np.asarray(outline) - np.asarray(p, float))
        so = relo @ t
        po = relo[np.abs(so) < dz] - np.outer(so[np.abs(so) < dz], t)
        nbin = 16
        rin = np.zeros(nbin)
        if len(po):
            bo = ((np.arctan2(po @ v, po @ u) + np.pi) / (2 * np.pi) * nbin).astype(int) % nbin
            np.maximum.at(rin, bo, np.linalg.norm(po, axis=1))
        rin = np.maximum(rin, np.maximum(np.roll(rin, 1), np.roll(rin, -1)))
        bc = ((np.arctan2(perp @ v, perp @ u) + np.pi) / (2 * np.pi) * nbin).astype(int) % nbin
        sel &= rn > rin[bc] - 0.5 * conn.ref.scale
    if not sel.any():
        return None
    n = len(conn.pts)
    ok = np.zeros(n, bool); ok[near[sel]] = True
    r = np.full(n, np.inf); r[near] = rn
    th = np.zeros(n); th[near] = np.arctan2(perp @ v, perp @ u)
    s = np.full(n, np.inf); s[near] = s_
    rm = 0.5 * (r_min + r_max)
    seed = int(np.argmin(np.where(ok, (r - rm) ** 2 + (rm * th) ** 2 + s ** 2, np.inf)))
    # BFS in the double cover: node 2k + sheet for the k-th allowed cell; a step flips the sheet when it
    # crosses the half-plane theta = pi
    offs, data = _csr(conn)
    sub = np.flatnonzero(ok)
    loc = -np.ones(len(ok), np.int64)
    loc[sub] = np.arange(len(sub))
    deg = offs[sub + 1] - offs[sub]
    a = np.repeat(sub, deg)
    start = np.repeat(offs[sub] - (np.cumsum(deg) - deg), deg)
    b = data[start + np.arange(len(a))]
    keep = (b >= 0)
    a, b = a[keep], b[keep]
    keep = ok[b]
    a, b = a[keep], b[keep]
    flip = (np.abs(th[a] - th[b]) > np.pi).astype(np.int64)
    la, lb = loc[a], loc[b]
    rows = np.concatenate([2 * la, 2 * la + 1])
    cols = np.concatenate([2 * lb + flip, 2 * lb + 1 - flip])
    m = 2 * len(sub)
    G = csr_matrix((np.ones(len(rows), np.int8), (rows, cols)), shape=(m, m))
    order, pred = breadth_first_order(G, 2 * loc[seed], directed=True, return_predecessors=True)
    goal = 2 * loc[seed] + 1
    if pred[goal] < 0:
        return None
    walk, st = [], pred[goal]
    while st >= 0:
        walk.append(int(sub[st // 2]))
        st = pred[st]
    walk = walk[::-1]
    dth = np.diff(th[walk + walk[:1]])
    wind = np.sum((dth + np.pi) % (2 * np.pi) - np.pi)
    if abs(abs(wind) - 2 * np.pi) > 1e-6:
        return None
    if wind < 0:
        walk = walk[:1] + walk[1:][::-1]
    # the walk must not wrap around a periodic direction
    steps = conn.minimage(np.diff(conn.pts[walk + walk[:1]], axis=0))
    if np.abs(steps.sum(0)).max() > 1e-6:
        return None
    return walk


def defect_clusters(conn, min_size=8, mask=None, hops=2):
    """Clusters of cells that are NOT kept, connected through chains of at most `hops` face steps
    (hops=2 merges fragments of a fat, noisy core separated by single perfect cells). Cells touching a
    container wall, or outside `mask`, are ignored. Returns a list of index arrays (bad cells only)."""
    n = len(conn.pts)
    bad = np.zeros(n, bool)
    for i, c in enumerate(conn.cells):
        if c is not None and not conn.kept[i] and len(c.nbr) and min(c.nbr) >= 0:
            bad[i] = True
    if mask is not None:
        bad &= mask
    seen = np.zeros(n, bool)
    out = []
    for i in np.flatnonzero(bad):
        if seen[i]:
            continue
        comp, dq = [], deque([i])
        seen[i] = True
        while dq:
            a = dq.popleft()
            comp.append(a)
            front = {a}
            reach = set()
            for _ in range(hops):
                nxt = set()
                for x in front:
                    for b in conn.cells[x].nbr:
                        if b >= 0 and b not in reach:
                            nxt.add(b)
                reach |= nxt
                front = nxt
            for b in reach:
                if bad[b] and not seen[b]:
                    seen[b] = True
                    dq.append(b)
        if len(comp) >= min_size:
            out.append(np.array(comp))
    return out


def _unwrap(conn, idx):
    """Positions of a cluster, unwrapped across periodic boundaries (BFS over members within a few
    face steps), and the periodic translations through which the cluster connects to its own images
    (non-empty if it wraps around the box, e.g. a straight dislocation along a periodic axis)."""
    members = set(idx.tolist())
    # member pairs within two face steps, and their minimum-image displacements (one vectorised call)
    pairs = []
    for a in members:
        near = set()
        for b in conn.cells[a].nbr:
            if b < 0:
                continue
            near.add(b)
            near.update(c2 for c2 in conn.cells[b].nbr if c2 >= 0)
        pairs.extend((a, b) for b in near if b in members and b != a)
    nbr = {}
    if pairs:
        P = np.array(pairs)
        D = conn.minimage(conn.pts[P[:, 1]] - conn.pts[P[:, 0]])
        for (a, b), dv in zip(pairs, D):
            nbr.setdefault(a, []).append((b, dv))
    start = int(idx[0])
    pos = {start: conn.pts[start].copy()}
    dq = deque([start])
    wraps = []
    while dq:
        a = dq.popleft()
        for b, dv in nbr.get(a, ()):
            q = pos[a] + dv
            if b not in pos:
                pos[b] = q
                dq.append(b)
            else:
                w = q - pos[b]
                if np.abs(w).max() > 1e-6:
                    wraps.append(w)
    keep = [i for i in idx if i in pos]
    return np.array([pos[i] for i in keep]), wraps


def trace_lines(conn, min_cluster=8, spacing=2.5, rho=4.0, annulus=(0.8, 8.0), slab=2.0, mask=None, hops=1,
                retries=True, loops=None, fill=0.5):
    """Trace dislocation lines through the defect tubes. Lengths in units of the reference spacing d.
    Returns a list of dicts, one per cluster: 'cells', 'points' (ordered samples along the tube),
    'tangents', 'b' (per sample: dict from Connection.burgers, or None), 'b_mode' (most common
    lattice vector), 'agreement' (fraction of successful samples with b_mode), 'b_counts' (every lattice
    vector returned on the line, with counts) and 'conflict' (True if two different non-zero vectors were
    returned by samples of the same line and not explained by segments, sums, sense flips or sub-enclosures;
    'sense_flips' and 'sub_enclosures' count the latter two), and 'loops' (per sample: the loop's cells, or None).
    With `retries`, samples whose
    loop did not close are retried with wider annuli, thicker slabs, tilted planes and shifted samples.
    With loops="layer", each loop is the shortest walk around the line through the first layer of the cluster
    (kept cells sharing a face with it; theory note, Lemma "one layer suffices"), tried in slabs of 1, 1.5 and 2
    times `slab`; samples where none closes fall back to the annulus loops (counted in trace_lines.stats)."""
    if loops is None:                 # default from the environment, so that any analysis script can switch
        import os
        loops = os.environ.get("VB_LOOPS", "annulus")
    d = conn.ref.scale
    n_retry, n_rescued = [0], [0]
    n_layer, n_fallback, n_sliver = [0], [0], [0]
    offs, data = _csr(conn)
    out, pending = [], []
    clusters = defect_clusters(conn, min_cluster, mask, hops=hops)
    owner = -np.ones(len(conn.pts), np.int64)
    for ci, idx in enumerate(clusters):
        owner[idx] = ci
    for ci, idx in enumerate(clusters):
        X, wraps = _unwrap(conn, idx)
        wrap_dir = None
        if wraps:                       # the tube connects to its periodic image: that's the line
            w = max(wraps, key=np.linalg.norm)
            wrap_dir = w / np.linalg.norm(w)
        if loops == "layer":
            # fill the cluster's pores (kept cells with at least `fill` of their faces on the cluster, twice), then
            # take its first layer: the kept cells sharing a face with it. Filling only moves loops outwards
            inD = np.zeros(len(conn.pts), bool); inD[idx] = True
            for _ in range(2):
                cand = data[np.concatenate([np.arange(offs[i], offs[i + 1]) for i in np.flatnonzero(inD)])]
                cand = np.unique(cand[(cand >= 0)]); cand = cand[~inD[cand]]
                frac = np.array([inD[data[offs[c]:offs[c + 1]]].mean() for c in cand]) if len(cand) else np.zeros(0)
                inD[cand[frac >= fill]] = True
            nb_ = data[np.concatenate([np.arange(offs[i], offs[i + 1]) for i in np.flatnonzero(inD)])]
            nb_ = nb_[nb_ >= 0]
            layer = np.zeros(len(conn.pts), bool)
            layer[nb_[conn.kept[nb_] & ~inD[nb_]]] = True
        if loops == "layer" and len(X) == len(idx):
            # depth of each cluster cell: face steps to the nearest cell outside the filled cluster. The core of
            # a dislocation is the thickest part of the cluster, thermal blobs attached to it are thin
            depth = np.full(len(conn.pts), -1, np.int64)
            front = [i for i in np.flatnonzero(inD) if (~inD[data[offs[i]:offs[i + 1]][data[offs[i]:offs[i + 1]] >= 0]]).any()]
            depth[front] = 0
            dq = deque(front)
            while dq:
                a_ = dq.popleft()
                for b_ in data[offs[a_]:offs[a_ + 1]]:
                    if b_ >= 0 and inD[b_] and depth[b_] < 0:
                        depth[b_] = depth[a_] + 1; dq.append(b_)
            dX = depth[idx]
        # samples by farthest-point sampling (layer loops: over the core, the cells at least one step deep, so that
        # samples do not sit on thin thermal blobs attached to it)
        def fps(Y):
            samp = [0]
            dist = np.linalg.norm(Y - Y[0], axis=1)
            while dist.max() > spacing * d:
                k = int(dist.argmax())
                samp.append(k)
                dist = np.minimum(dist, np.linalg.norm(Y - Y[k], axis=1))
            return Y[samp]
        S = None
        if loops == "layer" and len(X) == len(idx) and (dX >= 1).sum() >= 4:
            from scipy.spatial import cKDTree
            Sd = fps(X[dX >= 1])
            # the deep cells follow the whole line: they span the cluster, or the line crosses the periodic box,
            # or they give enough samples to trace it (thin cores, e.g. a small loop at moderate T, have few)
            if wrap_dir is not None or len(Sd) >= 6 or cKDTree(X[dX >= 1]).query(X)[0].max() < 1.5 * spacing * d:
                S = Sd
        if S is None:
            S = fps(X)
        # local tangents by PCA
        T = []
        for q in S:
            nb = X[np.linalg.norm(X - q, axis=1) < rho * d]
            w, V = np.linalg.eigh(np.cov((nb - nb.mean(0)).T) if len(nb) > 3 else np.eye(3))
            T.append(V[:, -1] * (1 if w[-1] > 0 else 0))
        T = np.array(T)
        ring = False
        if wrap_dir is not None:
            T[:] = wrap_dir
            order = np.argsort(S @ wrap_dir)
        else:
            # closed loop? the samples encircle their centroid in their best-fit plane, far from it
            c0 = S.mean(0)
            w_, V_ = np.linalg.eigh(np.cov((S - c0).T)) if len(S) > 3 else (np.zeros(3), np.eye(3))
            nrm = V_[:, 0]
            rr = np.linalg.norm((S - c0) - np.outer((S - c0) @ nrm, nrm), axis=1)
            ring = len(S) >= 6 and np.mean(rr > 0.5 * np.median(rr)) >= 0.85 and np.median(rr) > 2 * spacing * d
            if ring:
                # robust ring samples: angular bins of ALL cluster cells, median position per bin
                u_ = V_[:, 2]; v_ = np.cross(nrm, u_)
                angX = np.arctan2((X - c0) @ v_, (X - c0) @ u_)
                nb_ = max(8, int(round(2 * np.pi * np.median(rr) / (spacing * d))))
                edges_ = np.linspace(-np.pi, np.pi, nb_ + 1)
                Sb = []
                for k in range(nb_):
                    m_ = (angX >= edges_[k]) & (angX < edges_[k + 1])
                    if m_.sum() >= 3:
                        Sb.append(np.median(X[m_], axis=0))
                S = np.array(Sb)
                T = np.zeros_like(S)
                order = np.arange(len(S))
            else:
                order = _chain_order(S)
        S, T = S[order], T[order]
        if ring:      # tangents along the ring (central differences), sense: normal has a positive largest component
            T = np.roll(S, -1, axis=0) - np.roll(S, 1, axis=0)
            T /= np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-12)
            if nrm[np.argmax(np.abs(nrm))] < 0:
                S, T = S[::-1], -T[::-1]
        # consistent sense: the first tangent along the chain, then each continuous with the last
        if not ring:
            if len(S) > 1 and T[0] @ (S[1] - S[0]) < 0:
                T[0] = -T[0]
            for k in range(1, len(T)):
                if T[k] @ T[k - 1] < 0:
                    T[k] = -T[k]
        # canonical line sense (xi and -xi with b and -b are the same dislocation): the mean
        # tangent's largest component is positive (rings: fixed above by their normal)
        mt = T.sum(0)
        if not ring and mt[np.argmax(np.abs(mt))] < 0:
            T = -T
            S, T = S[::-1], T[::-1]
        # smooth along the chain (fat, irregular cores at high temperature)
        if len(T) > 2:
            Ts = T.copy()
            for k in range(len(T)):
                Ts[k] = T[max(k - 1, 0)] + 2 * T[k] + T[min(k + 1, len(T) - 1)]
            T = Ts / np.maximum(np.linalg.norm(Ts, axis=1, keepdims=True), 1e-12)
        bs, loops_ = [], []
        for k in range(len(S) if not ring else 0):   # re-centre samples on the local tube centroid
            for _ in range(2):
                m_ = np.linalg.norm(X - S[k], axis=1) < rho * d
                if loops == "layer" and len(X) == len(idx) and m_.any():
                    m_ &= dX >= dX[m_].max() - 1          # the deepest cells nearby: the core, not blobs
                nb = X[m_]
                if len(nb):
                    S[k] = nb.mean(0)
        for q, t in zip(S, T):
            if not np.any(t):
                bs.append(None)
                loops_.append((q, t, None))
                continue
            rel = X - q
            sc = rel @ t
            near = (np.abs(sc) < rho * d / 2) & (np.linalg.norm(rel, axis=1) < rho * d)   # local only
            perp = rel[near] - np.outer(sc[near], t)
            rt = np.linalg.norm(perp, axis=1).max() if near.any() else 0.0
            r0, r1 = rt + annulus[0] * d, rt + annulus[1] * d
            # another cluster within reach: first try a loop that stays closer than it, so that the loop
            # encloses this cluster alone (e.g. one of two dissociated partials); else the full annulus
            lp = None
            if loops == "layer":
                # a shell loop must wind around the whole local core: the cluster cells in the slab (a loop around a
                # sliver of a flat core, e.g. the edge of a partial's ribbon, returns 0 and is rejected); if the
                # shortest one does not, retry with shell cells outside the core's outline, then fall back
                m_ = np.linalg.norm(X - q, axis=1) < rho * d
                if len(X) == len(idx) and m_.any():
                    m_ &= dX >= dX[m_].max() - 1            # the core near q (its deepest cells), not attached blobs
                local = conn.pts[idx[m_]]
                for fz in (1.0, 1.5, 2.0):
                    for outl in (None, local):
                        lp = loop_around_axis(conn, q, t, 0.0, r1, fz * slab * d, allowed=layer, outline=outl)
                        if lp is not None and _enclosed(conn, lp, q, t, local, slab * d) < 0.9:
                            n_sliver[0] += 1; lp = None
                        if lp is not None:
                            break
                    if lp is not None:
                        break
                if lp is not None:
                    lp = _tighten(conn, lp)                   # a geodesic ring: locally shortest, same holonomy
                n_layer[0] += lp is not None; n_fallback[0] += lp is None
            if lp is None and len(clusters) > 1:
                oth = np.flatnonzero((owner >= 0) & (owner != ci))
                relo = conn.minimage(conn.pts[oth] - q)
                so = relo @ t
                inslab = np.abs(so) < 2 * slab * d
                if inslab.any():
                    do = np.linalg.norm(relo[inslab] - np.outer(so[inslab], t), axis=1).min()
                    if rt + 0.5 * d < do < r1:
                        lp = loop_around_axis(conn, q, t, rt, do, slab * d)
            if lp is None:
                lp = loop_around_axis(conn, q, t, r0, r1, slab * d)
            # retry ladder for samples whose loop did not close (wider annulus, thicker slab, tilted plane,
            # shifted sample). Safe: by homotopy invariance any closed loop of kept cells returns the exact
            # content it encloses; retries only change which loops are found, never what a loop returns
            q_used, t_used = q, t
            if lp is None and retries:
                for fr, fz, tilt, shift in _RETRIES:
                    q_used, t_used = q + shift * spacing * d * t, _tilted(t, tilt)
                    lp = loop_around_axis(conn, q_used, t_used, r0, r0 + fr * (r1 - r0), fz * slab * d)
                    if lp is not None:
                        break
                n_retry[0] += 1; n_rescued[0] += lp is not None
            bs.append(None if lp is None else conn.burgers(lp))
            loops_.append((q_used, t_used, lp))
        pending.append(dict(idx=idx, S=S, T=T, bs=bs, loops=loops_, ring=ring))
    _attribute(conn, pending, slab * d)
    for P in pending:
        idx, S, T, bs, ring = P["idx"], P["S"], P["T"], P["bs"], P["ring"]
        keys = [tuple(float(v) for v in np.round(b["lattice"], 3)) for b in bs if b is not None and b["R_identity"]]
        nz = [k for k in keys if any(abs(v) > 1e-6 for v in k)]
        mode = max(set(nz), key=nz.count) if nz else (tuple([0.0] * 3) if keys else None)
        # segments: samples grouped by their (non-zero) Burgers vector and then by proximity, so that a
        # cluster holding e.g. two partial dislocations is reported as segments with their own b
        groups = {}
        for k, b_ in enumerate(bs):
            if b_ is not None and b_["R_identity"]:
                kk = tuple(float(v) for v in np.round(b_["lattice"], 3))
                if any(abs(v) > 1e-6 for v in kk):
                    groups.setdefault(kk, []).append(k)
        segs = []
        for key, ks in groups.items():
            ks = list(ks)
            while ks:
                comp, front = [ks.pop(0)], None
                front = list(comp)
                while front:
                    a_ = front.pop()
                    near = [k for k in ks if np.linalg.norm(conn.minimage(S[k] - S[a_])) < 3 * spacing * d]
                    for k in near:
                        ks.remove(k); comp.append(k); front.append(k)
                segs.append(dict(b=key, samples=sorted(comp), points=S[sorted(comp)]))
        # conflict check: different non-zero vectors on one line are expected only when they form separate
        # segments (e.g. two partials) or when one is the sum of others (a loop enclosing several); anything
        # else would contradict exactness and is flagged
        cnt = {}
        for kk in keys:
            cnt[kk] = cnt.get(kk, 0) + 1
        # a vector whose negative is a segment's vector is a sense flip (the local tangent reversed, e.g. at a sharp
        # bend): b and -b with opposite line senses are the same dislocation. A vector v with k - v equal to v or to
        # another vector of the line, for a segment vector k, is a sub-enclosure: a loop that encircled one of the
        # dislocations of a dissociated segment (e.g. one superpartial of a superdislocation). Both are correct
        # loop contents and are counted, not flagged
        distinct = [np.array(k_) for k_ in set(nz)]
        conflict, flips, subs = False, 0, 0
        if len(distinct) > 1:
            seg_keys = [sg["b"] for sg in segs if len(sg["samples"]) >= 2]
            for v_ in distinct:
                # sums of two segment vectors, the same one twice included (a loop around two parallel
                # dislocations with equal b, e.g. both superpartials of a superdislocation)
                sk = [np.array(k_) for k_ in seg_keys]
                sums = [a_ + b2 for i_, a_ in enumerate(sk) for b2 in sk[i_:]]
                if any(np.allclose(k_, v_) for k_ in seg_keys) or any(np.allclose(s_, v_) for s_ in sums):
                    continue
                if any(np.allclose(np.array(k_), -v_) for k_ in seg_keys):
                    flips += cnt[tuple(float(t_) for t_ in v_)]
                elif any(np.allclose(np.array(k_) - v_, v_) or any(np.allclose(np.array(k_) - v_, w_) for w_ in distinct)
                         for k_ in seg_keys):
                    subs += cnt[tuple(float(t_) for t_ in v_)]
                else:
                    conflict = True
        out.append(dict(cells=idx, points=S, tangents=T, b=bs, b_mode=mode, ring=bool(ring), segments=segs,
                        b_counts={str(list(k_)): v for k_, v in cnt.items()}, conflict=conflict, sense_flips=flips, sub_enclosures=subs,
                        loops=[lp_ for (_, _, lp_) in P["loops"]],
                        agreement=(nz.count(mode) / len(nz)) if nz else 0.0,
                        zero_fraction=(1 - len(nz) / len(keys)) if keys else 0.0, n_samples=len(S)))
    trace_lines.stats = dict(retried=n_retry[0], rescued=n_rescued[0], layer=n_layer[0], layer_fallback=n_fallback[0],
                             layer_rejected=n_sliver[0])
    return out


# retry ladder: (annulus width factor, slab factor, tilt in degrees, shift along the line in sample spacings)
_RETRIES = [(1.5, 1.0, 0, 0), (2.0, 1.5, 0, 0), (1.0, 1.0, 0, 0.5), (1.0, 1.0, 0, -0.5),
            (1.5, 1.0, 12, 0), (1.5, 1.0, -12, 0), (2.5, 2.0, 0, 0)]


def _tighten(conn, lp):
    """Shorten a closed walk of kept cells by local moves until none applies: remove back-and-forth steps
    (a -> b -> a) and replace a -> b -> c by a -> c when a and c are neighbours and the elementary loop a, b, c has
    trivial holonomy (by flatness, every elementary loop of kept cells around a Voronoi edge does), so the walk's
    holonomy is unchanged. The result is a locally shortest (geodesic) ring in its class."""
    w = list(lp)
    changed = True
    while changed and len(w) > 3:
        changed = False
        n = len(w)
        for i in range(n):
            a, b, c = w[i - 1], w[i], w[(i + 1) % n]
            if a == c:                                              # a -> b -> a
                w = [x for k, x in enumerate(w) if k not in (i, (i + 1) % n)]
                changed = True; break
            if c in conn.cells[a].nbr:
                h = conn.holonomy([a, b, c])
                if h is not None and np.abs(h[0] - np.eye(3)).max() < 1e-6 and np.abs(h[1]).max() < 1e-6:
                    del w[i]
                    changed = True; break
    return w


def _enclosed(conn, lp, q, t, pts, dz):
    """Fraction of the points within dz of the plane through q (normal t) about which the loop, projected on that
    plane, winds non-trivially."""
    t, u, v = _perp_basis(t)
    X = conn.minimage(np.asarray(pts) - q)
    X = X[np.abs(X @ t) < dz]
    if not len(X):
        return 1.0
    P = conn.minimage(conn.pts[lp] - q)
    P, X = np.c_[P @ u, P @ v], np.c_[X @ u, X @ v]
    w = np.zeros(len(X))
    for a, b in zip(P, np.roll(P, -1, axis=0)):
        aa, bb = a[None] - X, b[None] - X
        w += np.arctan2(aa[:, 0] * bb[:, 1] - aa[:, 1] * bb[:, 0], (aa * bb).sum(1))
    return float(np.mean(np.abs(w) > np.pi))


def _tilted(t, deg):
    """t rotated by deg degrees about a fixed axis perpendicular to it."""
    if not deg:
        return t
    t, u, v = _perp_basis(t)
    a = np.radians(deg)
    return np.cos(a) * t + np.sin(a) * u


def _winding(poly, pts):
    """Winding numbers of the closed 2D polygon `poly` around each of `pts`."""
    a = poly[None, :, :] - pts[:, None, :]
    b = np.roll(a, -1, axis=1)
    ang = np.arctan2(a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0], (a * b).sum(-1))
    return np.rint(ang.sum(1) / (2 * np.pi)).astype(int)


def _attribute(conn, pending, half_slab):
    """A loop placed around one cluster may also enclose another (a thermal blob next to the dislocation,
    or the other partial of a dissociated dislocation). Holonomies add, so:
      - a loop enclosing only clusters whose own vector is zero keeps its vector;
      - a loop enclosing (nearly all of) another cluster with a non-zero vector has that vector subtracted;
      - a loop enclosing only part of such a cluster cannot be attributed, and is dropped.
    A cluster's own vector is the most common non-zero vector of its attributable loops (zero if there is
    none); clusters are resolved iteratively."""
    if len(pending) < 2:
        return
    enc = []            # per cluster, per sample: {other cluster: fraction of its nearby traced points enclosed}
    for a, P in enumerate(pending):
        row = []
        for q, t, lp in P["loops"]:
            hit = {}
            if lp is not None:
                t_, u, v = _perp_basis(t)
                poly = conn.minimage(conn.pts[lp] - q)
                poly2 = np.column_stack([poly @ u, poly @ v])
                for c, Q in enumerate(pending):
                    if c == a:
                        continue
                    rel = conn.minimage(Q["S"] - q)
                    sl = np.abs(rel @ t_) < 2 * half_slab
                    if not sl.any():
                        continue
                    f = float(np.mean(_winding(poly2, np.column_stack([rel[sl] @ u, rel[sl] @ v])) != 0))
                    if f > 0:
                        hit[c] = f
            row.append(hit)
        enc.append(row)

    def lat(b):
        return np.asarray(b["lattice"], float)

    def attributed(b, h, own):
        """The loop's own vector, None if it cannot be attributed yet (or ever). Conservative: a loop that
        returns zero stays zero, and a subtraction is accepted only if it shortens the vector (the enclosed
        cluster explains part of the loop's holonomy); anything else is left unattributed."""
        raw = lat(b)
        if not np.any(np.abs(raw) > 1e-6):
            return raw
        sub = np.zeros(3)
        for c, f in h.items():
            if own[c] is None:
                return None
            if np.any(np.abs(own[c]) > 1e-6):
                if f < 0.9:
                    return None
                sub = sub + own[c]
        if not np.any(np.abs(sub) > 1e-6):
            return raw
        v = raw - sub
        return v if np.linalg.norm(v) < np.linalg.norm(raw) - 1e-6 else None

    own = [None] * len(pending)
    for _ in range(4):
        changed = False
        for a, P in enumerate(pending):
            if own[a] is not None:
                continue
            vals, waiting = [], False
            for b, h in zip(P["bs"], enc[a]):
                if b is None or not b["R_identity"]:
                    continue
                if any(own[c] is None for c in h):
                    waiting = True
                    continue
                v_ = attributed(b, h, own)
                if v_ is not None:
                    vals.append(tuple(np.round(v_, 3)))
            nz = [v for v in vals if any(abs(x) > 1e-6 for x in v)]
            if nz or (vals and not waiting):
                own[a] = np.array(max(set(nz), key=nz.count) if nz else (0.0, 0.0, 0.0))
                changed = True
        if not changed:
            break
    # still unresolved: clusters whose loops enclose each other (e.g. two partials too close for separate
    # loops). Group them; the largest member reports the group's total through its own loops, which enclose
    # the others, and the loops of the other members are dropped.
    unres = [a for a in range(len(pending)) if own[a] is None]
    grp = {a: a for a in unres}

    def root(a):
        while grp[a] != a:
            a = grp[a]
        return a
    for a in unres:
        for h in enc[a]:
            for c in h:
                if c in grp:
                    grp[root(c)] = root(a)
    groups = {}
    for a in unres:
        groups.setdefault(root(a), []).append(a)
    for members in groups.values():
        big = max(members, key=lambda a: len(pending[a]["idx"]))
        for a in members:
            if a != big:
                pending[a]["bs"] = [None] * len(pending[a]["bs"])
                own[a] = np.zeros(3)
        ks = [tuple(np.round(lat(b), 3)) for b in pending[big]["bs"] if b is not None and b["R_identity"]]
        nz = [k for k in ks if any(abs(x) > 1e-6 for x in k)]
        own[big] = np.array(max(set(nz), key=nz.count) if nz else (0.0, 0.0, 0.0))
        for k, h in enumerate(enc[big]):        # the group's members are already counted in big's loops
            enc[big][k] = {c: f for c, f in h.items() if c not in members}
    for a, P in enumerate(pending):
        for k, (b, h) in enumerate(zip(P["bs"], enc[a])):
            if b is None or not h:
                continue
            v_ = attributed(b, h, own)
            if v_ is None:
                P["bs"][k] = None
                continue
            if np.allclose(v_, lat(b)):
                continue
            b = dict(b)
            b["lattice"] = v_
            b["encloses_other"] = len(h)
            P["bs"][k] = b


def _chain_order(S):
    """Order sample points along a curve: minimum spanning tree, then the longest path through it
    (remaining points attached at their nearest position)."""
    n = len(S)
    if n < 3:
        return np.arange(n)
    D = np.linalg.norm(S[:, None] - S[None], axis=2)
    inT = np.zeros(n, bool); inT[0] = True
    best = D[0].copy(); parent = np.zeros(n, int)
    adj = {i: [] for i in range(n)}
    for _ in range(n - 1):
        cand = np.where(inT, np.inf, best)
        k = int(cand.argmin())
        inT[k] = True
        adj[k].append(parent[k]); adj[parent[k]].append(k)
        upd = (D[k] < best) & ~inT
        best[upd] = D[k][upd]; parent[upd] = k

    def far(a):
        dist = {a: 0.0}; prev = {a: None}; dq = deque([a])
        while dq:
            x = dq.popleft()
            for y in adj[x]:
                if y not in dist:
                    dist[y] = dist[x] + D[x, y]; prev[y] = x; dq.append(y)
        z = max(dist, key=dist.get)
        return z, prev

    a, _ = far(0)
    b, prev = far(a)
    path = []
    x = b
    while x is not None:
        path.append(x); x = prev[x]
    rest = [i for i in range(n) if i not in set(path)]
    for i in rest:                       # insert remaining points next to their nearest path point
        j = int(np.argmin(D[i, path]))
        path.insert(j + 1, i)
    return np.array(path)
