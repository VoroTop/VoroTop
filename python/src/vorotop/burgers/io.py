"""Input: LAMMPS text dumps."""

import numpy as np


def read_lammps_dump(fn):
    """Yield (timestep, box, cols) per frame; box is 3x2 [lo, hi], cols a dict of numpy arrays keyed by column name,
    sorted by id. Triclinic dumps (BOX BOUNDS xy xz yz) are converted from the bounding box to the true lo/hi with the
    LAMMPS formulas, and cols["H"] holds the box as a restricted-triclinic matrix (rows a = (lx,0,0), b = (xy,ly,0),
    c = (xz,yz,lz)); see box_matrix. Files ending in .gz are read directly."""
    import gzip
    opener = (lambda p: gzip.open(p, "rt")) if str(fn).endswith(".gz") else open
    with opener(fn) as f:
        while True:
            line = f.readline()
            if not line:
                return
            assert line.startswith("ITEM: TIMESTEP"), line
            step = int(f.readline())
            f.readline()
            n = int(f.readline())
            hdr = f.readline()
            rows = [list(map(float, f.readline().split())) for _ in range(3)]
            box = np.array([r[:2] for r in rows])
            tri = "xy xz yz" in hdr
            if tri:
                xy, xz, yz = rows[0][2], rows[1][2], rows[2][2]
                box[0] -= [min(0.0, xy, xz, xy + xz), max(0.0, xy, xz, xy + xz)]
                box[1] -= [min(0.0, yz), max(0.0, yz)]
            names = f.readline().split()[2:]
            a = np.loadtxt([f.readline() for _ in range(n)], ndmin=2)
            a = a[np.argsort(a[:, names.index("id")])]
            cols = {k: a[:, j] for j, k in enumerate(names)}
            cols["periodic"] = tuple(p == "pp" for p in hdr.split()[-3:]) if "BOUNDS" in hdr else (True,) * 3
            if tri:
                L = box[:, 1] - box[:, 0]
                cols["H"] = np.array([[L[0], 0.0, 0.0], [xy, L[1], 0.0], [xz, yz, L[2]]])
            yield step, box, cols


def box_matrix(box, cols):
    """The box to pass to the engine: the restricted-triclinic matrix for a triclinic dump, else the edge lengths."""
    return cols["H"] if "H" in cols else box[:, 1] - box[:, 0]


def positions(cols):
    for keys in (("x", "y", "z"), ("xu", "yu", "zu"), ("xs", "ys", "zs")):
        if all(k in cols for k in keys):
            return np.stack([cols[k] for k in keys], 1)
    raise KeyError("no position columns")
