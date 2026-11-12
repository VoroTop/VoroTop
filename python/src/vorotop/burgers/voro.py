"""Voronoi tessellation: faces of every cell, via the small voro++ driver `voro_faces`.

The driver (tools/voro_faces.cc) prints, per particle, its faces as cycles of local vertex
indices (consistently oriented by voro++) and the neighbour across each face. Build it with
    g++ -O2 -fopenmp -I<voro++>/src voro_faces.cc -L<voro++>/src -lvoro++ -o voro_faces
It is built and installed with the package (bin/voro_faces); $VOROTOP_BURGERS_VORO_FACES overrides it.
"""

import os
import subprocess
import tempfile

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
VORO_FACES = os.environ.get("VOROTOP_BURGERS_VORO_FACES", os.path.join(_HERE, "bin", "voro_faces"))   # built and installed with the package


def default_threads():
    """Threads for voro++ and the engine: VB_THREADS if set; otherwise every CPU of a SLURM job; otherwise 4
    (so that interactive use on a shared machine stays modest)."""
    if os.environ.get("VB_THREADS"):
        return int(os.environ["VB_THREADS"])
    if os.environ.get("SLURM_CPUS_PER_TASK"):
        return len(os.sched_getaffinity(0))
    return 4


def faces(pts, L, periodic, threads=None):
    """Faces of every cell in the box [0, L) (periodic flags per axis).
    Returns a list (one entry per particle) of [(neighbour, [local vertex ids]), ...], or None for
    a cell voro++ could not compute. Negative neighbours are container walls."""
    pts = np.asarray(pts, float)
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        L = np.asarray(L, float)
        if L.shape == (3, 3):                       # fully periodic restricted-triclinic box, rows a, b, c
            f.write(f"{len(pts)} T {L[0,0]} {L[1,0]} {L[1,1]} {L[2,0]} {L[2,1]} {L[2,2]}\n")
        else:
            f.write(f"{len(pts)} {L[0]} {L[1]} {L[2]} {int(periodic[0])} {int(periodic[1])} {int(periodic[2])}\n")
        np.savetxt(f, pts, fmt="%.10f")
        name = f.name
    try:
        env = dict(os.environ, OMP_NUM_THREADS=str(threads or default_threads()))
        out = subprocess.run([VORO_FACES, name], capture_output=True, text=True, env=env, check=True).stdout
    finally:
        os.unlink(name)
    cells = [None] * len(pts)
    for line in out.splitlines():
        t = list(map(int, line.split()))
        i, nf, k = t[0], t[1], 2
        if nf == 0:
            continue
        fs = []
        for _ in range(nf):
            nb, m = t[k], t[k + 1]
            fs.append((nb, t[k + 2:k + 2 + m]))
            k += 2 + m
        cells[i] = fs
    return cells
