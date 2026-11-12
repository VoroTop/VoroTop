"""vorotop-burgers: Burgers vectors of dislocations from Voronoi topology.

    vorotop-burgers DUMP [DUMP ...] --cif CRYSTAL.cif [--spacegroup N] [--types 1=Al,2=Cu] [--out PREFIX]
                    [--frames START:STOP:STEP] [--threads N] [--orient "1 1 -2, 1 1 1, 1 -1 0"]

Every frame of every LAMMPS text dump (multi-frame, optionally gzipped, orthogonal or triclinic, fully periodic) is
analyzed whole. The crystal is given as a CIF file of the ideal structure; it is snapped onto its space group (the one
given with --spacegroup, else the one spglib finds), and the output says so. Outputs (PREFIX defaults to the first
dump's name):
  PREFIX.clusters.jsonl  one record per frame: crystalline fraction, warnings, and every defect cluster that carries Burgers
                         vectors (atoms, genus; each vector in laboratory coordinates and as its family,
                         and, with --orient, exactly in the cell of the CIF)
  PREFIX.atoms.dump      the frames again, with per-atom columns crystalline (0/1), cluster (id, -1 if crystalline) and
                         carries (1 if the atom's cluster carries a Burgers vector), for OVITO
Nothing needs tuning; see the paper for what the results mean."""
import argparse
import json
import sys

import numpy as np

from . import api, io


def _types(spec):
    out = {}
    if spec:
        from pymatgen.core import Element
        for item in spec.split(","):
            t, sym = item.split("=")
            out[int(t)] = Element(sym.strip()).Z
    return out


def _orient(spec):
    if not spec:
        return None
    rows = [[float(x) for x in part.split()] for part in spec.split(",")]
    if len(rows) != 3 or any(len(r) != 3 for r in rows):
        raise argparse.ArgumentTypeError('--orient needs three directions, e.g. "1 1 -2, 1 1 1, 1 -1 0"')
    return rows


def _frames(spec):
    if not spec:
        return None
    parts = [int(x) if x else None for x in spec.split(":")]
    return slice(*parts)


def _write_dump(f, step, box, cols, res):
    n = len(res["crystalline"])
    carry = np.zeros(n, int)
    for c in res["clusters"]:
        if c["vectors"]:
            carry[c["atoms"]] = 1
    tri = "H" in cols
    f.write(f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{n}\n")
    if tri:
        H = cols["H"]; xy, xz, yz = H[1, 0], H[2, 0], H[2, 1]
        lo, hi = box[:, 0].copy(), box[:, 1].copy()
        lo_b = [lo[0] + min(0.0, xy, xz, xy + xz), lo[1] + min(0.0, yz), lo[2]]
        hi_b = [hi[0] + max(0.0, xy, xz, xy + xz), hi[1] + max(0.0, yz), hi[2]]
        f.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
        f.write(f"{lo_b[0]} {hi_b[0]} {xy}\n{lo_b[1]} {hi_b[1]} {xz}\n{lo_b[2]} {hi_b[2]} {yz}\n")
    else:
        f.write("ITEM: BOX BOUNDS pp pp pp\n" + "".join(f"{lo} {hi}\n" for lo, hi in box))
    x = io.positions(cols)
    typ = cols.get("type", np.ones(n))
    ids = cols.get("id", np.arange(1, n + 1))
    f.write("ITEM: ATOMS id type x y z crystalline cluster carries\n")
    for i in range(n):
        f.write(f"{int(ids[i])} {int(typ[i])} {x[i, 0]:.6f} {x[i, 1]:.6f} {x[i, 2]:.6f} "
                f"{int(res['crystalline'][i])} {int(res['cluster'][i])} {carry[i]}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="vorotop-burgers", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dumps", nargs="+", help="LAMMPS text dumps (multi-frame, .gz allowed)")
    ap.add_argument("--cif", required=True, help="the ideal crystal structure")
    ap.add_argument("--spacegroup", type=int, help="declared space-group number (default: as found by spglib)")
    ap.add_argument("--types", help="LAMMPS atom types as elements, e.g. 1=Al,2=Cu (needed for several species)")
    ap.add_argument("--out", help="output prefix (default: the first dump's name)")
    ap.add_argument("--frames", help="frames to analyze, START:STOP:STEP over all dumps (default: all)")
    ap.add_argument("--threads", type=int, help="threads for the tessellation and the engine")
    ap.add_argument("--orient", type=_orient, help="crystal directions along the box's x, y and z axes, in lattice "
                    "indices of the CIF cell (as LAMMPS orient), e.g. \"1 1 -2, 1 1 1, 1 -1 0\"; then every vector "
                    "is also given exactly in the CIF cell")
    a = ap.parse_args(argv)
    crystal = api.crystal_from_cif(a.cif, spacegroup=a.spacegroup)
    if a.orient is not None:
        try:
            crystal.axes(a.orient)
        except ValueError as e:
            sys.exit(f"--orient: {e}")
    print(f"crystal: {crystal.spacegroup}, {crystal.snap_note}; relative Voronoi gap {max(crystal.gap, 0.0):.3g}", file=sys.stderr)
    if crystal.kind in ("fcc", "hcp"):
        print(f"close-packed route ({crystal.kind}): first neighbors from the table of perturbed cell topologies",
              file=sys.stderr)
    elif crystal.kind == "degenerate":
        print("degenerate cells: every perturbed cell topology is a reference cell of its own (VoroTypes)", file=sys.stderr)
    types = _types(a.types)
    multi = len(set(np.asarray(crystal.numbers).tolist())) > 1
    if multi and not types:
        sys.exit("the crystal has several species: give --types, e.g. --types 1=Al,2=Cu")
    prefix = a.out or a.dumps[0].rsplit(".gz", 1)[0].rsplit(".", 1)[0]
    sel = _frames(a.frames)
    k = -1
    with open(prefix + ".clusters.jsonl", "w") as fj, open(prefix + ".atoms.dump", "w") as fd:
        for fn in a.dumps:
            for step, box, cols in io.read_lammps_dump(fn):
                k += 1
                if sel is not None and k not in range(*sel.indices(10 ** 9)):
                    continue
                if not all(cols["periodic"]):
                    sys.exit(f"{fn}, step {step}: the box must be periodic in all three directions")
                x = io.positions(cols) - box[:, 0]
                species = None
                if multi:
                    species = np.array([types[int(t)] for t in cols["type"]])
                res = api.analyze(x, io.box_matrix(box, cols), crystal, species=species, threads=a.threads,
                                  orient=a.orient)
                rec = dict(file=fn, step=int(step), n_atoms=len(x), crystalline_fraction=round(res["crystalline_fraction"], 6),
                           n_clusters=len(res["clusters"]), warnings=res["warnings"], crystal=res["crystal"],
                           clusters=[dict(id=c["id"], n_atoms=c["n_atoms"], genus=c["genus"],
                                          vectors=[{k: v[k] for k in ("family", "lab", "length", "cell", "region") if k in v}
                                                   for v in c["vectors"]])
                                     for c in res["clusters"] if c["vectors"]])
                fj.write(json.dumps(rec) + "\n"); fj.flush()
                _write_dump(fd, step, box, cols, res)
                print(f"{fn} step {step}: {len(rec['clusters'])} clusters carry Burgers vectors; "
                      f"{100 * res['crystalline_fraction']:.1f}% crystalline", file=sys.stderr)


if __name__ == "__main__":
    main()
