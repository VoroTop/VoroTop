# VoroTop in Python: Burgers vectors

The Python interface of [VoroTop](https://www.vorotop.org). Its first analysis, `vorotop.burgers`, finds the
dislocations in an atomistic configuration and returns their exact Burgers vectors. The only input besides the atoms
is the ideal crystal structure: there are no structure templates, no distance cutoffs and nothing to tune. It reads
raw thermal snapshots as they are, so every frame of a molecular-dynamics trajectory can be analyzed.

How it works, in brief: each atom's Voronoi cell is matched, face by face, to the ideal crystal's cell; neighboring
matched cells are related by exact crystal symmetries; the atoms whose cells do not match form defect clusters; and
the Burgers vectors of each cluster are read from closed loops through the layer of crystal that encases it. See the
paper for details and guarantees.

> E. A. Lazar, *Burgers vectors from Voronoi topology*, [arXiv:2611.22488](https://arxiv.org/abs/2611.22488) (2026).

## Install

From this folder (`python/` of the VoroTop repository):

```
pip install .
```

This builds the bundled [voro++](https://math.lbl.gov/voro++/) and the package's C++ programs (the analysis engine,
the close-packed route, the shortest-circuit search, and a tessellation tool; a C++17 compiler with OpenMP is needed)
and installs the Python package with its dependencies (numpy, scipy, spglib, pymatgen).

## Use

```
vorotop-burgers dump.lammpstrj.gz --cif crystal.cif --types 1=Al,2=Cu
```

- `dump...`: one or more LAMMPS text dumps, multi-frame, optionally gzipped, orthogonal or triclinic. The box must
  be periodic in all three directions; every atom of every frame is analyzed.
- `--cif`: the ideal crystal structure. It is snapped onto its space group (the one given with `--spacegroup N`,
  otherwise the one spglib finds), and the output says how.
- `--types`: the element of each LAMMPS atom type; needed when the crystal has several species.
- `--orient "1 1 -2, 1 1 1, 1 -1 0"`: optional, for a single-crystal sample whose orientation is known: the crystal
  directions along the box's x, y and z axes, in lattice indices of the CIF cell (as LAMMPS `orient`).
- `--frames START:STOP:STEP`, `--out PREFIX`, `--threads N`: optional.

Outputs:

- `PREFIX.clusters.jsonl`: one record per frame, with the fraction of atoms whose cells are crystalline, any warnings, and
  every defect cluster that carries a Burgers vector: its number of atoms, the genus of its boundary surface, and its
  vectors. Each vector is given in the laboratory frame (`lab`, in the units of the dump) and as its family
  (`family`, e.g. `1/2<110>`, in the conventional cell of the space group; hexagonal crystals in three-index notation).
  With `--orient`, it is also given exactly in the CIF cell (`cell`). Without an orientation the family is all the
  atoms determine: rotating a crystal by one of its own symmetries leaves the same atoms, so which member of the
  family a vector is depends on how the sample's axes were labeled. Every vector is given up to sign: a cluster
  carries no preferred line direction. Every connected region of crystalline cells (a grain, or crystal enclosed by
  liquid) is analyzed in its own frame, and `region` says which one a vector was read in; vectors from different
  regions are compared in the laboratory frame.
- `PREFIX.atoms.dump`: the frames again, with per-atom columns `crystalline`, `cluster` and `carries`, for viewing in
  OVITO.

From Python:

```python
from vorotop.burgers import api, io
crystal = api.crystal_from_cif("crystal.cif")
step, box, cols = next(io.read_lammps_dump("dump.lammpstrj.gz"))
result = api.analyze(io.positions(cols) - box[:, 0], io.box_matrix(box, cols), crystal, species=None)
for c in result["clusters"]:
    if c["vectors"]:
        print(c["n_atoms"], [(v["family"], v["lab"]) for v in c["vectors"]])
```

## Which crystals

- **Crystals whose Voronoi cells are topologically stable** (every Voronoi vertex shared by exactly four cells): bcc,
  B2, Laves phases, θ-Al₂Cu and many more. The results are exact. The relative Voronoi gap g of the crystal is
  printed; above about 0.2 the method is robust at room temperature, and below it thermal noise and strain blur the
  topology (wider cores, merged or missing lines, never wrong vectors).
- **fcc and hcp metals**, through shipped tables of all perturbed cell topologies (hcp at the c/a ratio of Mg and at
  the ideal ratio). Orthogonal boxes only for now. A dissociated dislocation and its stacking fault form one cluster
  carrying the full vector.
- **Other crystals with degenerate cells**, such as ω-Ti: every cell topology that small perturbations can produce
  becomes a reference cell of its own. The enumeration needs [VoroTypes](https://github.com/VoroTop/VoroTypes)
  installed (or `$VOROTYPES` pointing to a source checkout).

## Run time and memory

With 16 threads (`--threads 16`) on current Xeon nodes, a snapshot takes about 6--10 s per million atoms for bcc
metals and theta-Al2Cu, 8--11 s for fcc Cu, and 13--18 s for hcp Mg, with about 1.4 GB of memory per million atoms
(up to 2 GB for hcp metals and dislocation networks in complex crystals). Times grow in proportion to the number of
atoms, tested up to 79 million. More than 16 threads currently gives little or no further gain. They grow faster where one defect cluster holds a large part of the sample: close to melting, where the
non-crystalline atoms percolate, and around large faults attached to dislocation networks.

## Warnings the output can carry

- the crystal's Voronoi gap is below 0.2, so thermal noise and strain may destroy the ideal cell topology (expect wide
  cores and merged or missing lines, never wrong vectors);
- the box is too short (less than three times the longest neighbor distance) in some direction, so a smallest loop can
  wrap around it;
- some cells were labeled by geometric fit where topology could not decide;
- some smallest loops do not close (the warning gives how many, and the largest net translation among them);
- the crystal is topologically pseudo-symmetric, and the best labeling of some region beat the runner-up by less than
  0.01 Å of rms misfit, so that region's labeling (and the vectors read from it) is nearly a tie and worth checking.

## License

BSD 3-clause (see `LICENSE`). The bundled voro++ is distributed under its own license (`extern/voro/LICENSE`).
