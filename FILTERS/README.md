# Filter Files

This directory contains pre-built filter files for common crystal structures that can be used with VoroTop's `-f` option.

## Included Filters

| File | Structure | Description |
|------|-----------|-------------|
| `BCC.filter` | Body-centered cubic | Truncated octahedron topology (6 square + 8 hexagonal faces) |
| `FCC.filter` | Face-centered cubic | 6,294 topological types, excludes HCP types |
| `HCP.filter` | Hexagonal close-packed | 21,611 topological types, excludes FCC types |
| `FCC-both-HCP.filter` | FCC and HCP | 26,530 types: FCC-only, shared FCC/HCP, and HCP-only |
| `ICOS.filter` | Icosahedral | Pentagonal dodecahedron topology (12 pentagonal faces) |
| `crystal_2d.filter` | 2D hexagonal crystal | Single type: hexagon surrounded by hexagons |
| `crystal-and-defects_2d.filter` | 2D crystal with defects | Five types: crystal, grain boundary, dislocation, vacancy, interstitial |
| `fcc-perfect.filter` | Perfect FCC | Single topology of the ideal (unperturbed) FCC Voronoi cell |

## Usage

Specify a filter file when running VoroTop using the `-f` option:

```bash
VoroTop input.dump -f FILTERS/FCC.filter -l      # classify particles and output LAMMPS dump
VoroTop input.dump -f FILTERS/BCC.filter -c       # cluster analysis using BCC filter
```

## Filter File Format

Filter files are plain text with three types of lines:

### Comment lines

Lines beginning with `#` are comments and are ignored.

```
#   This is a comment describing the filter
```

### Structure type lines

Lines beginning with `*` define the structure types used in the filter. Each line contains `*`, followed by an integer index and a plain-text name. Types must be numbered consecutively starting from 1.

```
*   1   Crystal
*   2   GrainBoundary
*   3   Dislocation
```

### Topology entry lines

All other lines associate a Voronoi topology with a structure type. Each line begins with the structure type index, followed by a topology vector enclosed in parentheses.

In **two dimensions**, topology vectors are p-vectors describing the arrangement of neighbor edge counts around a particle. In **three dimensions**, topology vectors are Weinberg vectors encoding the complete combinatorial structure of the Voronoi cell.

```
1   (6,6,6,6,6,6,6)
2   (5,6,6,7,6,7)
```

### Complete example

Here is a minimal 2D filter with two structure types:

```
#   Simple 2D crystal filter
*   1   Crystal
*   2   GrainBoundary
1   (6,6,6,6,6,6,6)
2   (5,6,6,7,6,7)
2   (7,5,6,6,5,6,6,6)
```

Particles whose topology matches a listed entry are assigned the corresponding structure type. Particles with topologies not in the filter are assigned type 0 (unclassified), which typically indicates defects.

## Creating Your Own Filters

VoroTop can generate filter files automatically from crystal structures using the `-mf` option:

```bash
VoroTop crystal.dump -mf                          # exact topologies only
VoroTop crystal.dump -mf 1000000 0.02             # with Gaussian perturbations
```

The perturbation mode replicates the unit cell into a supercell, applies random displacements, and collects all observed topologies into a filter. This produces filters that are robust to thermal vibrations. A target sample count of 1,000,000 is recommended for thorough sampling.

Additional filter files for various structures can be found at [vorotop.org](https://www.vorotop.org).
