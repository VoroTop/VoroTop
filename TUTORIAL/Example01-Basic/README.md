# Example 1: Basic Voronoi topology calculations

**Learning objectives:**

* After completing this example, you will be able to:
    * Run *VoroTop* to calculate the Voronoi topologies of individual particles in two and three dimensions.
    * Understand the mathematical concepts behind Voronoi cells and their topological descriptions.
    * Interpret the computed Voronoi topology data including p-vectors, Weinberg vectors, symmetry, and chirality information.
    * Analyze the structural meaning of different Voronoi topologies in ordered and disordered systems.
    * Use Voronoi topology as a foundation for more advanced structural analysis.
<br><br>

## Background:

The Voronoi cell of a particle is the region of space closer to that particle than to any other particle in the system. This fundamental geometric concept provides a natural way to define neighborhoods and analyze local structure in particle systems. In two dimensions, Voronoi cells are polygons; in three dimensions, they are polyhedra. The *topology* of a Voronoi cell describes not just its shape, but the combinatorial arrangement of its neighbors—information that remains invariant under rotations, translations, and small perturbations.

**Voronoi topology in two dimensions:**
In 2D systems, we characterize each particle by its number of neighbors and the number of neighbors that each of its neighbors has, arranged in a canonical sequence called a *p-vector*. For example, a particle in a perfect hexagonal crystal has six neighbors, each of which also has six neighbors, giving the p-vector (6,6,6,6,6,6,6). The first number counts the central particle's neighbors, and subsequent numbers count the neighbors of each neighboring particle in order.

**Voronoi topology in three dimensions:**
In 3D systems, the topology is more complex and is described by *Weinberg vectors*—sequences that completely characterize the three-dimensional arrangement of faces, edges, and vertices of the Voronoi polyhedron. These vectors capture not only the number of faces and their shapes, but also how those faces are connected to each other.

**Why topology matters:**
Topological descriptions are remarkably robust—they typically don't change under small perturbations of particle positions, making them ideal for analyzing systems with thermal vibrations, measurement errors, or other sources of noise. This robustness, combined with their ability to distinguish different types of local structural environments, makes Voronoi topology a powerful tool for automated structure analysis.
<br><br>


----

**Input data:**

* **Files:** `IdealGas2D.dump`, `IdealGas3D.dump`
* **Format:** LAMMPS dump format with system dimensions, boundary conditions, and particle coordinates
* **Description:** 
  - `IdealGas2D.dump`: 1000 particles uniformly distributed in a unit square, representing a two-dimensional ideal gas
  - `IdealGas3D.dump`: 1000 particles uniformly distributed in a unit cube, representing a three-dimensional ideal gas
* **System characteristics:** These random configurations provide examples of the diverse Voronoi topologies that arise in disordered systems, contrasting with the uniform topologies found in crystalline systems.
<br><br>

**Steps:**

1.  **Navigate to the example directory:** Open your terminal or command prompt and navigate to the example directory within the *VoroTop* tutorial repository.

    ```bash
    cd VoroTop/Tutorial/Example01-Basic
    ```

2.  **Execute the *VoroTop* commands:** Run the following *VoroTop* commands to perform the analysis:

    ```bash
    VoroTop IdealGas3D.dump -vt
    VoroTop IdealGas2D.dump -2 -vt
    ```

    **Command explanation:**
    - First argument: Input data file name
    - `-vt`: Compute and output Voronoi topology vectors for each particle
    - `-2`: Specify two-dimensional analysis (only for 2D systems)

3.  **Examine the output files:** Files named `IdealGas2D.vectors` and `IdealGas3D.vectors` will be created in this directory. These contain the complete Voronoi topology information for each particle.
<br><br>

----

## Expected output and explanation:

### Two-dimensional output (`IdealGas2D.vectors`)

The file contains one line per particle with the following structure:

```text
1     7     (0,1,3,0,1,1,0,1)    (7,4,8,7,5,5,5,10)      1      1
2     5     (0,0,3,0,1,0,1)      (5,5,5,7,5,9)           1      1
3     8     (0,1,2,2,2,1)        (8,4,6,7,5,6,5,8,7)     1     -1
4     5     (0,0,1,2,2)          (5,5,6,7,6,7)           1      1
5     5     (0,0,0,2,1,2)        (5,6,7,6,8,8)           2      0
```

**Column descriptions:**
* **Column 1 - Particle ID:** Unique identifier from the input file
* **Column 2 - Edge count:** Number of edges of the Voronoi cell (= number of neighbors)
* **Column 3 - Neighbor distribution:** Count of neighbors with each number of sides (3-sided, 4-sided, 5-sided, etc.)
* **Column 4 - p-vector:** Canonical sequence describing the topology: (central_neighbors, neighbor1_count, neighbor2_count, ...)
* **Column 5 - Symmetry order:** Order of the symmetry group of the local arrangement
* **Column 6 - Chirality:** Handedness of the arrangement (-1=left-handed, 0=non-chiral, +1=right-handed)

**Interpreting 2D results:**
- **Particle 1:** Has 7 neighbors arranged as (7,4,8,7,5,5,5,10), meaning one neighbor has 4 edges, the next has 8 edges, the next has 7 edges, the next three each have 5 edges, and the last neighbor has 10.
- **Particle 3:** Shows chirality (-1), indicating its neighborhood arrangement has a definite handedness
- **Particle 5:** Has symmetry order 2 and is non-chiral (0), indicating an arrangement with mirror symmetry.

----

### Three-dimensional output (`IdealGas3D.vectors`)

The file contains more complex information reflecting 3D polyhedra:

```text
688    6     (2,2,2)       (1,2,3,1,3,4,5,1,5,6,7,2,7,8,4,8,6,8,7,6,5,4,3,2,1)                   4      0
135    7     (2,3,0,2)     (1,2,3,1,3,4,5,1,5,6,7,8,2,8,9,10,4,10,6,10,9,7,9,8,7,6,5,4,3,2,1)    4      0       
133    8     (1,3,3,1)     (1,2,3,1,3,4,5,6,1,6,7,8,2,8,9,10,4,10,11,5,11,12,7,12,9,12,11,10,9,8,7,6,5,4,3,2,1)  2  0       
422    8     (0,4,4)       (1,2,3,4,1,4,5,6,1,6,7,8,2,8,9,10,3,10,11,5,11,12,7,12,9,12,11,10,9,8,7,6,5,4,3,2,1)  8  0
966    8     (2,3,1,1,1)   (1,2,3,1,3,4,5,1,5,6,7,2,7,8,9,10,4,10,11,12,6,12,8,12,11,9,11,10,9,8,7,6,5,4,3,2,1)  1 -1
```

**Column descriptions:**
* **Column 1 - Particle ID:** Unique identifier from the input file
* **Column 2 - Face count:** Number of faces of the Voronoi polyhedron (= number of neighbors)
* **Column 3 - Face distribution:** Count of faces with each number of sides (3-sided, 4-sided, 5-sided, etc.)
* **Column 4 - Weinberg vector:** Complete topological description of the 3D Voronoi cell
* **Column 5 - Symmetry order:** Order of the symmetry group of the Voronoi polyhedron
* **Column 6 - Chirality:** Handedness of the polyhedron (-1=left-handed, 0=non-chiral, +1=right-handed)

**Interpreting 3D results:**
- **Particle 688:** Has 6 faces: two triangular (3-sided), two quadrilateral (4-sided), and two pentagonal (5-sided)
- **Particle 966:** Shows chirality (-1), indicating the polyhedron cannot be superimposed on its mirror image
- **Symmetry orders** range from 1 (no symmetry) to 8 (high symmetry), reflecting the regularity of the local arrangement

----

## Understanding the data patterns

**Distribution characteristics:**
In ideal gas systems, you'll observe:
- **Diverse topologies:** Many different p-vectors and Weinberg vectors, reflecting the random nature of the arrangements
- **Average coordination:** The average number of neighbors is 6 in 2D and varies in 3D, consistent with mathematical expectations
- **Chirality distribution:** Roughly equal numbers of left-handed and right-handed arrangements, with some non-chiral configurations
- **Low symmetry:** Most arrangements have low symmetry orders due to the random positioning

**Contrast with ordered systems:**
In crystalline systems (which you'll explore in later examples), you would see:
- **Uniform topologies:** All or most particles have identical Voronoi topologies
- **High symmetry:** Higher symmetry orders reflecting the regular crystal structure
- **Predictable patterns:** Specific, reproducible p-vectors or Weinberg vectors characteristic of the crystal structure

**Physical significance:**
- **Structural fingerprints:** Each topology serves as a "fingerprint" of the local structural environment
- **Defect identification:** Deviations from expected topologies can indicate structural defects
- **Phase characterization:** Different phases (crystal, liquid, glass) have characteristic topology distributions
- **Robustness:** The topological description remains stable under thermal vibrations and small perturbations
<br><br>



## Next steps and connections:

This basic topology calculation forms the foundation for all subsequent *VoroTop* analyses:

- **Example 2** will show how to use filters to classify these topologies into structural categories
- **Example 3** will demonstrate how to analyze topology distributions to characterize different phases
- **Example 4** will explore spatial correlations between particles with different topologies
- **Examples 5-6** will show how to visualize and cluster particles based on their topologies

The rich topological information computed here enables automated, quantitative analysis of structure in particle systems—from identifying individual defects to characterizing entire phases and their transitions.
<br><br>

## References

* Lazar, E.A. "*VoroTop: Voronoi Cell Topology Visualization and Analysis Toolkit*", [arXiv](https://arxiv.org/abs/1804.04221), [Model. Simul. Mater. Sci. Eng. 26:1](https://iopscience.iop.org/article/10.1088/1361-651X/aa9a01), 2017.

* Lazar, E.A., Lu, J., Rycroft, C.H., and Schwarcz, D., "*Characterizing structural features of two-dimensional particle systems through Voronoi topology*", [arXiv](https://arxiv.org/abs/2406.00553), [Model. Simul. Mater. Sci. Eng. 32:085022](https://iopscience.iop.org/article/10.1088/1361-651X/ad8ad9), 2024.

## Acknowledgments
This README was written with assistance from Claude (Anthropic).
