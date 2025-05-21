# Example 1: Basic Voronoi topology calculations

**Learning objectives:**

* After completing this example, you will be able to:
    * Run *VoroTop* to calculate the Voronoi topologies of individual particles.
    * Interpret the computed Voronoi topology for individual particles in two and three dimensions.
    * Understand additional recorded data about symmetry and chirality.
<br><br>

**Background:**

The Voronoi cell of a particle is the region of space closer to the given particle than to any other particle.   In two dimensions, this region is a polygon, and in three dimensions it is a polyhedron.   The *topology* of the Voronoi cell is a combinatorial desription of how its neighbors are arranged. In two dimensions, the Voronoi cell topology of a particle is a list of natural numbers, specifying the number  of neighbors of the central particle and the numbers of neighbors of each of its neighbors, ordered in a canonical way
<br><br>

**Input data:**

* **File:** IdealGas2D.dump, IdealGas3D.dump
* **Format:** LAMMPS dump format; initial lines describing systems dimensions and boundary conditions.  Subsequent lines include particle information: particle ID, x, y, and z coordinates.
* **Description:** Two systems of 1000 particles, uniformly chosen in the unit square (2D) or cube (3D).  This random configuration models an ideal gas.
<br><br>

**Steps:**

1.  **Navigate to the example directory:** Open your terminal or command prompt and navigate to the `Example01-Basic/` directory within the *VoroTop* tutorial repository.

    ```bash
    cd VoroTop/Tutorial/Example01-Basic
    ```

2.  **Execute the *VoroTop* command:** Run the following *VoroTop* command (or execute the provided script) to perform the analysis:

    ```bash
    VoroTop IdealGas3D.dump -vt
    VoroTop IdealGas2D.dump -2 -vt
    ```


3.  **Examine the output file:** Files named `IdealGas2D.vectors` and `IdealGas3D.vectors` will be created in this directory. Open them with a text editor to view the results.
<br><br>

**Expected output and explanation:**

The file `IdealGas2D.vectors` contains one line for each particle, such as the following:

```text
1     7     (0,1,3,0,1,1,0,1)    (7,4,8,7,5,5,5,10)      1      1
2     5     (0,0,3,0,1,0,1)      (5,5,5,7,5,9)           1      1
3     8     (0,1,2,2,2,1)        (8,4,6,7,5,6,5,8,7)     1     -1
4     5     (0,0,1,2,2)          (5,5,6,7,6,7)           1      1
5     5     (0,0,0,2,1,2)        (5,6,7,6,8,8)           2      0
```

* **Column 1** lists the particle ID, as specified in the input file.  
* **Column 2** lists the number of edges of the Voronoi cell of the given particle; this can be understood as a count of neighbors.
* **Column 3** lists the number of Voronoi neighbors with a given number of sides.  For example, the first particle has no neighbors with three-sided Voronoi cells, one neighbor with a four-sided Voronoi cell, three neighbors with five-sided Voronoi cells, and so forth.
* **Column 4** is the canonical *p*-vector, describing the arrangement of its neighbors.
* **Column 5** is order of the symmetry group of the arrangement of the Voronoi cells of the particle and its neighbors.
* **Column 6** is the chirality of the arrangement of the Voronoi cells of the particle and its neighbors; this number can be -1, 0, or 1.
<br><br>

The file `IdealGas3D.vectors` contains one line for each particle, such as the following:

```text
688    6     (2,2,2)       (1,2,3,1,3,4,5,1,5,6,7,2,7,8,4,8,6,8,7,6,5,4,3,2,1)                                     4      0
135    7     (2,3,0,2)     (1,2,3,1,3,4,5,1,5,6,7,8,2,8,9,10,4,10,6,10,9,7,9,8,7,6,5,4,3,2,1)                      4      0       
133    8     (1,3,3,1)     (1,2,3,1,3,4,5,6,1,6,7,8,2,8,9,10,4,10,11,5,11,12,7,12,9,12,11,10,9,8,7,6,5,4,3,2,1)    2      0       
422    8     (0,4,4)       (1,2,3,4,1,4,5,6,1,6,7,8,2,8,9,10,3,10,11,5,11,12,7,12,9,12,11,10,9,8,7,6,5,4,3,2,1)    8      0
966    8     (2,3,1,1,1)   (1,2,3,1,3,4,5,1,5,6,7,2,7,8,9,10,4,10,11,12,6,12,8,12,11,9,11,10,9,8,7,6,5,4,3,2,1)    1     -1
```

* **Column 1** lists the particle ID, as specified in the input file.  
* **Column 2** lists the number of faces of the Voronoi cell of the given particle; this can be understood as a count of neighbors.
* **Column 3** lists the number of faces with a given number of sides, beginning with three.  For example, particle 688 has two three-sided faces, two four-sided faces, and two five-sided faces.
* **Column 4** is the canonical Weinberg vector, describing the three-dimensional Voronoi cell topology.
* **Column 5** is order of the symmetry group of the Voronoi cell of the particle.
* **Column 6** is the chirality of the Voronoi cell of the particle; this number can be -1, 0, or 1.



## References

* Lazar, E.A. "*VoroTop: Voronoi Cell Topology Visualization and Analysis Toolkit*", [arXiv](https://arxiv.org/abs/1804.04221), [Model. Simul. Mater. Sci. Eng. 26:1](https://iopscience.iop.org/article/10.1088/1361-651X/aa9a01), 2017.

* Lazar, E.A., Lu, J., Rycroft, C.H., and Schwarcz, D., "*Characterizing structural features of two-dimensional particle systems through Voronoi topology*", [arXiv](https://arxiv.org/abs/2406.00553), [Model. Simul. Mater. Sci. Eng. 32:085022](https://iopscience.iop.org/article/10.1088/1361-651X/ad8ad9), 2024. 

