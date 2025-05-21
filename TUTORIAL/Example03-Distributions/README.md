# Example 3: Calculate the distribution of Voronoi topologies


**Learning objectives:**

* After completing this example, you will:
  * Be able to automate the computation of the distribution of Voronoi topologies in a system
  * Understand the structure of a distribution file
  * Use the results of this analysis to analyze disordered systems
<br><br>

**Background:**

The structure of particle systems can be studied by considering the distribution of Voronoi topologies in the system.  In many crystals, for example, the Voronoi cells of all particles are identical and thus have the same topology.  In many systems, however, particles are distributed in a more `random' manner.  However, certain arrangements of particles might still be more common than others for both entropic and energetic reasons.
<br><br>


## Two dimensions

**Input data:**

* **Data files:** We include three different 'disordered' systems: `IdealGas2D`, `PerturbedSquareLattice2D`, and `Polycrystal2D`.  Each of these files are in the LAMMPS dump format; initial lines describe system dimensions and boundary conditions.  Subsequent lines include particle ID, x, y, and z coordinates.
<br>


**Steps:**

1.  **Navigate to the example directory:** Open your terminal or command prompt and navigate to the `Example03-Distributions/` directory.

    ```bash
    cd VoroTop/Tutorial/Example03-Distributions
    ```

2.  **Execute the *VoroTop* command:** Run the following *VoroTop* command (or execute the provided script) to perform the analysis:

    ```bash
    VoroTop IdealGas2D -2 -d
    VoroTop Polycrystal2D -2 -d
    VoroTop PerturbedSquareLattice2D -2 -d
    ```
* The first argument of the command always specifies the data file.
* The -2 switch specifies that the system and its Voronoi topology analysis should be performed as two-dimensional.
* The -d switch indicates that *VoroTop* should compute the distribution of Voronoi topologies in a system.

3.  **Examine the output file:** A file named `<datafile>.distribution` will be created in this directory. It can be opened with a text editor to view the results and be further analyzed.
<br><br>

**Expected output and explanation:**

Each distribution file `<datafile>.distribution` begins with an optional section of comments, which can record information about the source of the distribution.  Each subsequent line records one Voronoi topology that appears at least once in the system, as well as additional information.

Here are several lines from each of the examples.

#### `IdealGas2D.distribution`
   ```text
   #       DISTRIBUTION CREATED FROM IdealGas2D
   #       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
   (4,6,6,7,8)     66      41      0       25
   (4,5,6,7,7)     61      33      0       28
   (4,5,7,6,8)     58      27      0       31
   (4,6,7,7,8)     58      28      0       30
   (5,5,6,7,6,7)   56      28      0       28
   ```

#### `PerturbedSquareLattice2D.distribution`
   ```text
   #       DISTRIBUTION CREATED FROM PerturbedSquareLattice2D
   #       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
   (5,5,6,7,6,7)   72      40      0       32
   (4,6,6,7,7)     62      0       62      0
   (4,6,6,7,8)     62      34      0       28
   (4,6,6,6,7)     54      0       54      0
   (5,6,6,7,6,8)   52      27      0       25
   ```

#### `Polycrystal2D.distribution`
   ```text
   #       DISTRIBUTION CREATED FROM Polycrystal2D
   #       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
   (6,6,6,6,6,6,6) 15707   0       15707   0
   (6,6,6,6,6,6,7) 427     0       427     0
   (6,5,6,6,6,6,7) 336     173     0       163
   (6,5,6,6,6,6,6) 205     0       205     0
   (5,6,6,6,6,7)   138     0       138     0
   ```


* **Column 1** lists a Voronoi topology.  In two dimensions, this is denoted by a *p*-vector; in three dimensions it is denoted by a Weinberg vector.  
* **Column 2** lists the total count of particles with the given Voronoi topology; this is the sum of Columns 3-5.
* **Column 3** lists the number of particles with a left-handed version of this Voronoi topology; this can only be non-zero in chiral arrangements.
* **Column 4** lists the number of particles with no chirality.  Particles whose Voronoi cell topology has no chirality 
* **Column 5** lists the number of particles with a right-handed version of this Voronoi topology; this can only be non-zero in chiral arrangements.

<br>



## Three dimensions

**Input data:**

* **Data files:** We include three different systems, each containing 108,000 atoms.
   * **HighT-Cu** A single crystal of coppper (Cu), heated to 85% of its bulk melting temperature.
   * **Liquid-Cu** The above crystal was superheated to 120% of its bulk melting temperature and allowed to melt.
   * **IdealGas3D** This was created by randomly distributing particles in space, a model of an ideal gas.

Each file is in the LAMMPS dump format; initial lines describe system dimensions and boundary conditions.  Subsequent lines include particle ID, x, y, and z coordinates.
<br>


**Steps:**

1.  **Navigate to the example directory:** Open your terminal or command prompt and navigate to the `Example03-Distributions/` directory.

    ```bash
    cd VoroTop/Tutorial/Example03-Distributions
    ```

2.  **Execute the *VoroTop* command:** Run the following *VoroTop* command (or execute the provided script) to perform the analysis:

    ```bash
    VoroTop HighT-Cu -d
    VoroTop Liquid-Cu -d
    VoroTop IdealGas3D -d
    ```
* The first argument of the command always specifies the data file.
* The -d switch indicates that *VoroTop* should compute the distribution of Voronoi topologies in a system.

3.  **Examine the output file:** A file named `<datafile>.distribution` will be created in this directory. It can be opened with a text editor to view the results and be further analyzed.
<br><br>

**Expected output and explanation:**

Each distribution file `<datafile>.distribution` begins with an optional section of comments, which can record information about the source of the distribution.  Each subsequent line records one Voronoi topology that appears at least once in the system, as well as additional information.

Here are several lines from each of the examples.

#### `HighT-Cu.distribution`
   ```text
#       DISTRIBUTION CREATED FROM dump.melt.010000
#       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
(1,2,3,4,1,4,5,6,7,1,7,8,9,2,9,10,11,12,3,12,13,14,5,14,15,16,17,6,17,18,19,8,19,20,21,10,21,22,11,22,23,24,13,24,15,24,23,25,16,25,26,18,26,
          20,26,25,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)      9868    4946    0       4922
(1,2,3,4,1,4,5,6,7,1,7,8,9,2,9,10,11,12,3,12,13,14,5,14,15,16,6,16,17,18,8,18,19,20,10,20,21,11,21,22,13,22,23,15,23,24,17,24,19,24,23,22,21,
          20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)                        7272    3662    0       3610
(1,2,3,4,1,4,5,6,7,1,7,8,9,10,2,10,11,12,3,12,13,14,5,14,15,16,6,16,17,18,8,18,19,9,19,20,21,11,21,22,13,22,23,24,15,24,17,24,23,20,23,22,21,
          20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)                        6329    0       6329    0
(1,2,3,4,1,4,5,6,7,1,7,8,9,10,2,10,11,12,13,3,13,14,15,5,15,16,17,6,17,18,19,8,19,20,9,20,21,11,21,22,23,12,23,14,23,22,24,16,24,18,24,22,21,
          20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)                        6114    3108    0       3006
(1,2,3,4,1,4,5,6,7,1,7,8,9,2,9,10,11,12,3,12,13,14,5,14,15,16,6,16,17,18,8,18,19,10,19,20,21,11,21,13,21,20,22,15,22,17,22,20,19,18,17,16,15,
          14,13,12,11,10,9,8,7,6,5,4,3,2,1)                                          4961    2498    0       2463
   ```

#### `Liquid-Cu.distribution`
   ```text
#       DISTRIBUTION CREATED FROM dump.melt.040000
#       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
(1,2,3,4,1,4,5,6,7,1,7,8,9,10,2,10,11,12,3,12,13,14,5,14,15,16,6,16,17,8,17,18,19,9,19,20,11,20,21,13,21,22,15,22,18,22,21,20,19,18,17,16,15,
          14,13,12,11,10,9,8,7,6,5,4,3,2,1)                       1457    0       1457    0
(1,2,3,4,1,4,5,6,7,1,7,8,9,2,9,10,11,12,3,12,13,14,5,14,15,16,6,16,17,18,8,18,19,10,19,20,21,11,21,13,21,20,22,15,22,17,22,20,19,18,17,16,15,
          14,13,12,11,10,9,8,7,6,5,4,3,2,1)                       1349    655     0       694
(1,2,3,4,1,4,5,6,7,1,7,8,9,2,9,10,11,12,3,12,13,14,5,14,15,16,6,16,17,8,17,18,10,18,19,20,11,20,13,20,19,15,19,18,17,16,15,14,13,12,11,10,9,
          8,7,6,5,4,3,2,1)                                         924    0       924     0
(1,2,3,4,1,4,5,6,7,1,7,8,9,10,2,10,11,12,13,3,13,14,15,5,15,16,17,6,17,18,8,18,19,20,9,20,21,11,21,22,23,12,23,14,23,22,24,16,24,19,24,22,21,
          20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)      779    380     0       399
(1,2,3,4,1,4,5,6,7,1,7,8,9,2,9,10,11,12,3,12,13,14,5,14,15,16,6,16,17,18,8,18,19,20,10,20,21,11,21,22,13,22,23,15,23,24,17,24,19,24,23,22,21,
          20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)      775    388     0       387
   ```

#### `IdealGas3D.distribution`
   ```text
   #       DISTRIBUTION CREATED FROM PV-108000-0-0.data
   #       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
   (1,2,3,1,3,4,5,6,1,6,7,8,9,2,9,10,11,4,11,12,5,12,13,7,13,14,8,14,10,14,13,12,11,10,9,8,7,6,5,4,3,2,1)                     320     162     0    158
   (1,2,3,1,3,4,5,6,1,6,7,8,2,8,9,10,4,10,11,5,11,12,7,12,9,12,11,10,9,8,7,6,5,4,3,2,1)                                       172       0    172     0
   (1,2,3,4,1,4,5,6,1,6,7,8,9,2,9,10,11,3,11,12,13,5,13,14,7,14,15,8,15,16,10,16,12,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)   156      79      0    77
   (1,2,3,4,1,4,5,6,1,6,7,8,2,8,9,10,3,10,11,12,5,12,13,7,13,14,9,14,11,14,13,12,11,10,9,8,7,6,5,4,3,2,1)                     131       0    131     0
   (1,2,3,4,1,4,5,6,1,6,7,8,2,8,9,10,3,10,11,12,5,12,13,14,7,14,15,9,15,16,11,16,13,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)   121      71      0    50
   ```

* **Column 1** lists a Voronoi topology.  In three dimensions, this is denoted by a Weinberg vector.  
* **Column 2** lists the total count of particles with the given Voronoi topology; this is the sum of Columns 3-5.
* **Column 3** lists the number of particles with a left-handed version of this Voronoi topology; this can only be non-zero in chiral arrangements.
* **Column 4** lists the number of particles with no chirality.  Particles whose Voronoi cell topology has no chirality 
* **Column 5** lists the number of particles with a right-handed version of this Voronoi topology; this can only be non-zero in chiral arrangements.

<br><br>



## References

* Lazar, E.A., Lu, J., Rycroft, C.H., and Schwarcz, D., "*Characterizing structural features of two-dimensional particle systems through Voronoi topology*", [arXiv](https://arxiv.org/abs/2406.00553), [Model. Simul. Mater. Sci. Eng. 32:085022](https://iopscience.iop.org/article/10.1088/1361-651X/ad8ad9), 2024. 

* Lazar, E.A. "*VoroTop: Voronoi Cell Topology Visualization and Analysis Toolkit*", [arXiv](https://arxiv.org/abs/1804.04221), [Model. Simul. Mater. Sci. Eng. 26:1](https://iopscience.iop.org/article/10.1088/1361-651X/aa9a01), 2017.
