# Example 2: Voronoi topology analysis through filter files


**Learning objectives:**

* After completing this example, you will:
  * Understand how a filter file is structured
  * Know how to use a filter file to identify structural defects in crystalline systems.
  * Understand the LAMMPS dump file created by *VoroTop*, which contains structure classification information.
<br><br>

**Background:**

A filter is a list of one or more families of Voronoi topologies used by *VoroTop* to identify crystalline and defect structure.  As a simple example, a filter can enumerate only the unique $p$-vector $(6,6,6,6,6,6,6)$ associated with the ideal hexagonal lattice in two dimensions, or a the family of Voronoi topologies associated with FCC crystals in three dimensions.
<br>

Filter files are divided into three parts: 
1. The first part consists of optional comments about the filter, such as its source and number of Voronoi topologies; all lines that begin with a ‘#’ are treated as comments. 
2. Lines in the second part begin with an astrisk * and specify user-defined structure types. Each such line, after the *, includes an index and a name for the structure type. Indices of structure types are listed in increasing order and begin with 1.
3. Remaining lines record Voronoi cell topologies.  Each line begins with a structure type index, and is followed by either a *p*-vector in two dimensions or else a Weinberg vector in three dimensions, specifying the Voronoi topology. 
<br><br>




We use a filter for hexagonal crystals to identify defects in a two-dimensional polycrystal.


**Input data:**

* **Data file:** `Polycrystal2D`.  LAMMPS dump format; initial lines describe system dimensions and boundary conditions.  Subsequent lines include particle ID, x, y, and z coordinates.
* **Filter file:** `crystal_2d.filter`.  Specifies a single structure type called Crystal and associated with it a single two-dimensional Voronoi topology. 
<br>

   ```text
   #       Filter with unique type in two-dimensional defect-free
   #       triagonal crystal, in which all particles are hexagons 
   #       surrounded by hexagons.
   *       1       Crystal         
   1       (6,6,6,6,6,6,6)
   ```
<br>

**Steps:**

1.  **Navigate to the example directory:** Open your terminal or command prompt and navigate to the `Example02-FilterAnalysis/` directory.

    ```bash
    cd VoroTop/Tutorial/Example02-FilterAnalysis
    ```

2.  **Execute the *VoroTop* command:** Run the following *VoroTop* command (or execute the provided script) to perform the analysis:

    ```bash
    VoroTop Polycrystal2D -l -2 -f crystal_2d.filter
    ```
* The first argument of the command always specifies the data file; in this case it is `Polycrystal2D`.
* The -l switch specifies that *VoroTop* should output a LAMMPS dump file after performing Voronoi topology analysis.
* The -2 switch specifies that the system and its Voronoi topology analysis should be performed as two-dimensional.
* The -f switch specifies a filter file that should be applied to the data.

3.  **Examine the output file:** A file named `Polycrystal2D.dump` will be created in this directory. It can be opened with a text editor to view the results, or else in OVITO, in which results of the Voronoi topology analysis can be visualized.
<br><br>

**Expected output and explanation:**

The file `Polycrystal2D.dump` is a LAMMPS dump file, containing system information as well as one line per atom, including information about coordinates and structure type. Here are several example lines from the file.

```text
ITEM: ATOMS id x y z vt
16210   64.538    6.69981   0   1       
1898    66.8254   5.30426   0   0       
16205   67.4699   6.20108   0   0       
16530   65.0532   22.4752   0   1       
13199   45.799    6.69556   0   1       
```

* **Column 1** lists the particle ID, as specified in the input file.  
* **Columns 2-4** lists the x, y, and z coordinates of the particles; notice that all z-coordinates are zero, a result of the system being two-dimensional.
* **Column 5** lists the structure types of the particles, classified through Voronoi topology (vt) using the provided filter.  In this example, particles with structure type 1 are Crystal; all other particles are assigned structure type 0.  If a Voronoi topology is associated with multiple structure types, then only the lowest-index one is provided.  Further analysis is possible using an option to resolve indeterminate types; this feature is not currently implemented.
<br><br>



## References

* Lazar, E.A., Lu, J., Rycroft, C.H., and Schwarcz, D., "*Characterizing structural features of two-dimensional particle systems through Voronoi topology*", [arXiv](https://arxiv.org/abs/2406.00553), [Model. Simul. Mater. Sci. Eng. 32:085022](https://iopscience.iop.org/article/10.1088/1361-651X/ad8ad9), 2024. 

* Lazar, E.A. "*VoroTop: Voronoi Cell Topology Visualization and Analysis Toolkit*", [arXiv](https://arxiv.org/abs/1804.04221), [Model. Simul. Mater. Sci. Eng. 26:1](https://iopscience.iop.org/article/10.1088/1361-651X/aa9a01), 2017.
