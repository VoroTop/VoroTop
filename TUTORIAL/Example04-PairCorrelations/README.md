# Example 4: Compute the Voronoi pair correlation function


**Learning objectives:**

* After completing this example, you will:
  * Be able to automate the computation of the Voronoi pair correlation function
  * Understand how to interpret the results of this analysis 
<br><br>

**Background:**

Roughly speaking, the classical [pair](https://faculty.college.emory.edu/sites/weeks/idl/gofr.html) [correlation](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Advanced_Statistical_Mechanics_(Tuckerman)/07%3A_Distribution_Functions_and_Liquid_Structure/7.06%3A_The_Pair_Correlation_Function) [function](https://homepage.univie.ac.at/franz.vesely/simsp/dx/node22.html) measures the probability of finding neighbors at various distances from a central particle.  Voronoi topology can be used to define a discrete version of this classical tool.  This discretized approach has advantages in studying many ordered and disordered systems; further details can be found in 
* Worlitzer, V.M., Ariel, G., and Lazar, E.A., "*Pair correlation function based on Voronoi topology*", [arXiv](https://arxiv.org/abs/2210.09731), [Phys. Rev. E 108:064115](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.108.064115), 2023. 

*VoroTop* automates the calculation of the Voronoi pair correlation function, in both its normalized and unnormalized forms.
<br><br>


## The unnormalized Voronoi pair correlation function -u

**Input data:**

* **Data files:** We include three different 'disordered' systems: `IdealGas2D`, `PerturbedSquareLattice2D`, and `Polycrystal2D` used previously in Example 3.  Each of these files are in the LAMMPS dump format; initial lines describe system dimensions and boundary conditions.  Subsequent lines include particle ID, x, y, and z coordinates.
<br>


**Steps:**

1.  **Navigate to the example directory:** Open your terminal or command prompt and navigate to the `Example04-PairCorrelations/` directory.

    ```bash
    cd VoroTop/Tutorial/Example04-PairCorrelations
    ```

2.  **Execute the *VoroTop* command:** Run the following *VoroTop* command (or execute the provided script) to perform the analysis:

    ```bash
    VoroTop IdealGas2D -2 -u 
    VoroTop Polycrystal2D -2 -u 10
    VoroTop PerturbedSquareLattice2D -2 -u 15
    ```
* The first argument of the command always specifies the data file.
* The -2 switch specifies that the system and its Voronoi topology analysis should be performed as two-dimensional.
* The -u switch indicates that *VoroTop* should compute the unnormalized Voronoi pair correlation function.  More precisely, for each natural number *k*, the function *u(k)* is the average number of neighbors at a distance *k* from a central particle.
* By default, *VoroTop* computes the Voronoi pair correlation function up to a Voronoi distance of 20.  However, the user can provide an integer after the -v switch to change this distance.

3.  **Examine the output:** Output of this analysis is printed to the screen.
<br><br>

**Expected output and explanation:**

The Voronoi pair correlation function is printed to the screen.  By default, the maximum distance considered is 20, though this can be changed by including 
an integer after the -v command.

Here is the output for each of the examples.

#### `IdealGas2D`
   ```text
1    6         1.78958   1.33775   0.000405093
2    13.7155   8.29476   2.88006   0.00119881
3    23.0147   15.8079   3.97591   0.00253068
4    33.0396   24.0349   4.90254   0.0044427
5    43.5542   31.6067   5.62198   0.00696319
6    54.4166   38.9857   6.24386   0.0101123
7    65.5527   45.7799   6.76608   0.0139059
8    76.8186   51.8274   7.19913   0.0183514
9    88.1383   61.0108   7.81094   0.023452
10   99.5878   70.5753   8.40091   0.0292152
11   111.267   79.1882   8.89877   0.0356542
12   122.916   86.4384   9.29722   0.0427674
13   134.537   94.6384   9.72822   0.0505531
14   146.356   102.147   10.1068   0.0590228
15   158.103   108.141   10.3991   0.0681723
16   170.071   111.257   10.5479   0.0780143
17   182.096   115.773   10.7598   0.0885523
18   194.123   124.43    11.1548   0.0997863
19   206.103   131.746   11.4781   0.111714
20   218.054   136.685   11.6912   0.124332
   ```

#### `Polycrystal2D`
   ```text
1    6         0.0234954   0.153282   0.000405093
2    12.0235   0.0687767   0.262253   0.0011009
3    18.0677   0.11405     0.337713   0.00214648
4    24.1259   0.19074     0.436738   0.00354266
5    30.2028   0.278904    0.528114   0.0052905
6    36.2984   0.38215     0.618183   0.0073911
7    42.409    0.534548    0.731128   0.00984533
8    48.5384   0.689959    0.830637   0.0126543
9    54.6845   0.858093    0.926333   0.0158189
10   60.8471   1.06748     1.03319    0.0193401
   ```

#### `PerturbedSquareLattice2D`
   ```text
1    6         1.79947   1.34144   0.000401745
2    13.705    8.61523   2.93517   0.0011883
3    22.976    19.0124   4.36032   0.00250695
4    33.0972   38.9632   6.24205   0.00440646
5    43.7693   78.826    8.8784    0.00691848
6    54.9392   151.45    12.3065   0.0100715
7    66.6108   265.459   16.2929   0.0138945
8    78.6103   428.582   20.7022   0.0184061
9    90.9167   634.009   25.1795   0.023624
10   103.483   880.291   29.6697   0.0295631
11   116.281   1167.25   34.165    0.0362367
12   129.149   1470.31   38.3446   0.0436489
13   142.062   1784.48   42.2432   0.0518021
14   155.056   2115.57   45.9953   0.0607011
15   168.055   2442.03   49.4169   0.0703461
   ```

* **Column 1** lists a discrete Voronoi distance *k*, measuring the number of Voronoi cells that needs to traversed from one particle to another.  
* **Column 2** lists the average number of *k*-neighbors of particles in the system.
* **Column 3** lists the sample variance of numbers of *k*-neighbors.
* **Column 4** lists the sample standard deviation of numbers of *k*-neighbors.
* **Column 5** lists the fraction of particles in the system that are at most a Voronoi distance of *k* from a central particle.  As this fraction approaches 1, such as in small systems, the reported data is less reliable.
<br>

Notice that *u(1)=6* in all systems. This reflects the fact that when particles are in *general position*, Voronoi cells have on average six edges (assuming periodic boundary conditions).


## The (unnormalized) Voronoi pair correlation function -v

**Input data:**

* **Data files:** We include three different systems, each containing 108,000 atoms.
   * **HighT-Cu** A single crystal of coppper (Cu), heated to 85% of its bulk melting temperature.
   * **Liquid-Cu** The above crystal was superheated to 120% of its bulk melting temperature and allowed to melt.
   * **IdealGas3D** This was created by randomly distributing particles in space, a model of an ideal gas.

Each file is in the LAMMPS dump format; initial lines describe system dimensions and boundary conditions.  Subsequent lines include particle ID, x, y, and z coordinates.
<br>


**Steps:**

1.  **Navigate to the example directory:** Open your terminal or command prompt and navigate to the `Example04-PairCorrelations/` directory.

    ```bash
    cd VoroTop/Tutorial/Example04-PairCorrelations
    ```

2.  **Execute the *VoroTop* command:** Run the following *VoroTop* command (or execute the provided script) to perform the analysis:

    ```bash
    VoroTop HighT-Cu -v
    VoroTop Liquid-Cu -v 10
    VoroTop IdealGas3D -v 15
    ```
* The first argument of the command always specifies the data file.
* The -v switch indicates that *VoroTop* should compute the normalized Voronoi pair correlation function.
* By default, *VoroTop* computes the unnormalized Voronoi pair correlation function up to a Voronoi distance of 20.  However, the user can provide an integer after the -u switch to change this distance.

3.  **Examine the output:** Output of this analysis is printed to the screen.
<br><br>

**Expected output and explanation:**

The Voronoi pair correlation function is printed to the screen.  By default, the maximum distance considered is 20, though this can be changed by including an integer after the -u command.

Here is the output for each of the examples.

#### `HighT-Cu`
   ```text
1    0.90766    0.00344719    0.0587127    0.000139824
2    0.768431   0.00140585    0.0374947    0.000637544
3    0.674016   0.000471621   0.0217168    0.00178535
4    0.61372    0.000190122   0.0137885    0.00386689
5    0.573068   9.68224e-05   0.00983984   0.00716678
6    0.54388    5.60522e-05   0.0074868    0.0119691
7    0.521936   3.57028e-05   0.00597518   0.0185582
8    0.504838   2.42606e-05   0.0049255    0.0272189
9    0.491124   1.79054e-05   0.00423148   0.0382359
10   0.479887   1.39825e-05   0.00373932   0.051895
11   0.470501   1.12704e-05   0.00335714   0.0684821
12   0.462551   9.32594e-06   0.00305384   0.0882841
13   0.45572    7.93464e-06   0.00281685   0.111588
14   0.449785   6.91703e-06   0.00263003   0.138681
15   0.444596   6.14393e-06   0.0024787    0.169852
16   0.44001    5.54393e-06   0.00235456   0.20539
17   0.435933   5.08204e-06   0.00225434   0.245585
18   0.432287   4.71137e-06   0.00217057   0.290726
19   0.428881   4.2715e-06    0.00206676   0.34109
20   0.414947   1.30096e-05   0.00360688   0.395541
   ```

#### `Liquid-Cu`
   ```text
1    0.935769   0.00862073    0.0928479   0.000143867
2    0.817385   0.00333674    0.0577645   0.000673295
3    0.744292   0.00168795    0.0410847   0.00194078
4    0.707263   0.000993      0.0315119   0.00433959
5    0.684275   0.000628425   0.0250684   0.00827983
6    0.667666   0.00042497    0.0206148   0.0141752
7    0.655054   0.000306642   0.0175112   0.0224448
8    0.645217   0.000230284   0.0151751   0.0335137
9    0.637356   0.000180822   0.013447    0.047811
10   0.630891   0.000145285   0.0120534   0.0657682
   ```

#### `IdealGas3D`
   ```text
1    1.00021    0.0463382     0.215263    0.000153136
2    1.0007     0.0331548     0.182084    0.000801297
3    1.00092    0.0181241     0.134626    0.0025058
4    1.00101    0.0109        0.104403    0.00590092
5    1.00098    0.00713556    0.0844722   0.0116648
6    1.00073    0.00501364    0.0708071   0.0205011
7    1.00059    0.00369884    0.0608181   0.0331329
8    1.00047    0.00281837    0.0530884   0.0502961
9    1.00038    0.00220681    0.0469767   0.072737
10   1.0003     0.00177857    0.0421731   0.101209
11   1.00017    0.001468      0.0383145   0.136469
12   1.00012    0.00123532    0.0351471   0.179284
13   0.999987   0.00106263    0.032598    0.23042
14   0.999921   0.000925713   0.0304255   0.290651
15   0.999864   0.00081281    0.0285098   0.360753
   ```

* **Column 1** lists a discrete Voronoi distance *k*, measuring the number of Voronoi cells that needs to traversed from one particle to another.  
* **Column 2** lists the average number of *k*-neighbors of particles in the system.
* **Column 3** lists the sample variance of numbers of *k*-neighbors.
* **Column 4** lists the sample standard deviation of numbers of *k*-neighbors.
* **Column 5** lists the fraction of particles in the system that are at most a Voronoi distance of *k* from a central particle.  As this fraction approaches 1, such as in small systems, the reported data is less reliable.

<br>

Further details about this method can be found in the references below, in particular the first.
<br>

## References

* Worlitzer, V.M., Ariel, G., and Lazar, E.A., "*Pair correlation function based on Voronoi topology*", [arXiv](https://arxiv.org/abs/2210.09731), [Phys. Rev. E 108:064115](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.108.064115), 2023. 

* Lazar, E.A., Lu, J., Rycroft, C.H., and Schwarcz, D., "*Characterizing structural features of two-dimensional particle systems through Voronoi topology*", [arXiv](https://arxiv.org/abs/2406.00553), [Model. Simul. Mater. Sci. Eng. 32:085022](https://iopscience.iop.org/article/10.1088/1361-651X/ad8ad9), 2024. 

* Lazar, E.A. "*VoroTop: Voronoi Cell Topology Visualization and Analysis Toolkit*", [arXiv](https://arxiv.org/abs/1804.04221), [Model. Simul. Mater. Sci. Eng. 26:1](https://iopscience.iop.org/article/10.1088/1361-651X/aa9a01), 2017.
