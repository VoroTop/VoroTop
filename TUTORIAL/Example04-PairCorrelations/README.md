# Example 4: Compute the Voronoi pair correlation function

**Learning objectives:**

* After completing this example, you will:
  * Understand the theoretical foundation of pair correlation functions and their role in structural analysis.
  * Learn how the Voronoi pair correlation function differs from classical radial distribution functions.
  * Compute and interpret both normalized and unnormalized Voronoi pair correlation functions.
  * Analyze spatial correlations and structural ordering in diverse particle systems.
  * Use pair correlation data to distinguish between ordered, disordered, and intermediate structural states.
  * Connect pair correlation results to physical properties and structural evolution processes.
<br><br>


---

## Background:


The pair correlation function is one of the most fundamental tools in statistical mechanics and condensed matter physics for characterizing spatial correlations in particle systems. Roughly speaking, the classical [pair correlation function](https://faculty.college.emory.edu/sites/weeks/idl/gofr.html) quantifies the probability of finding neighbors at various distances from a central particle.  The Voronoi pair correlation function provides a discrete, topology-based alternative that offers unique advantages for analyzing both ordered and disordered systems.  Details can be found in 
* Worlitzer, V.M., Ariel, G., and Lazar, E.A., "*Pair correlation function based on Voronoi topology*", [arXiv](https://arxiv.org/abs/2210.09731), [Phys. Rev. E 108:064115](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.108.064115), 2023. 


Despite similarities, the the classical and Voronoi pair correlation functions differ in important respects:

**Classical pair correlation function g(r):**
- **Continuous:** Depends on Euclidean distance *r*
- **Geometry-dependent:** Sensitive to particle positions and thermal vibrations  
- **Integration-based:** Requires careful choice of bin sizes and integration limits
- **Density-sensitive:** Strongly affected by local density fluctuations
<br><br>

**Voronoi pair correlation function v(k):**
- **Discrete:** Based on a discrete distance *k* - the number of Voronoi cell boundaries to traverse
- **Topology-based:** Robust against thermal vibrations and small coordinate errors 
- **Naturally quantized:** No binning artifacts or resolution issues
- **Universal normalization:** Meaningful comparison across different densities
- **Structure-focused:** Emphasizes connectivity rather than metric distance
<br><br>


---

**Physical interpretation:**

The unnormalized Voronoi pair correlation function *u*(*k*) measures the average number of particles at topological distance *k* from a central particle. Here, "topological distance *k*" means the minimum number of Voronoi cell boundaries that must be crossed to reach from one particle to another through the tessellation network.

- ***k* = 1:** Direct Voronoi neighbors (share a cell boundary)
- ***k* = 2:** Neighbors of neighbors (second shell)
- ***k* = 3:** Third coordination shell, and so on

The classical pair correlation *g*(*r*) is normalized so that for the ideal gas it is equal to 1 for all distances *r*>0.  In a similar fashion, the Voronoi pair correlation is normalized so that it is equal to 1 for all *k* for the ideal gas.

This topological perspective reveals structural ordering that may be obscured in traditional distance-based analyses, particularly in systems with significant thermal motion or structural disorder.
<br><br>


----

**Input data:**

We analyze solids, liquids, and gases to demonstrate the versatility of Voronoi pair correlation analysis:

**Two-dimensional systems:**
* **`Polycrystal2D`** - Lennard-Jones polycrystal with grain boundaries 
* **`Liquid2D`** - Lennard-Jones liquid
* **`IdealGas2D`** - randomly distributed particles

**Three-dimensional systems:**
* **`HighTemperatureCu`** - Copper crystal, heated to 85% of bulk melting temperature
* **`LiquidCu`** - Liquid copper, heated to 120% of bulk melting temperature
* **`IdealGas3D`** - randomly distributed particles

**System characteristics:**
- **Size range:** 10,000-108,000 particles for statistical reliability
- **Structural diversity:** from highly ordered to completely disordered
- **Physical relevance:** Realistic materials systems with practical applications
<br><br>

**Steps:**

1.  **Navigate to the example directory:**

    ```bash
    cd VoroTop/Tutorial/Example04-PairCorrelations
    ```

2.  **Compute unnormalized pair correlation functions for 2D systems:**

    ```bash
    VoroTop Polycrystal2D -2 -u 10
    VoroTop Liquid2D -2 -u 10
    VoroTop IdealGas2D -2 -u 10
    ```

3.  **Compute normalized pair correlation functions for 3D systems:**

    ```bash
    VoroTop HighTemperatureCu -v 10
    VoroTop LiquidCu -v 10
    VoroTop IdealGas3D -v 10
    ```

    **Command explanation:**
    - `-u [max_distance]`: Compute unnormalized Voronoi pair correlation function
    - `-v [max_distance]`: Compute normalized Voronoi pair correlation function  
    - `max_distance`: Maximum topological distance to analyze (default: 10)
    - `-2`: Specify two-dimensional analysis (for 2D systems only)

4.  **Analyze the output:** Results are printed directly to the terminal for immediate analysis.
<br><br>

**Expected output and comprehensive analysis:**

The pair correlation function output provides rich statistical data about structural correlations. 
* **Column 1:** Topological distance *k*, (number of Voronoi boundaries to traverse)
* **Column 2:** Average number of *k*-neighbors per particle (normalized if using -v)
* **Column 3:** Sample variance of *k*-neighbor counts
* **Column 4:** Sample standard deviation of *k*-neighbor counts  
* **Column 5:** Cumulative fraction of system sampled at distance *k*; when this approaches 1.0, finite-size effects become significant and data reliability decreases.
<br><br>

---

## Results

Let's analyze representative results from each system type:

### Two-dimensional, unnormalized


#### **`VoroTop Polycrystal2D -2 -u`** 
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

**Analysis of polycrystalline structure:**
- **Perfect k=1:** Exactly 6 first neighbors (hexagonal crystal signature)
- **Crystal structure:** *k*=2 shows very close to 12 neighbors (second shell of hexagonal lattice)
- **Regular progression:** Values close to 18, 24, 30... neighbors follow hexagonal lattice geometry
- **Low variance:** Small fluctuations indicate structural regularity with deviations from perfect values due to defects
<br>

----

#### **`VoroTop Liquid2D -2 -u`** 
```text
     6         0.692675  0.832271  0.000384088
2    12.6927   2.14594   1.4649    0.00108053
3    19.9419   3.70308   1.92434   0.00217474
4    27.5332   5.36133   2.31545   0.00368548
5    35.3607   6.97912   2.6418    0.0056257
6    43.319    8.53064   2.92073   0.0080026
7    51.4057   10.1501   3.18593   0.0108232
8    59.5385   11.9619   3.45859   0.0140901
9    67.7558   13.5607   3.68248   0.0178078
10   76.0365   15.1113   3.88733   0.0219799
```

**Analysis of two-dimensional LJ liquid behavior:**
- **k=1 baseline:** Exactly 6 neighbors on average (topological requirement for 2D general position)
- **Order and growth**: Monotonic increase with topological distance, but slower than in ideal gas
- **Random structure:** No oscillations or preferred distances
- **Cumulative sampling:** Fraction sampled grows steadily, indicating good statistics
<br>

----

#### **`VoroTop IdealGas2D -2 -u`** 
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
```

**Analysis of ideal gas behavior:**
- **k=1 baseline:** Exactly 6 neighbors on average (topological requirement for 2D general position)
- **Smooth growth:** Monotonic increase with topological distance
- **Random structure:** No oscillations or preferred distances
- **Statistical regularity:** Small variance reflects large system size
- **Cumulative sampling:** Fraction sampled grows steadily, indicating good statistics
<br><br>

----
----

### Three-dimensional, normalized


#### **`VoroTop HighTemperatureCu -v`** 
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
```

**Analysis of high-temperature crystalline copper:**
- **Below unity:** v(1) = 0.91 indicates fewer than ideal neighbors due to thermal effects
- **Crystalline decay:** Monotonic decrease characteristic of correlated structure
- **Thermal broadening:** Non-zero variances reflect thermal vibrations
- **Structural persistence:** Clear correlation extends to large distances
- **FCC signature:** Pattern consistent with face-centered cubic structure
<br>

----

#### **`VoroTop LiquidCu -v`** 
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

**Analysis of liquid copper:**
- **Near-unity k=1:** v(1) ≈ 0.94 close to random expectation
- **Liquid characteristics:** Faster approach to unity than crystalline system
- **Short-range order:** Some correlation visible in first few shells
- **Reduced structure:** Less long-range correlation than crystalline phase
- **Thermal motion:** Higher variances due to liquid-state dynamics
<br>

----

#### **`VoroTop IdealGas3D -v`** 
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
```

**Analysis of 3D ideal gas:**
- **Unity normalization:** v(k) ≈ 1.0 for all k (random expectation)
- **No correlations:** Minimal deviations from unity indicate complete disorder
- **Statistical fluctuations:** Small variances due to large system size
- **Baseline behavior:** Provides reference for comparison with structured systems
<br><br>


----
## Structural signatures in pair correlation functions:

Here are sevaral observations about the Voronoi pair correlaton functions as applied to various systems.

**Perfect crystals:**
- **Exact coordination:** Integer values at specific k distances
- **Geometric progression:** Predictable neighbor counts following lattice geometry
- **Low variance:** Minimal fluctuations due to structural regularity
- **Long-range order:** Correlations persist to large distances

**Defective crystals:**
- **Near-integer values:** Small deviations from perfect coordination
- **Moderate variance:** Fluctuations due to defects and thermal effects
- **Preserved trends:** Overall pattern maintained despite imperfections
- **Defect signatures:** Systematic deviations reveal defect concentrations

**Liquids:**
- **Approach to unity:** v(k) → 1 for large k (random limit)
- **Short-range structure:** Some correlation in first 2-3 shells
- **Higher variance:** Increased fluctuations due to dynamic disorder
- **Rapid decay:** Correlations disappear quickly with distance

**Glasses and amorphous systems:**
- **Intermediate behavior:** Between crystalline and liquid characteristics
- **Extended short-range order:** Correlations extend further than simple liquids
- **Structural heterogeneity:** Large variances indicate local environment diversity

----

## Conclusions

The Voronoi pair correlation function represents a powerful complement to traditional structural analysis methods, providing topology-based insights that remain robust under thermal fluctuations and coordinate uncertainties. 

Through the systematic analysis of systems ranging from perfect crystals to ideal gases, we have demonstrated how this discrete, connectivity-based approach reveals structural signatures that may be obscured in classical distance-based correlation functions. The clear distinctions between crystalline decay patterns (v(k) < 1, monotonic decrease), liquid short-range order (rapid approach to unity), and ideal gas behavior (v(k) ≈ 1 for all k) provide quantitative fingerprints for automated phase identification and structural characterization. 


<br>

---

## References

* Worlitzer, V.M., Ariel, G., and Lazar, E.A., "*Pair correlation function based on Voronoi topology*", [arXiv](https://arxiv.org/abs/2210.09731), [Phys. Rev. E 108:064115](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.108.064115), 2023. 

* Lazar, E.A., Lu, J., Rycroft, C.H., and Schwarcz, D., "*Characterizing structural features of two-dimensional particle systems through Voronoi topology*", [arXiv](https://arxiv.org/abs/2406.00553), [Model. Simul. Mater. Sci. Eng. 32:085022](https://iopscience.iop.org/article/10.1088/1361-651X/ad8ad9), 2024. 

* Lazar, E.A. "*VoroTop: Voronoi Cell Topology Visualization and Analysis Toolkit*", [arXiv](https://arxiv.org/abs/1804.04221), [Model. Simul. Mater. Sci. Eng. 26:1](https://iopscience.iop.org/article/10.1088/1361-651X/aa9a01), 2017.