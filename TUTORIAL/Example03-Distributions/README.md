# Example 3: Calculate the distribution of Voronoi topologies

**Learning objectives:**

* After completing this example, you will:
  * Understand how Voronoi topology distributions characterize the structural nature of particle systems.
  * Compute and interpret topology distributions for both ordered and disordered systems.
  * Analyze the relationship between system structure and topology distribution patterns.
  * Use distribution data to distinguish between different phases and structural states.
  * Understand the role of chirality and symmetry in structural characterization.
  * Apply statistical analysis to quantify structural order and disorder.
<br><br>


## Background:

The distribution of Voronoi topologies in a system provides a powerful "structural fingerprint" that characterizes the nature of particle arrangements. While individual topologies describe local neighborhoods, the statistical distribution reveals global structural properties and can distinguish between fundamentally different types of systems.

**Why distributions matter:**
- **Phase identification:** Different phases (crystal, liquid, glass) have characteristic topology distributions
- **Order quantification:** The breadth and shape of distributions reflect structural order vs. disorder
- **Thermal effects:** Temperature influences topology distributions in predictable ways
- **Comparative analysis:** Distributions enable quantitative comparison between different systems
- **Structural evolution:** Changes in distributions track structural transformations over time

**Distribution characteristics:**

**Crystalline systems:**  
In perfect crystals, all particles have identical local environments, producing distributions with sharp peaks at specific topologies. Even with thermal vibrations, crystalline systems show narrow distributions dominated by a few related topologies from the same structural family.

**Liquid systems:**  
Liquids exhibit broad distributions with many different topologies present. However, certain topologies may still be more common due to local packing preferences and energetic factors. The distribution reflects the balance between entropy (which favors diversity) and local energy minimization.

**Disordered systems:**  
Truly random systems (like ideal gases) show very broad distributions where many topologies appear with relatively similar frequencies. The distribution shape reflects the underlying geometric constraints and statistical mechanical principles governing random packing.

**Glass systems:**  
Glasses often show intermediate behavior—broader than crystals but with some preferred local arrangements that create peaks in the distribution, reflecting their "frustrated" structure between order and disorder.
<br><br>

----

**Input data:**

We analyze six different systems that span the spectrum from perfect order to complete disorder:

### Two-dimensional systems:
* **`IdealGas2D`** - Random point distribution (complete disorder)
* **`Liquid2D`** - A Lennard-Jones liquid 
* **`Polycrystal2D`** - Hexagonal polycrystal with grain boundaries (imperfect order)

### Three-dimensional systems:
* **`HighTemperatureCu`** - Copper crystal at 85% of melting temperature (thermally excited crystal)
* **`LiquidCu.dump`** - Copper heated to 120% of melting temperature (liquid state)
* **`IdealGas3D.dump`** - Random point distribution in 3D (complete disorder)

**System sizes:** All systems contain 10,000-108,000 particles to ensure statistical significance in the topology distributions.
<br><br>

**Steps:**

1.  **Navigate to the example directory:**

    ```bash
    cd VoroTop/Tutorial/Example03-Distributions
    ```

2.  **Execute distribution analysis for 2D systems:**

    ```bash
    VoroTop IdealGas2D -2 -d
    VoroTop Liquid2D -2 -d
    VoroTop Polycrystal2D -2 -d
    ```

3.  **Execute distribution analysis for 3D systems:**

    ```bash
    VoroTop HighTemperatureCu -d
    VoroTop LiquidCu.dump -d
    VoroTop IdealGas3D.dump -d
    ```

    **Command explanation:**
    - `-2`: Specify two-dimensional analysis (only for 2D systems)
    - `-d`: Calculate and output the distribution of Voronoi topologies

4.  **Examine the output files:** Distribution files named `<datafile>.distribution` will be created for analysis.
<br><br>

**Expected output and detailed analysis:**

Each distribution file contains a header followed by topology data with frequency and chirality information as follows:
* **Column 1:** Voronoi topology (p-vector for 2D, Weinberg vector for 3D)
* **Column 2:** Total count of particles with this topology
* **Column 3:** Count with left-handed chirality  
* **Column 4:** Count with no chirality (symmetric arrangements)
* **Column 5:** Count with right-handed chirality

**Note:** Columns 3 + 4 + 5 = Column 2 (total count)
<br><br>



----

### Two-dimensional results analysis:

#### **`IdealGas2D.distribution`** - Complete disorder
```text
#       DISTRIBUTION CREATED FROM IdealGas2D
#       Total particles sampled: 17280
#       Total Voronoi topologies observed: 7809
#       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
(4,6,6,7,8)     66      41      0       25
(4,5,6,7,7)     61      33      0       28
(4,5,7,6,8)     58      27      0       31
(4,6,7,7,8)     58      28      0       30
(5,5,6,7,6,7)   56      28      0       28
(5,5,7,6,7,7)   56      30      0       26
(4,6,6,7,7)     55      0       55      0
(4,6,6,6,7)     52      0       52      0
```

**Analysis of ideal gas distribution:**
- **Broad diversity:** Many different topologies with similar frequencies
- **No dominant topology:** Most common topology appears only 66 times out of 17280 total particles
- **Chirality balance:** Left-handed and right-handed forms appear in roughly equal numbers, and most arrangements of particles lack mirror symmetry
- **Statistical behavior:** Distribution reflects random geometric arrangements

#### **`Liquid2D.distribution`** 
```text
#       DISTRIBUTION CREATED FROM Liquid2D
#       Total particles sampled: 18225
#       Total Voronoi topologies observed: 1802
#       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
(6,5,6,6,6,6,7) 528     262     0       266
(5,6,6,6,6,7)   499     0       499     0
(5,6,6,7,6,7)   404     0       404     0
(6,5,6,6,6,7,6) 379     187     0       192
(6,6,6,6,6,6,7) 379     0       379     0
(6,5,6,6,6,6,6) 348     0       348     0
(5,6,6,6,7,7)   344     0       344     0
```

**Analysis of two-dimensional LJ:**
- **Intermediate diversity:** Fewer unique topologies than ideal gas
- **Emerging preferences:** Some topologies significantly more common
- **Symmetry effects:** Many topologies are non-chiral (thermal perturbations of symmetric arrangements)

#### **`Polycrystal2D.distribution`** - Crystalline order with defects
```text
#       DISTRIBUTION CREATED FROM Polycrystal2D
#       Total particles sampled: 17280
#       Total Voronoi topologies observed: 31
#       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
(6,6,6,6,6,6,6) 15707   0       15707   0
(6,6,6,6,6,6,7) 427     0       427     0
(6,5,6,6,6,6,7) 336     173     0       163
(6,5,6,6,6,6,6) 205     0       205     0
(5,6,6,6,6,7)   138     0       138     0
```

**Analysis of hexagonal polycrystal:**
- **Dominant topology:** (6,6,6,6,6,6,7) represents ~87% of particles (perfect hexagonal crystal)
- **Defect topologies:** All other topologies represent grain boundaries, point defects, or thermal perturbations
- **Structural quality:** High crystallinity evident from single dominant peak
- **Defect characterization:** Minority topologies correspond to defect types
<br><br>


----

### Three-dimensional results analysis:

#### **`HighTemperatureCu.distribution`** - High-temperature crystal
```text
#       DISTRIBUTION CREATED FROM HighTemperatureCu
#       Total particles sampled: 108000
#       Total Voronoi topologies observed: 775
#       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
(1,2,3,4,1,4,5,6,7,1,7,8,9,2,9,10,11,12,3,12,13,14,5,14,15,16,17,6,17,18,19,8,19,20,21,10,21,22,11,22,23,24,13,24,15,24,23,25,16,25,26,18,26,20,26,25,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)      9868    4946    0       4922
(1,2,3,4,1,4,5,6,7,1,7,8,9,2,9,10,11,12,3,12,13,14,5,14,15,16,6,16,17,18,8,18,19,20,10,20,21,11,21,22,13,22,23,15,23,24,17,24,19,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)        7272    3662    0       3610
(1,2,3,4,1,4,5,6,7,1,7,8,9,10,2,10,11,12,3,12,13,14,5,14,15,16,6,16,17,18,8,18,19,9,19,20,21,11,21,22,13,22,23,24,15,24,17,24,23,20,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)        6329    0       6329    0
```

**Analysis of high-temperature copper:**
- **FCC signature:** Top topologies correspond to face-centered cubic crystal structure
- **Thermal broadening:** Multiple related topologies due to thermal vibrations
- **Crystalline character:** Despite high temperature, clear crystalline peaks dominate
- **Chirality patterns:** Some topologies show chirality, others are symmetric

#### **`LiquidCu.distribution`** - Liquid copper
```text
#       DISTRIBUTION CREATED FROM LiquidCu
#       Total particles sampled: 108000
#       Total Voronoi topologies observed: 30764
#       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
(1,2,3,4,1,4,5,6,7,1,7,8,9,10,2,10,11,12,3,12,13,14,5,14,15,16,6,16,17,8,17,18,19,9,19,20,11,20,21,13,21,22,15,22,18,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)  1457    0       1457    0
(1,2,3,4,1,4,5,6,7,1,7,8,9,2,9,10,11,12,3,12,13,14,5,14,15,16,6,16,17,18,8,18,19,10,19,20,21,11,21,13,21,20,22,15,22,17,22,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)  1349    655     0       694
(1,2,3,4,1,4,5,6,7,1,7,8,9,2,9,10,11,12,3,12,13,14,5,14,15,16,6,16,17,8,17,18,10,18,19,20,11,20,13,20,19,15,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)    924     0       924     0
```

**Analysis of liquid copper:**
- **Broader distribution:** Many more topologies with smaller individual counts
- **Liquid characteristics:** Loss of long-range crystalline order
- **Local structure:** Some topologies still preferred due to local packing constraints
- **Reduced crystallinity:** No single topology dominates like in crystalline systems

#### **`IdealGas3D.distribution`** - Three-dimensional disorder
```text
#       DISTRIBUTION CREATED FROM IdealGas3D
#       Total particles sampled: 108000
#       Total Voronoi topologies observed: 80998
#       Columns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types
(1,2,3,1,3,4,5,6,1,6,7,8,9,2,9,10,11,4,11,12,5,12,13,7,13,14,8,14,10,14,13,12,11,10,9,8,7,6,5,4,3,2,1)  320     162     0       158
(1,2,3,1,3,4,5,6,1,6,7,8,2,8,9,10,4,10,11,5,11,12,7,12,9,12,11,10,9,8,7,6,5,4,3,2,1)    172     0       172     0
(1,2,3,4,1,4,5,6,1,6,7,8,9,2,9,10,11,3,11,12,13,5,13,14,7,14,15,8,15,16,10,16,12,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)        156     79      0       77
```

**Analysis of 3D ideal gas:**
- **Maximum diversity:** Each topology appears very infrequently
- **Random character:** Distribution reflects purely geometric constraints
- **Chirality balance:** Chiral topologies show roughly equal left/right distributions
- **Statistical significance:** Large system size ensures reliable statistics

<br>

## Conclusions

The distribution of Voronoi topologies provides a powerful statistical fingerprint that transforms complex particle arrangements into quantifiable structural signatures. 

Through the analysis of systems ranging from perfect crystals to ideal gases, we have seen how the breadth, shape, and dominant peaks of topology distributions directly reflect the underlying physical processes governing particle organization. In crystalline systems, sharp peaks at specific topologies reveal short- and long-range order, while broad, diverse distributions in disordered systems reflect the balance between entropic randomness and local energetic constraints. 

These statistical measures enable automated phase identification, quantitative assessment of structural quality, and direct comparison between different materials, processing conditions, and simulation parameters. Most importantly, topology distributions provide a robust, parameter-free approach to structural characterization that remains meaningful across temperature ranges where traditional geometric methods fail, making them essential tools for understanding structure-property relationships in realistic materials systems.
<br><br>


## References

* Lazar, E.A. "*VoroTop: Voronoi Cell Topology Visualization and Analysis Toolkit*", [arXiv](https://arxiv.org/abs/1804.04221), [Model. Simul. Mater. Sci. Eng. 26:1](https://iopscience.iop.org/article/10.1088/1361-651X/aa9a01), 2017.

* Lazar, E.A., Lu, J., Rycroft, C.H., and Schwarcz, D., "*Characterizing structural features of two-dimensional particle systems through Voronoi topology*", [arXiv](https://arxiv.org/abs/2406.00553), [Model. Simul. Mater. Sci. Eng. 32:085022](https://iopscience.iop.org/article/10.1088/1361-651X/ad8ad9), 2024.