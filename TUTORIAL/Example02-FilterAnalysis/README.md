# Example 2: Voronoi topology analysis through filter files

**Learning objectives:**

* After completing this example, you will:
  * Understand how filter files enable automated structural classification in particle systems.
  * Learn the three-part structure of filter files and how to interpret their contents.
  * Apply filters to distinguish between crystalline regions and structural defects.
  * Interpret *VoroTop* output files containing structural classification information.
  * Recognize how filters transform raw topological data into meaningful structural categories.
  * Understand the concept of families of Voronoi topologies and their relationship to physical structures.
<br><br>

----

## Background:

Raw Voronoi topology data, while mathematically precise, requires interpretation to become physically meaningful. A *filter* provides this interpretation by defining which Voronoi topologies correspond to which structural features. Think of a filter as a "dictionary" that translates the mathematical language of p-vectors into the physical language of crystals, defects, and other structural elements.

**The concept of structural families:**
In perfect crystals, the Voronoi cell of every particle is identical, or else one of several kinds corresponding to different particles in the repeating unit cell. In real systems, however, thermal vibrations and small strains mean that a single structural feature (like a crystal or a specific type of defect) can be associated with multiple related Voronoi topologies. 

For example, in a perfect square lattice, the Voronoi cell of each particle is a square.  When particle positions are perturbed, as might occur in a finite-temperature simulation, the geometry and topology of Voronoi cells change.  However, even in such systems, the topologies of the particles belong to a related "family" that we can associate with the crystalline state. Filters capture these families of topologies and assign them to meaningful structural categories.

**Why filters are essential:**
- **Automation:** Filters enable automated analysis of large systems without manual inspection
- **Consistency:** They provide reproducible structural classifications across different systems
- **Physical insight:** They connect mathematical topology to physical structure
- **Flexibility:** Different filters can emphasize different aspects of the same system
- **Comparison:** They enable quantitative comparison between different systems or simulation conditions
<br><br>

**Filter applications:**
Filters are particularly powerful for identifying:
- **Crystal phases:** Different crystal structures have characteristic Voronoi topology families
- **Structural defects:** Vacancies, dislocations, and interstitials each have recognizable topological signatures
- **Interfaces:** Grain boundaries and phase boundaries create distinctive local topologies
- **Disorder:** Non-crystalline regions can be identified by exclusion from crystalline families
<br><br>

----

**Input data:**

* **Data file:** `Polycrystal2D` - LAMMPS dump format with particle coordinates
* **Filter file:** `crystal_2d.filter` - Defines structural classifications for 2D hexagonal crystals
* **System description:** A two-dimensional polycrystalline system created by cooling a Lennard-Jones liquid. The system contains multiple crystal grains separated by grain boundaries, with occasional point defects (vacancies). This realistic system demonstrates how filters work on complex, heterogeneous structures rather than idealized perfect crystals.
<br><br>

----

**Understanding the filter file format:**

Before running the analysis, let's examine the filter file structure:

```text
#       Filter with unique type in two-dimensional defect-free
#       triagonal crystal, in which all particles are hexagons 
#       surrounded by hexagons.
*       1       Crystal         
1       (6,6,6,6,6,6,6)
```

A filter file is composed of three parts:

**Part 1 - Comments (lines beginning with #):**
- Provide documentation about the filter's purpose and contents
- Describe the source data or theoretical basis
- Include metadata like creation date, author, or version information
- Help users understand when and how to apply the filter

**Part 2 - Structure type definitions (lines beginning with \*):**
- Define user-readable names for each structural category
- Format: `* [index] [name]`
- Indices must be consecutive integers starting from 1
- Names should be descriptive (Crystal, GrainBoundary, Vacancy, etc.)

**Part 3 - Topology assignments:**
- Map specific Voronoi topologies to structure types
- Format: `[structure_type_index] [p-vector_or_Weinberg_vector]`
- Each line assigns one topology to one structure type
- Multiple topologies can be assigned to the same structure type (families)

**Understanding this specific filter:**
This simple filter defines only one structure type called "Crystal" and associates it with the p-vector `(6,6,6,6,6,6,6)`. This topology corresponds to particles in perfect hexagonal crystals—particles with six neighbors, each of which also has six neighbors.
<br><br>

----

**Steps:**

1.  **Navigate to the example directory:** Open your terminal or command prompt and navigate to the example directory.

    ```bash
    cd VoroTop/Tutorial/Example02-FilterAnalysis
    ```

2.  **Examine the filter file:** Before running the analysis, look at the filter contents:

    ```bash
    cat crystal_2d.filter
    ```

3.  **Execute the *VoroTop* command:** Run the following command to apply the filter:

    ```bash
    VoroTop Polycrystal2D -2 -f crystal_2d.filter -l
    ```

    **Command explanation:**
    - `Polycrystal2D`: Input file containing particle coordinates
    - `-2`: Specify two-dimensional analysis
    - `-f crystal_2d.filter`: Load and apply the specified filter file
    - `-l`: Output results in LAMMPS dump format for easy visualization

4.  **Examine the output file:** A file named `Polycrystal2D.dump` will be created containing the original data plus structural classifications.
<br><br>

----

**Expected output and explanation:**

The output file `Polycrystal2D.dump` is an enhanced LAMMPS dump file with additional structural information:

```text
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
18000
ITEM: BOX BOUNDS pp pp pp
0 100
0 100
0 1
ITEM: ATOMS id x y z vt
16210   64.538    6.69981   0   1       
1898    66.8254   5.30426   0   0       
16205   67.4699   6.20108   0   0       
16530   65.0532   22.4752   0   1       
13199   45.799    6.69556   0   1       
```

**New column - Voronoi Topology (vt):**
* **vt = 1:** Particles classified as "Crystal" (matching the p-vector (6,6,6,6,6,6,7) from the filter)
* **vt = 0:** Particles with Voronoi topologies NOT found in the filter (structural defects, grain boundaries, etc.)


-----

**Physical interpretation:**

**Crystal particles (vt = 1):**
- Have local environments closely matching perfect hexagonal crystals
- Form the bulk regions within individual crystal grains
- Represent the majority of particles in well-crystallized regions
- Maintain their classification despite thermal vibrations due to topological robustness

**Defect particles (vt = 0):**
- Include all particles whose local topology differs from perfect crystal
- **Grain boundary particles:** Located at interfaces between different crystal orientations
- **Point defects:** Vacancies, interstitials, and their immediate neighbors
- **Thermal defects:** Particles temporarily perturbed from ideal positions
- **Edge effects:** Particles near system boundaries with incomplete coordination

**Visual analysis using OVITO or similar tools:**
When visualized, you should observe:
- **Large regions of vt=1 particles** forming coherent crystal grains
- **Networks of vt=0 particles** outlining grain boundaries
- **Small clusters of vt=0 particles** marking point defects
- **Clear structural organization** with physical meaning

**Quantitative analysis:**
You can analyze the output to determine:
- **Crystallinity fraction:** Ratio of vt=1 to total particles
- **Grain sizes:** Size distribution of connected vt=1 regions  
- **Defect density:** Concentration of vt=0 particles
- **Structural quality:** How closely the system approximates perfect crystal structure
<br><br>

----

**Understanding filter limitations and extensions:**

**Current filter limitations:**
This simple filter has several limitations that illustrate important concepts:

1. **Binary classification:** Only distinguishes "crystal" vs "everything else"
2. **Strict topology matching:** Requires exact match to the ideal p-vector
3. **No defect characterization:** Cannot distinguish between different types of defects
4. **Single crystal structure:** Only recognizes one type of crystalline order
<br><br>

----

**Advanced filter concepts:**

**Multi-type filters:**
More sophisticated filters can include multiple structure types:
```text
*       1       Crystal        
*       2       GrainBoundary  
*       3       Dislocation    
*       4       Vacancy        
*       5       Interstitial   
1       (6,6,6,6,6,6,6)
2       (5,6,6,7,6,7)
2       (7,5,6,6,5,6,6,6)
2       (5,6,6,6,6,7)
2       (7,5,6,6,6,6,6,6)
3       (5,6,6,6,6,7)
3       (7,5,6,6,6,6,6,6)
4       (5,6,6,6,7,7)           
4       (7,5,6,6,6,5,7,7)       
4       (6,5,6,6,6,7,7)         
4       (7,5,6,6,6,6,7,6)       
4       (5,6,6,6,6,7)           
4       (5,6,6,6,6,8)           
4       (8,5,6,6,6,5,6,6,6)     
4       (6,5,6,6,6,6,8)         
4       (6,6,6,6,6,6,8)         
5       (6,5,7,5,7,5,7)
```

**Topology families:**
Notice that a single structure type can be associated with multiple Voronoi topologies.

**Indeterminate types:**
Notice also that some topologies belong to multiple structure types, requiring special handling for disambiguation.
<br><br>

----

**Building more comprehensive filters:**

**Step 1: Identify target structures**
- Determine which structural features you want to classify
- Consider both bulk phases and defect types
- Decide on the level of detail needed

**Step 2: Determine characteristic topologies**
- Use *VoroTop* with `-vt` option to survey topology distributions
- Identify topologies associated with each structural feature
- Consider thermal effects and topology families

**Step 3: Construct and test the filter**
- Create filter file with appropriate structure types
- Test on known systems to validate classifications
- Refine based on results and physical expectations
<br><br>

----

**Connecting to subsequent examples:**

This filter-based analysis provides the foundation for more advanced techniques:

- **Example 3** will show how to analyze distributions of classified structures
- **Example 4** will explore spatial correlations between different structure types  
- **Example 5** will demonstrate visualization of filter-classified structures
- **Example 6** will use filter classifications as input for cluster analysis
- **Example 7** will create filters from ideal structures through perturbation analysis.

The structural classifications generated here transform raw topological data into physically meaningful categories that enable quantitative analysis of structure-property relationships in complex materials systems.
<br><br>

**Exercise suggestions:**

1. **Modify the filter** to include additional structure types for grain boundaries
2. **Apply the filter** to different input systems and compare crystallinity fractions
3. **Analyze the spatial distribution** of classified particles using visualization tools
4. **Calculate quantitative metrics** like average grain size or defect density
5. **Compare results** with other structural analysis methods
<br><br>

## Conclusions

Filter-based analysis represents a fundamental transformation from descriptive to predictive structural characterization, converting the mathematical language of Voronoi topology into the physical language of crystals, defects, and material properties. 

Through the systematic application of filters, we have demonstrated how complex polycrystalline systems can be automatically parsed into constituent structural elements—crystalline regions, grain boundaries, point defects, and disordered phases—without subjective parameter choices or manual intervention. The robustness of topological classification ensures that these structural assignments remain stable under thermal vibrations, small coordinate perturbations, and measurement uncertainties that typically plague geometric approaches. Most significantly, filters enable the transition from qualitative visual inspection to quantitative statistical analysis. 

The filter framework's flexibility allows adaptation to any crystal system, defect type, or materials application, while its automation capabilities make possible the high-throughput analysis of large datasets required for modern materials discovery and process optimization. This approach establishes the foundation for all subsequent *VoroTop* analyses, transforming particle coordinates into physically meaningful structural information that directly connects atomic-scale arrangements to macroscopic material behavior.


## References

* Lazar, E.A. "*VoroTop: Voronoi Cell Topology Visualization and Analysis Toolkit*", [arXiv](https://arxiv.org/abs/1804.04221), [Model. Simul. Mater. Sci. Eng. 26:1](https://iopscience.iop.org/article/10.1088/1361-651X/aa9a01), 2017.

* Lazar, E.A., Lu, J., Rycroft, C.H., and Schwarcz, D., "*Characterizing structural features of two-dimensional particle systems through Voronoi topology*", [arXiv](https://arxiv.org/abs/2406.00553), [Model. Simul. Mater. Sci. Eng. 32:085022](https://iopscience.iop.org/article/10.1088/1361-651X/ad8ad9), 2024.