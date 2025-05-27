# *VoroTop* Tutorial Through Examples

**Welcome to the comprehensive hands-on tutorial for *VoroTop*—the leading software package for Voronoi topology analysis of particle systems.**

Through carefully designed examples progressing from fundamental concepts to advanced applications, you will master the complete workflow of structural analysis using Voronoi topology. This tutorial teaches both the theoretical foundations and practical skills needed to analyze complex materials systems, from simple crystals to realistic polycrystalline structures with defects.

## What is *VoroTop*?

*VoroTop* is a powerful open-source command-line program that uses **Voronoi cell topology** to analyze local structure in particle systems. Unlike traditional approaches that rely on geometric measurements sensitive to thermal vibrations and coordinate errors, Voronoi topology provides a robust, mathematically rigorous framework for:

- **Automated defect identification** in crystalline materials
- **Quantitative structural characterization** of disordered systems  
- **Grain boundary analysis** in polycrystalline materials
- **Phase identification** and transition studies
- **High-temperature system analysis** where conventional methods fail

### Why Voronoi Topology?

**Robustness:** Topological descriptions remain stable under small coordinate perturbations, making them ideal for analyzing realistic systems with thermal vibrations or measurement noise.

**Universality:** The same mathematical framework applies to any particle system—crystals, polycrystals, liquids, glasses, colloidal systems, and granular materials.

**Automation:** Enables systematic, reproducible analysis of large datasets without manual intervention or subjective parameter choices.

**Physical insight:** Directly connects local atomic arrangements to macroscopic material properties and behavior.

## Getting Started

### Prerequisites

Before beginning this tutorial, ensure you have:

- ***VoroTop* installed** on your system ([installation instructions](../README.md))
- **Basic command-line familiarity** for running programs and navigating directories
- **Understanding of materials science concepts** (atoms, crystals, defects) at an introductory level
- **Optional:** Visualization software like [OVITO](https://www.ovito.org/) for viewing results

### Tutorial Philosophy

This tutorial follows a **progressive learning approach**:

1. **Build foundational understanding** through simple, clear examples
2. **Introduce complexity gradually** with realistic materials systems
3. **Connect theory to practice** through hands-on analysis of actual data
4. **Develop expertise** in advanced techniques and research applications

Each example is **self-contained** but builds upon previous knowledge, allowing you to either work through the complete sequence or jump to specific topics of interest.

## Tutorial Examples

### Core Methodology (Examples 1-4)

These foundational examples teach the essential concepts and techniques underlying all Voronoi topology analysis:

**[Example 1: Basic Voronoi topology calculations](Example01-Basic/README.md)**
* **Foundation concepts:** Learn what Voronoi cells are and how their topology describes local structure
* **Hands-on computation:** Calculate p-vectors and Weinberg vectors for individual particles
* **Data interpretation:** Understand symmetry, chirality, and topological descriptors
* **System comparison:** Analyze ordered vs. disordered arrangements
* **Key skills:** Basic *VoroTop* operation, output interpretation, structural insight

**[Example 2: Voronoi topology analysis through filter files](Example02-FilterAnalysis/README.md)**
* **Automated classification:** Transform raw topology data into meaningful structural categories  
* **Filter concepts:** Understand how topology families represent physical structures
* **Practical application:** Distinguish crystals from defects in realistic systems
* **File formats:** Master filter file structure and creation
* **Key skills:** Structural classification, defect identification, automation

**[Example 3: Calculate the distribution of Voronoi topologies](Example03-Distributions/README.md)**
* **Statistical analysis:** Characterize entire systems through topology distributions
* **Phase identification:** Distinguish crystals, liquids, and glasses by their signatures
* **Quantitative metrics:** Calculate structural order parameters and crystallinity measures
* **Comparative analysis:** Compare different materials, temperatures, and processing conditions
* **Key skills:** Statistical characterization, phase analysis, quantitative metrics

**[Example 4: Compute the Voronoi pair correlation function](Example04-PairCorrelations/README.md)**
* **Spatial correlations:** Analyze how structural features organize across length scales
* **Advanced theory:** Understand topology-based correlation functions vs. classical approaches
* **Structure-property connections:** Relate local topology to material behavior
* **Temperature effects:** Study how thermal energy affects structural correlations
* **Key skills:** Correlation analysis, length scale characterization, advanced theory

### Visualization and Analysis Tools (Examples 5-6)

These examples teach essential techniques for visualizing results and identifying structural patterns:

**[Example 5: Generate EPS files for visualization](Example05-PostScriptImages/README.md)**
* **Publication-quality graphics:** Create high-resolution vector graphics for papers and presentations
* **Multiple coloring schemes:** Highlight different structural features effectively
* **Visual analysis:** Extract insights through strategic visualization approaches
* **Figure preparation:** Professional workflows for scientific communication
* **Key skills:** Scientific visualization, figure preparation, visual analysis

**Example 6: Cluster analysis**
* **Multi-particle defects:** Identify extended structural features like grain boundaries and dislocations
* **Automated grouping:** Systematically organize particles into meaningful structural units
* **Quantitative characterization:** Measure defect sizes, distributions, and concentrations
* **Crystal nucleation:** Analyze phase transformation and solidification processes
* **Key skills:** Cluster identification, defect analysis, phase transformation studies

### Advanced Techniques (Examples 7-10)

These examples cover sophisticated analysis methods for research applications and complex materials systems:

**Example 7: Create filters from ideal structures**
* **Systematic filter generation:** Use perturbation analysis to create comprehensive, robust filters
* **Theoretical foundation:** Bridge ideal structures and realistic finite-temperature systems
* **Parameter optimization:** Develop filters tailored to specific materials and conditions
* **Quality validation:** Ensure filter accuracy and completeness through systematic testing
* **Key skills:** Filter development, perturbation analysis, systematic methodology

**Example 8: Identify defects in crystals using filters**
* **Comprehensive defect analysis:** Apply advanced filters to identify and characterize all defect types
* **Quantitative assessment:** Measure defect concentrations, spatial distributions, and interactions
* **Materials quality control:** Evaluate crystal perfection and processing effectiveness
* **Structure-property relationships:** Connect defect populations to material performance
* **Key skills:** Defect characterization, quality assessment, materials evaluation

**Example 9: Characterize grain boundaries using Voronoi analysis**
* **Interface structure:** Analyze the atomic-scale structure of grain boundaries
* **Boundary classification:** Distinguish different grain boundary types and characters
* **Polycrystalline analysis:** Characterize grain size distributions and boundary networks
* **Processing-structure relationships:** Connect synthesis conditions to final microstructure
* **Key skills:** Interface analysis, polycrystalline characterization, microstructure analysis

**Example 10: Advanced time-dependent analysis**
* **Dynamic processes:** Track structural evolution during phase transformations, deformation, or growth
* **Kinetic analysis:** Measure transformation rates and pathway mechanisms
* **Process optimization:** Identify critical processing parameters affecting final structure
* **Predictive modeling:** Develop structure-based models for material behavior
* **Key skills:** Dynamic analysis, process modeling, predictive capabilities



## Key Features and Capabilities

### What You'll Learn:

✅ **Fundamental theory** of Voronoi topology and its applications to materials science  
✅ **Practical skills** for analyzing real materials data using *VoroTop*  
✅ **Advanced techniques** for research-quality structural characterization  
✅ **Visualization methods** for scientific communication and insight generation  
✅ **Automation strategies** for high-throughput materials analysis  
✅ **Quality assessment** methods for validating results and ensuring accuracy  
✅ **Research applications** spanning crystallography, metallurgy, and materials engineering  

### System Types Covered:

- **Perfect crystals** (hexagonal, square, cubic lattices)
- **Defective crystals** (vacancies, interstitials, dislocations)
- **Polycrystalline materials** (grain boundaries, texture, size distributions)
- **Disordered systems** (liquids, glasses, ideal gases)
- **Interfaces and surfaces** (grain boundaries, phase boundaries)
- **Dynamic systems** (phase transformations, mechanical deformation)

## Running the Examples

Each example follows a consistent structure for easy navigation:

### Standard Example Format:

```bash
cd VoroTop/Tutorial/ExampleXX-TopicName
VoroTop input_file.dump [options]
# Analysis and interpretation of results
```

### File Organization:
- **README.md** - Complete tutorial with theory, instructions, and analysis
- **Input files** - Representative data for hands-on analysis
- **Reference outputs** - Expected results for validation
- **Scripts** (where applicable) - Automation tools and advanced workflows

### Getting Help:

- **Command-line help:** `VoroTop -h` for complete option reference
- **Example-specific guidance:** Each README.md contains detailed instructions
- **Troubleshooting:** Common issues and solutions documented in each example
- **Community support:** [*VoroTop* repository](https://gitlab.com/mLazar/VoroTop/) for questions and discussions

## Advanced Applications

Beyond the tutorial examples, *VoroTop* enables cutting-edge research in:

**🔬 Materials Discovery:** High-throughput screening of new materials and structures  
**🏭 Process Optimization:** Quantitative relationships between processing and microstructure  
**🔍 Quality Control:** Automated defect detection and materials characterization  
**🧪 Phase Behavior:** Advanced studies of phase transitions and stability  
**💻 Method Development:** New algorithms for structural analysis and pattern recognition  
**🌡️ Extreme Conditions:** Analysis of materials under high temperature, pressure, or deformation  

## Contributing and Community

This tutorial represents an evolving resource that grows with the *VoroTop* community:

- **Share your applications:** Examples of *VoroTop* use in your research
- **Suggest improvements:** Tutorial enhancements and additional examples  
- **Report issues:** Help improve tutorial clarity and accuracy
- **Contribute code:** Advanced analysis scripts and automation tools

## Citation and Acknowledgments

Further details, theoretical and practical, about Voronoi topology analysis of particle systems can be found in the following papers.  

* Lazar, E.A. "*VoroTop: Voronoi Cell Topology Visualization and Analysis Toolkit*", [arXiv](https://arxiv.org/abs/1804.04221), [Model. Simul. Mater. Sci. Eng. 26:1](https://iopscience.iop.org/article/10.1088/1361-651X/aa9a01), 2017.

* Lazar, E.A., Lu, J., Rycroft, C.H., and Schwarcz, D., "*Characterizing structural features of two-dimensional particle systems through Voronoi topology*", [arXiv](https://arxiv.org/abs/2406.00553), [Model. Simul. Mater. Sci. Eng. 32:085022](https://iopscience.iop.org/article/10.1088/1361-651X/ad8ad9), 2024. 

* Lazar, E.A., Han, J. and Srolovitz, D.J. "*A Topological Framework for Local Structure Analysis in Condensed Matter*", [arXiv](https://arxiv.org/abs/1508.05937), [Proc. Natl. Acad. Sci. 112:E5769](https://www.pnas.org/doi/10.1073/pnas.1505788112), 2015.

* Landweber, P.S., Lazar, E.A., Patel, N. "*On fiber diameters of continuous maps*", [arXiv](https://arxiv.org/abs/1503.07597), [Amer. Math. Monthly 123:4](https://www.jstor.org/stable/10.4169/amer.math.monthly.123.4.392), 2016.

* Lazar, E.A., Mason, J.K., MacPherson, R.D. and Srolovitz, D.J. "*Complete topology of cells, grains, and bubbles in three-dimensional microstructures*", [arXiv](https://arxiv.org/abs/1207.5054), [Phys. Rev. Lett. 109:095505](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.095505), 2012.

If you use *VoroTop* or this tutorial in your research, please consider citing one or some of these papers.


This README was written with assistance from Claude (Anthropic).

---

**Ready to begin your journey into advanced structural analysis?**  
**Start with [Example 1: Basic Voronoi topology calculations](Example01-Basic/README.md) and discover the power of topology-based materials characterization!**
