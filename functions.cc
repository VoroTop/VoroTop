////////////////////////////////////////////////////////
////                                                ////
////   ******************************************   ////
////   *                                        *   ////
////   *     VoroTop: Voronoi Cell Topology     *   ////
////   *   Visualization and Analysis Toolkit   *   ////
////   *             (Version 1.0)              *   ////
////   *                                        *   ////
////   *           Emanuel A. Lazar             *   ////
////   *          Bar Ilan University           *   ////
////   *               June 2024                *   ////
////   *                                        *   ////
////   ******************************************   ////
////                                                ////
////////////////////////////////////////////////////////

////   File: functions.cc


#include <iostream>
#include <stdexcept>

#include "filters.hh"
#include "variables.hh"



////////////////////////////////////////////////////
////
////   PRINT HELP MESSAGE WITH BASIC INSTRUCTIONS
////
////////////////////////////////////////////////////

void help_message(void) {
    std::cout <<
    "   ********************************************\n"
    "   *                                          *\n"
    "   *      VoroTop: Voronoi Cell Topology      *\n"
    "   *    Visualization and Analysis Toolkit    *\n"
    "   *              (Version 1.0)               *\n"
    "   *                                          *\n"
    "   *            Emanuel A. Lazar              *\n"
    "   *           Bar Ilan University            *\n"
    "   *                June 2024                 *\n"
    "   *                                          *\n"
    "   ********************************************\n"
    "                                               \n"
    "Syntax: VoroTop <filename> [options]           \n"
    "                                               \n"
    "         VoroTop reads an input file and computes the complete Voronoi cell      \n"
    "         topology for each point.  Output is determined by options.              \n"
    "                                                                                 \n"
    "Supported file formats (auto-detected):                                          \n"
    "         LAMMPS dump, LAMMPS data, Extended XYZ, VASP POSCAR/CONTCAR, and CIF.  \n"
    "         Both orthogonal and triclinic periodic cells are supported.  CIF files  \n"
    "         are processed with full symmetry expansion; lattice centering (I, F, A, \n"
    "         B, C, R) is automatically applied when not already included in the      \n"
    "         listed symmetry operations.                                             \n"
    "                                                                                 \n"
    "Available options:                                                               \n"
    "                                                                                 \n"
    " -f   : Specify filter file for classifying Voronoi topologies.  Option must be  \n"
    "        followed by the name of the filter file.                                 \n"
    "                                                                                 \n"
    " -2   : Specify that the system is two-dimensional.  If not specified, system is \n"
    "        treated as three-dimensional.                                            \n"
    "                                                                                 \n"
    " -vt  : Compute topology of each particle, save to <filename>.vectors.           \n"
    " -d   : Compute distribution of Voronoi topologies, represented as either        \n"
    "        p-vectors or Weinberg vectors, and save to <filename>.distribution.      \n"
    " -l   : Generate a LAMMPS <filename>.dump file using a filter provided by user.  \n"
    " -g   : Compute distribution of Voronoi topologies for perturbations of a system. \n"
    "        This option should be followed by two numbers: a target sample count     \n"
    "        and the perturbation magnitude (Gaussian standard deviation).  The       \n"
    "        system is replicated into a supercell to reach the target count.         \n"
    " -mf  : Generate a filter file from a crystal structure, saved to                \n"
    "        <filename>.filter.  If no additional arguments are given, the exact      \n"
    "        (unperturbed) topologies are computed.  If followed by two numbers       \n"
    "        (target sample count and perturbation size), the system is replicated    \n"
    "        into a supercell and a Gaussian perturbation is applied.  All observed   \n"
    "        topologies are collected into the filter.  A target of 1000000 is        \n"
    "        recommended for thorough sampling.                                       \n"
    "                                                                                 \n"
    " -c   : Cluster analysis.  Particles with topological types not included in the  \n"
    "        specified filter are considered to be defects.  Defect particles with    \n"
    "        adjacent Voronoi cells are considered to belong to the same defect; we   \n"
    "        call each Voronoi-contiguous set of defect particles a cluster and       \n"
    "        assign the constituent particles a unique index.  Similar cluster        \n"
    "        analysis is also performed for types in the filter.  This option         \n"
    "        requires loading a filter file for analysis.                             \n"
    " -r   : Resolve indeterminate types.                                             \n"
    "                                                                                 \n"
    " -e   : Output 2D system in EPS format.  This option can be followed by one or   \n"
    "        two arguments.  If one argument is given, it specifies the coloring      \n"
    "        scheme for the Voronoi cells.  If two arguments are given, the first     \n"
    "        specifies the coloring scheme and the second specifies the number of     \n"
    "        particles to be drawn.  The coloring scheme can be one of the following: \n"
    "           0: No particles drawn                                                   \n"       
    "           1: Particles drawn in black                                             \n"
    "           2: Particles colored according to number of edges                       \n"
    "           3: Particles colored according to filter classification                 \n"
    "           4: Particles colored according to Voronoi distance from central particle\n"
    "           5: Particles colored according to defect cluster ID                     \n"
    "           6: Particles colored according to crystal cluster ID                    \n"
    "           7: Particles colored according to defect cluster size                   \n"
    "           8: Particles colored according to crystal cluster size                  \n"
    " -n   : Do not draw Voronoi cells in EPS output; only draw particles.             \n"
    "                                                                                   \n"
    " -u   : Compute unnormalized Voronoi pair correlation function.  Takes an optional  \n"
    "        argument specifying the maximal Voronoi distance (default: 10).            \n"
    " -v   : Compute normalized Voronoi pair correlation function.  Takes an optional   \n"
    "        argument specifying the maximal Voronoi distance (default: 10).            \n"
    "        Normalization constants are available up to distance 20.                   \n"
    "                                                                                   \n"
    " -t   : Specify number of threads to use; if not specified, then this number will  \n"
    "        be set to one fewer than the number of processors available.               \n"
    "                                                                                   \n"
    " -h   : Display this help message.                                                 \n"
    "                                                                                   \n";
}



////////////////////////////////////////////////////
////
////   PARSE INPUT ARGUMENTS AND DETERMINE WHAT
////   SHOULD BE DONE
////
////////////////////////////////////////////////////

void parse_command_line_options(int argc, char *argv[])
{
    if(argc<2)
    {
        help_message();
        std::cout << "No input data file specified." << '\n';
        exit(0);
    }
    
    if(argc<3)
    {
        help_message();
        std::cout << "No command-line options specified." << '\n';
        exit(0);
    }
    
    // THE FILENAME OF THE INPUT DATA FILE IS THE FIRST ARGUMENT
    filename_data = argv[1];
    
    
    // DEFAULT THREADS IF NOT SPECIFIED BY USER
#ifdef _OPENMP
    if(omp_get_max_threads()>2) threads = omp_get_max_threads()-1;
    else                        threads = 1;
#else
    threads = 1;
#endif
    

    // DEFAULT DIMENSION 3; CAN BE CHANGED TO 2 THROUGH COMMAND-LINE FLAG -2
    dimension = 3;

    // FOR EPS COLORING SCHEME
    particle_coloring_scheme = -1;
    n_switch = 0;   // SWITCH USED FOR TURNING OFF DRAWING OF VORONOI CELLS

    // DEFAULT VALUE FOR PAIR CORRELATION FUNCTIONS
    const int DEFAULT_MAX_RADIUS = 10;

    // FIRST ARGUMENT (d=1) IS INPUT DATA FILENAME; PARSE REMAINING ARGUMENTS.
    for(int d=2; d<argc; d++)
    {
        if     (strcmp(argv[d],"-2")  ==0) dimension=2;    // SPECIFY FOR 2-DIMENSIONAL SYSTEMS
        else if(strcmp(argv[d],"-vt") ==0) vt_switch=1;    // PRINT OUT P-VECTORS OR W-VECTORS FOR EACH PARTICLE IN 2 AND 3 DIMENSIONS
        else if(strcmp(argv[d],"-c") ==0)   c_switch=1;    // CLUSTER ANALYSIS; TAKES OPTIONAL ARGUMENT
        else if(strcmp(argv[d],"-d")  ==0)  d_switch=1;    // PRINT OUT DISTRIBUTION OF VORONOI TOPOLOGIES
        else if(strcmp(argv[d],"-l")  ==0)  l_switch=1;    // OUTPUT LAMMPS DUMP FILE
        else if(strcmp(argv[d],"-r") ==0)   r_switch=1;    // RESOLVE INDETERMINATE TYPES

        // COMPUTE DISTRIBUTION FOR PERTURBATIONS OF A GIVEN SYSTEM
        else if(strcmp(argv[d],"-g")==0)
        {
            if(d+1 == argc || d+2 == argc)     // LAST OPTION
            {
                // WE HAVE TROUBLE SINCE THIS MEANS WE DON'T HAVE
                // TWO MORE VALUES, WHICH WOULD BE THE MIN AND
                // MAX VALUES WE ARE SEEKING
                help_message();
                std::cout << "-g option indicated; must be followed by two numbers, samples and perturbation.\n";
                exit(0);
            }
            else if(argv[d+1][0]=='-' || argv[d+2][0]=='-')
            {
                help_message();
                std::cout << "-g option indicated; must be followed by two numbers, samples and perturbation.\n";
                exit(0);
            }
            else  // WE READ IN THE TWO PARAMETERS FOLLOWING.  WE MAKE SURE THAT perturbation_samples IS AT 
                  // LEAST 1, AND MAKE SURE THAT perturbation_size >= 0.  IF THESE CONDITIONS ARE MET, THEN 
                  // WE USE THESE VALUES.  IF NOT, WE RETURN AN ERROR MESSAGE.                  
            {
                g_switch=1;
                
                perturbation_samples = atoi(argv[d+1]);
                perturbation_size    = atof(argv[d+2]);
                
                // WE CHECK HERE THAT perturbation_samples IS GREATER THAN 0 AND THAT perturbation_size IS AT LEAST 0
                if(0 < perturbation_samples && 0 <= perturbation_size)
                    d+=2;
                else
                {
                    help_message();
                    std::cout << "-g option indicated; must be followed by two numbers, samples and perturbation.\n";
                    exit(0);
                }
            }
        }
        
        // GENERATE FILTER FILE FROM CRYSTAL STRUCTURE
        else if(strcmp(argv[d],"-mf")==0)
        {
            mf_switch = 1;
            // CHECK IF TWO OPTIONAL NUMERIC ARGUMENTS FOLLOW (PERTURBATION MODE)
            if(d+2 < argc && argv[d+1][0]!='-' && argv[d+2][0]!='-')
            {
                perturbation_samples = atoi(argv[d+1]);
                perturbation_size    = atof(argv[d+2]);
                if(perturbation_samples < 1 || perturbation_size < 0)
                {
                    help_message();
                    std::cout << "-mf perturbation mode requires positive sample count and non-negative perturbation size.\n";
                    exit(0);
                }
                d += 2;
            }
        }

        // SPECIFY FILTER FILE
        else if(strcmp(argv[d],"-f")==0)                    
        {
            if(argc > d+1 && argv[d+1][0]!='-')
            {
                f_switch=1;
                filename_filter = argv[d+1];
                d++;
            }
            else
            {
                help_message();
                std::cout << "-f option indicated but no filter specified.\n";
                exit(0);
            }
        }

        // COMPUTE AND OUTPUT UNNORMALIZED VORONOI PAIR CORRELATION FUNCTION DATA
        // TAKES OPTIONAL ARGUMENT TO DETERMINE MAXIMAL RADIUS.
        else if(strcmp(argv[d],"-u")==0)                   
        {        
            // TAKES OPTIONAL ARGUMENT TO DETERMINE MAXIMAL RADIUS.
            u_switch = 1;
            if(argc > d+1 && argv[d+1][0]!='-')
            {
                max_radius = atoi(argv[d+1]);
                d++;
            }
            else                                           // DEFAULT VALUE max_radius = 10
                max_radius = DEFAULT_MAX_RADIUS;
        }

        // COMPUTE AND OUTPUT NORMALIZED VORONOI PAIR CORRELATION FUNCTION DATA
        // TAKES OPTIONAL ARGUMENT TO DETERMINE MAXIMAL RADIUS.
        else if(strcmp(argv[d],"-v")==0)                   
        {                                                  
            v_switch = 1;
            if(argc > d+1 && argv[d+1][0]!='-')
            {
                max_radius = atoi(argv[d+1]);
                d++;
            }
            else                                           // DEFAULT VALUE max_radius = 10
                max_radius = DEFAULT_MAX_RADIUS;
        }
        
        // OUTPUT EPS FILE; TAKES 0, 1, OR 2 ARGUMENTS.  IF NO ARGUMENTS FOLLOW, THEN DRAW ALL PARTICLES 
        // WITH DEFAULT COLORING ACCORDING TO NUMBER OF EDGES.  IF ONE ARGUMENT FOLLOWS, THEN DRAW ALL 
        // PARTICLES WITH THAT COLORING SCHEME. IF 2 ARGUMENTS FOLLOW, THEN THE FIRST CHOOSES THE COLOR 
        // SCHEME AND THE SECOND SPECIFIES THE ROUGH NUMBER OF PARTICLES DRAWN.
        else if(strcmp(argv[d],"-e")==0)                   
        {                                                 
            e_switch=1;                                    
                                                          
            if(argc > d+1 && argv[d+1][0]!='-')            
            {
                particle_coloring_scheme = atoi(argv[d+1]);  // FIRST OPTIONAL ARGUMENT GIVES COLOR SCHEME
                
                if(argc > d+2 && argv[d+2][0]!='-')        // SECOND ARGUMENT, DETERMINES NUMBER OF PARTICLES
                {
                    particles_in_eps = atoi(argv[d+2]);
                    d++;
                }                
                else particles_in_eps=0;                   // DRAW ALL PARTICLES

                d++;
            }
            else                                           // NO ARGUMENTS
            {
                if(particle_coloring_scheme==-1)           // IF COLORING SCHEME WAS SET BY -f SWITCH, THEN DON'T DEFAULT TO 2
                particle_coloring_scheme  = 2;             // DEFAULT particle_coloring_scheme (2) COLOR BY EDGE COUNT
                particles_in_eps = 0;                      // DRAW ALL PARTICLES
            }
            
            // CATCH ERROR IN COLOR SCHEME, OUTPUT INSTRUCTIONS
            if(particle_coloring_scheme<0 || particle_coloring_scheme>8)
            {
                help_message();
                std::cout << "Coloring scheme is given as " << particle_coloring_scheme << ", but must be between 0 and 8" << '\n';
                exit(0);
            }
        }

        // DO NOT DRAW VORONOI CELLS, ONLY PARTICLES
        else if(strcmp(argv[d],"-n") ==0)  n_switch=1;     

        // SPECIFIC NUMBER OF THREADS FOR MULTITHREADED VERSION
        else if(strcmp(argv[d],"-t") ==0)                  
        {
            if(argc > d+1 && argv[d+1][0]!='-')
                threads = atoi(argv[d+1]);
            d++;
        }
 
        // OUTPUT HELP MESSAGE
        else if(strcmp(argv[d],"-h")==0 || strcmp(argv[d],"-help")==0) 
        {
            help_message();
            exit(0);
        }

        // UNIDENTIFIED OPTION
        else
        {
            help_message();
            std::cout << "Unidentified option '" + std::string(argv[d]) + "' included.\n";
            exit(0);
        }
    }
        

    // CHECK COMPATIBILITY OF OPTIONS
    if(r_switch && !f_switch)
    {
        std::cout << "-r option requires filter file for resolution analysis\n";
        exit(0);
    }
    
    if(c_switch && !f_switch)
    {
        std::cout << "-c option requires filter file for cluster analysis\n";
        exit(0);
    }

    if(e_switch && l_switch)
    {
        std::cout << "-e and -l options not compatible.\n";
        exit(0);
    }

    if(e_switch && dimension==3)
    {
        std::cout << "-e option not compatible with 3D systems; use -2 option to specify treating system as 2D.\n";
        exit(0);
    }

    if(e_switch && !f_switch && (particle_coloring_scheme==3 || particle_coloring_scheme==5 || particle_coloring_scheme==6 || particle_coloring_scheme==7 || particle_coloring_scheme==8))
    {
        std::cout << "EPS coloring schemes 3, 5, 6, 7, and 8 require specifying a filter.\n";
        std::cout << "Changing coloring scheme to 2.\n";
        particle_coloring_scheme = 2; // DEFAULT COLORING SCHEME
    }

    // CLUSTER ANALYSIS REQUIRED FOR COLORING SCHEME 5, 6, 7, and8
    if(e_switch && (particle_coloring_scheme==5 || particle_coloring_scheme==6 || particle_coloring_scheme==7 || particle_coloring_scheme==8))
        c_switch = 1; 
    
}




