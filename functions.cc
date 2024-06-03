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
    "         VoroTop reads a LAMMPS dump input file and computes the complete Voronoi\n"
    "         cell topology for each point.  Output is determined by options.         \n"
    "                                                                                 \n"
    "                                                                                 \n"
    "Available options:                                                               \n"
    "                                                                                 \n"
    " -vt   : Compute topology of each particle, save to <filename>.vectors.          \n"
    " -d    : Compute distribution of Voronoi topologies, represented as either       \n"
    "         p-vectors or Weinberg vectors, and save to <filename>.distribution.     \n"
    " -f    : Generate a LAMMPS <filename>.dump file using a filter provided by user  \n"
    "                                                                                 \n"
    " -t    : Specify number of threads to use.                                       \n"
    " -g    : Compute distribution of w-vectors for perturbations of a system. This   \n"
    "         option should be followed by two numbers, one specifying the number of  \n"
    "         samples, the other specifying the size of the random perturbations.     \n"
    "                                                                                 \n"
    " -r    : Resolve indeterminate types.                                            \n"
    " -c    : Cluster analysis.  Particles with topological types not included in the \n"
    "         specified filter are considered to be defects.  Defect particles with   \n"
    "         adjacent Voronoi cells are considered to belong to the same defect; we  \n"
    "         call each Voronoi-contiguous set of defect particles a cluster and      \n"
    "         assign the constituent particles a unique index.  Similar cluster       \n"
    "         analysis is also performed for types in the filter.  This option        \n"
    "         requires loading a filter file for analysis.                            \n"
    "                                                                                 \n";
}



////////////////////////////////////////////////////
////
////   PARSE INPUT ARGUMENTS AND DETERMINE WHAT
////   SHOULD BE DONE
////
////////////////////////////////////////////////////

void parse_arguments(int argc, char *argv[])
{
    if(argc<2)
    {
        help_message();
        std::cerr << "No input data file specified.\n\n\n";
        exit(-1);
    }
    
    if(argc<3)
    {
        help_message();
        std::cerr << "No command-line options specified.\n\n\n";
        exit(-1);
    }
    
    // DETERMINE INPUT AND OUTPUT FILENAMES, FIX, SINCE WE DON'T NEED ALL OF THIS PARSING AND CONTROL
    filename_data = argv[1];
    //name_of_data_file = argv[1];
    //std::size_t found = name_of_data_file.find_last_of("/\\");
    //std::string path  = name_of_data_file.substr(0,found+1);
    //std::string file  = name_of_data_file.substr(found+1);
    
    
    // DEFAULT THREADS IF NOT SPECIFIED BY USER
#ifdef _OPENMP
    if(omp_get_max_threads()>2) threads = omp_get_max_threads()-2;
    else                        threads = 1;
#else
    threads = 1;
#endif
    

    // DEFAULT DIMENSION: 3
    dimension = 3;

    // FOR EPS COLORING SCHEME
    particle_coloring_scheme = -1;
    
    // FIRST ARGUMENT (d=1) IS INPUT DATA FILENAME; PARSE REMAINING ARGUMENTS.
    for(int d=2; d<argc; d++)
    {
        if     (strcmp(argv[d],"-2")  ==0) dimension=2;    // PRINT OUT W-VECTORS FOR EACH PARTICLE
        else if(strcmp(argv[d],"-vt") ==0) vt_switch=1;    // PRINT OUT W-VECTORS FOR EACH PARTICLE
        else if(strcmp(argv[d],"-d")  ==0)  d_switch=1;    // PRINT OUT DISTRIBUTION OF W-VECTORS OF SYSTEM

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
            else  // AT THIS POINT WE HAVE TWO OBJECTS FOLLOWING.  WE DETERMINE THESE VALUES
                  // AND MAKE SURE THAT THEY ARE VALUES BETWEEN 0 AND 1; IF THEY ARE WE USE
                  // THEM, OTHERWISE WE RETURN AN ERROR MESSAGE.
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
        
        else if(strcmp(argv[d],"-f")==0)                        // SPECIFY FILTER FILE
        {
            if(argc > d+1 && argv[d+1][0]!='-')
            {
                f_switch=1;
                filename_filter = argv[d+1];
                d++;
                particle_coloring_scheme = 3;
            }
            else
            {
                help_message();
                std::cout << "-f option indicated but no filter specified.\n";
                exit(0);
            }
        }

        else if(strcmp(argv[d],"-u")==0)                   // OUTPUT UNNORMALIZED VORONOI PAIR CORRELATION FUNCTION DATA
        {        
            // TAKES OPTIONAL ARGUMENT TO DETERMINE MAXIMAL RADIUS.
            u_switch = 1;
            if(argc > d+1 && argv[d+1][0]!='-')
            {
                max_radius = atoi(argv[d+1]);
                d++;
            }
            else                                           // DEFAULT VALUE max_radius = 20
                max_radius = 20;
        }

        else if(strcmp(argv[d],"-v")==0)                   // OUTPUT NORMALIZED VORONOI PAIR CORRELATION FUNCTION DATA
        {                                                  // TAKES OPTIONAL ARGUMENT TO DETERMINE MAXIMAL RADIUS.
            v_switch = 1;
            if(argc > d+1 && argv[d+1][0]!='-')
            {
                max_radius = atoi(argv[d+1]);
                d++;
            }
            else                                           // DEFAULT VALUE max_radius = 20
                max_radius = 20;
        }
        
        else if(strcmp(argv[d],"-r") ==0)                // RESOLVE INDETERMINATE TYPES
        {
            r_switch=1;
            if(argc > d+1 && argv[d+1][0]!='-')
            {
                resolve_trials = atoi(argv[d+1]);
                if(0 <= resolve_trials)
                    d++;
                else
                {
                    help_message();
                    std::cout << "-r option indicated; can be followed by optional integer for number of trials.\n";
                    exit(0);
                }
            }
            else
                resolve_trials = 5;
        }
        
        else if(strcmp(argv[d],"-c") ==0)                  // CLUSTER ANALYSIS; TAKES OPTIONAL ARGUMENT
        {
            c_switch=1;
            if(argc > d+1 && argv[d+1][0]!='-')
            {
                clustering_default = atoi(argv[d+1]);
                if(0 <= clustering_default)
                    d++;
                else
                {
                    help_message();
                    std::cout << "-c option indicated; can be followed by optional index for clustering.\n";
                    exit(0);
                }
                clustering_default_switch = 1;
            }
            else
            {
                clustering_default        = 0;
                clustering_default_switch = 0;
            }
        }
       
        // THIS SHOULD TAKE 0, 1, OR 2 ARGUMENTS.  IF 0, THEN DEFFAULT COLORING.  IF 1, THEN COLOR SCHEME
        // AND DRAW ENTIRE SYSTEM.  IF 2, THEN COLOR SCHEME AND DRAW THAT NUMBER OF PARTICLES.
        
        else if(strcmp(argv[d],"-e")==0)                   // OUTPUT EPS FILE; TAKES 0, 1, OR 2 ARGUMENTS.
        {                                                  // 0 ARGUMENTS: DEFAULT COLOR SCHEME, DRAW ENTIRE SYSTEM
            e_switch=1;                                    // 1 ARGUMENT : COLOR SCHEME, DRAW ENTIRE SYSTEM
            if(argc > d+1 && argv[d+1][0]!='-')            // 2 ARGUMENTS: COLOR SCHEME, DRAW GIVEN NUMBER OF PARTICLES
            {
                particle_coloring_scheme = atoi(argv[d+1]);         // FIRST OPTIONAL ARGUMENT GIVES COLOR SCHEME
                
                if(argc > d+2 && argv[d+2][0]!='-')        // SECOND ARGUMENT, DETERMINES NUMBER OF PARTICLES
                {
                    particles_in_eps = atoi(argv[d+2]);
                    d++;
                }
                
                d++;
            }
            else                                           // NO ARGUMENTS
            {
                if(particle_coloring_scheme==-1)           // IF COLORING SCHEME WAS SET BY -f SWITCH, THEN DON'T DEFAULT TO 2
                particle_coloring_scheme  = 2;             // DEFAULT particle_coloring_scheme (2) COLOR BY EDGE COUNT
                particles_in_eps = 0;                      // DRAW ALL PARTICLES
            }
            
            // CATCH ERROR IN COLOR SCHEME, OUTPUT INSTRUCTIONS
            if(particle_coloring_scheme<0 || particle_coloring_scheme>4)
            {
                help_message();
                std::cout << "Coloring scheme is given as " << particle_coloring_scheme << ", but must be between 0 and 4" << '\n';
                exit(0);
            }
        }
        else if(strcmp(argv[d],"-n") ==0)  n_switch=1;     // DO NOT DRAW VORONOI CELLS; ONLY DRAW PARTICLES

        else if(strcmp(argv[d],"-t") ==0)                  // SPECIFIC NUMBER OF THREADS FOR MULTITHREADED VERSION
        {
            if(argc > d+1 && argv[d+1][0]!='-')
                threads = atoi(argv[d+1]);
            d++;
        }
 
        else
        {
            help_message();
            std::cout << "Unidentified option \'" << argv[d] << "\' included.\n\n";
            exit(0);
        }
    }
        
    
    // CHECK COMBINATIONS OF OPTIONS
    if(r_switch && !f_switch)
    {
        std::cout << "-r option requires filter file for resolution analysis\n\n";
        exit(0);
    }
    
    if(c_switch && !f_switch)
    {
        std::cout << "-c option requires filter file for cluster analysis\n\n";
        exit(0);
    }
    
    if(vt_switch && f_switch)
    {
        std::cout << "-w option not compatible with other options\n\n";
        exit(0);
    }
    
    //filename_output = path + file;
}




