////////////////////////////////////////////////////////
////                                                ////
////   ******************************************   ////
////   *                                        *   ////
////   *     VoroTop: Voronoi Cell Topology     *   ////
////   *   Visualization and Analysis Toolkit   *   ////
////   *             (Version 0.4)              *   ////
////   *                                        *   ////
////   *           Emanuel A. Lazar             *   ////
////   *      University of Pennsylvania        *   ////
////   *           December 5, 2014             *   ////
////   *                                        *   ////
////   ******************************************   ////
////                                                ////
////////////////////////////////////////////////////////

////   File: functions.cc


#include <cstring>
#include <iostream>
#include <sys/stat.h>

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
    "   *              (Version 0.4)               *\n"
    "   *                                          *\n"
    "   *            Emanuel A. Lazar              *\n"
    "   *       University of Pennsylvania         *\n"
    "   *            December 5, 2014              *\n"
    "   *                                          *\n"
    "   ********************************************\n"
    "                                               \n"
    "Syntax: VoroTop <filename> [options]           \n"
    "                                               \n"
    "         VoroTop reads an input file and computes the complete Voronoi cell      \n"
    "         topology for each point.  Input formats currently accepted include      \n"
    "         LAMMPS dump and AtomEye cfg formats.  Output is determined by options.  \n"
    "                                                                                 \n"
    "                                                                                 \n"
    "Available options:                                                               \n"
    "                                                                                 \n"
    " -f    : Generate a LAMMPS <filename>.dump or AtomEye <filename>.cfg file using  \n"
    "         a filter provided by the user.                                          \n"
    "                                                                                 \n"
    " -w    : Compute w-vector for each particle, save to file <filename>.wvectors.   \n"
    " -d    : Compute distribution of w-vectors, save to file <filename>.distribution.\n"
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
    "                                                                                 \n"
    " -o    : Specify output file primary name <filename>.  If requested, w-vectors   \n"
    "         will be saved to <filename>.wvectors; distribution will be saved to     \n"
    "         <filename>.distribution; AtomEye data will be saved to <filename>.cfg.  \n"
    " -od   : Specify output directory.  If requested, w-vectors or AtomEye cfg files \n"
    "         will be saved to specified directory, using current filename; directory \n"
    "         must previously exist.                                                  \n"
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
    
    // DETERMINE INPUT AND OUTPUT FILENAMES
    name_of_data_file = argv[1];
    std::size_t found = name_of_data_file.find_last_of("/\\");
    std::string path  = name_of_data_file.substr(0,found+1);
    std::string file  = name_of_data_file.substr(found+1);
    
    // FIRST ARGUMENT (d=1) IS INPUT DATA FILENAME; PARSE REMAINING ARGUMENTS.
    for(int d=2; d<argc; d++)
    {
        if(strcmp(argv[d],"-f")==0)                        // SPECIFY FILTER FILE
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

        else if(strcmp(argv[d],"-w") ==0)  w_switch=1;    // PRINT OUT W-VECTORS FOR EACH PARTICLE
        else if(strcmp(argv[d],"-d") ==0)  d_switch=1;     // PRINT OUT DISTRIBUTION OF W-VECTORS OF SYSTEM

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
        
        // BY DEFAULT, INPUT FILE IS USED AS A BASENAME FOR OUTPUT FILES.  NEW
        // BASENAME OR OUTPUT DIRECTORIES CAN BE SPECIFIED.
        else if(strcmp(argv[d],"-od")==0)                       // SPECIFY OUTPUT DIRECTORY
        {
            if(argc > d+1 && argv[d+1][0]!='-')
            {
                path = argv[d+1];                               // WE HAVE CHANGED THE PATH
                if(path.back()!='\\' && path.back()!='/')
                    path.append("/");
                d++;
                
                // IF SPECIFIED OUTPUT DIRECTORY DOES NOT EXIST, CREATE IT
                struct stat st;
                if (stat(path.c_str(), &st) == -1)
                    mkdir(path.c_str(), 0744);
            }
            else
            {
                help_message();
                std::cout << "-od option indicated but no output directory specified.\n";
                exit(0);
            }
        }
        
        else if(strcmp(argv[d],"-o")==0)                        // SPECIFY OUTPUT FILE NAME
        {
            if(argc > d+1 && argv[d+1][0]!='-')
            {
                file = argv[d+1];                               // WE HAVE CHANGED THE FILE NAME
                d++;
            }
            else
            {
                help_message();
                std::cout << "-o option indicated but no output file specified.\n\n";
                exit(0);
            }
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
    
    if(w_switch && f_switch)
    {
        std::cout << "-w option not compatible with other options\n\n";
        exit(0);
    }
    
    filename_output = path + file;
}





