////////////////////////////////////////////////////////
////                                                ////
////   ******************************************   ////
////   *                                        *   ////
////   *     VoroTop: Voronoi Cell Topology     *   ////
////   *   Visualization and Analysis Toolkit   *   ////
////   *             (Version 0.3)              *   ////
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
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <sys/stat.h>

#include "filters.hh"
#include "vorotop.hh"



////////////////////////////////////////////////////
////
////   PRINT SYSTEM PROPERTIES TO STANDARD OUTPUT
////
////////////////////////////////////////////////////

void print_variables()
{
    std::cout << "Global variables\n";
    std::cout << "  Number of particles: " << number_of_particles << "\n";
    std::cout << "  Origin: (" << origin[0] << ", " << origin[1] << ", " << origin[2] << ")\n";
    std::cout << "  supercell edge A " << supercell_edges[0][0] << '\t' << supercell_edges[0][1] <<  '\t' << supercell_edges[0][2] << '\n';
    std::cout << "  supercell edge B " << supercell_edges[1][0] << '\t' << supercell_edges[1][1] <<  '\t' << supercell_edges[1][2] << '\n';
    std::cout << "  supercell edge C " << supercell_edges[2][0] << '\t' << supercell_edges[2][1] <<  '\t' << supercell_edges[2][2] << '\n';
    std::cout << "\n";
    
    std::cout << "LAMMPS variables\n";
    std::cout << "  timestep:    " << timestep << "\n";
    std::cout << "  entry count: " << attribute_labels.size() << "\n";
    for(int i=0; i<attribute_labels.size(); i++)
        std::cout << "  entry " << i << ":     " << attribute_labels[i] << "\n";
    std::cout << "scaled coords. " << scaled_coordinates << '\n';
    std::cout << "\n";
    
    std::cout << "CFG variables\n";
    std::cout << "  scale:       " << cfg_lscale << "\n";
    std::cout << "  atom mass:   " << cfg_atomic_mass << "\n";
    std::cout << "  chem symbol: " << cfg_chem_symbol << "\n";
    std::cout << "  entry count: " << attribute_labels.size() << "\n";
    for(int i=0; i<attribute_labels.size(); i++)
        std::cout << "  entry " << i << ":     " << attribute_labels[i] << "\n";
    std::cout << "\n";
}



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
         "   *              (Version 0.3)               *\n"
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
         " -w    : Compute w-vector for each particle, save to file <filename>.wvectors.   \n"
         " -d    : Compute distribution of w-vectors, save to file <filename>.distribution.\n"
         " -g    : Compute distribution of w-vectors for perturbations of a system. This   \n"
         "         option should be followed by two numbers, one specifying the number of  \n"
         "         samples, the other specifying the size of the random perturbations.     \n"
         "                                                                                 \n"
         " -ff   : Generates a cfg file <filename>.cfg using a filter provided by the user.\n"
         " -df   : Generates a cfg file <filename>.cfg using a distribution of types in    \n"
         "         the sample; distribution is saved to <filename>.distribution.  Each     \n"
         "         topological type is given a unique structure ID; structure ID's of more \n"
         "         common types have larger numbers than those of less common types.  This \n"
         "         can be useful for studying new structures.                              \n"
         "                                                                                 \n"
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
         "                                                                                 \n"
         " -p    : Print to screen system variables.                                       \n"
         " -h    : Print this information.                                                 \n"
         "                                                                                 \n";
}



////////////////////////////////////////////////////
////
////   COUNT WORDS IN A STRING
////
////////////////////////////////////////////////////

unsigned int countWordsInString(std::string const& str)
{
    std::stringstream stream(str);
    return std::distance(std::istream_iterator<std::string>(stream), std::istream_iterator<std::string>());
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
        if     (strcmp(argv[d],"-p") ==0)  p_switch=1;     // PRINT VARIABLES
        else if(strcmp(argv[d],"-h") ==0)  h_switch=1;     // PRINT HELP MENU
        
        else if(strcmp(argv[d],"-w") ==0)  w_switch=1;     // PRINT OUT W-VECTORS FOR EACH PARTICLE
        else if(strcmp(argv[d],"-d") ==0)  d_switch=1;     // PRINT OUT DISTRIBUTION OF W-VECTORS OF SYSTEM
        else if(strcmp(argv[d],"-df")==0) df_switch=1;     // DISTRIBUTION FILTER
        else if(strcmp(argv[d],"-c") ==0)  c_switch=1;     // CLUSTERING
        
        
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
                // THEM, OTHERWISE WE SEND AN ERROR MESSAGE.
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
        
        else if(strcmp(argv[d],"-ff")==0)                       // SPECIFY FILTER FILE
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
                std::cout << "-ff option indicated but no filter specified.\n";
                exit(0);
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
                file = argv[d+1];                                // WE HAVE CHANGED THE FILE NAME
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
    
    if(c_switch && !f_switch)
    {
        std::cout << "-c option requires filter file for cluster analysis\n\n";
        exit(0);
    }
    
    if(w_switch && (d_switch || df_switch || f_switch))
    {
        std::cout << "-w option not compatible with other options\n\n";
        exit(0);
    }
    
    filename_output = path + file;
}



////////////////////////////////////////////////////
////
////   COMPUTES W-VECTORS FOR ALL PARTICLES;
////   REPORTS NUMBER THAT ARE "BROKEN".
////
////////////////////////////////////////////////////

int calc_all_wvectors(voro::container_periodic& con, voro::particle_order& vo, bool extended)
{
    int broken = 0;     // NOTIFY USER IF SOME VORONOI CELLS DON'T COMPUTE
    
    // IF NEIGHBOR INFORMATION IS NECESSARY, FOR CLUSTERING
    if(c_switch)
    {
        voro::c_loop_order_periodic vlo(con,vo);
        voro::voronoicell_neighbor vcell_neighbors;
        
        if(vlo.start()) do
        {
            if(con.compute_cell(vcell_neighbors,vlo))
            {
                int pid = vlo.id[vlo.ijk][vlo.q];
                vcell_neighbors.neighbors(neighbors_list[pid]);
                neighbors_list_c[pid] = neighbors_list[pid].size();
                
                all_wvectors[pid] = calc_wvector(vcell_neighbors,extended);
            }
            else broken++;
        }
        while(vlo.inc());
    }
    
    else
    {
        voro::c_loop_order_periodic vlo(con,vo);
        voro::voronoicell vcell;

        if(vlo.start()) do
        {
            if(con.compute_cell(vcell,vlo))
            {
                int pid = vlo.id[vlo.ijk][vlo.q];
                all_wvectors[pid] = calc_wvector(vcell,extended);
            }
            else broken++;
        }
        while(vlo.inc());
    }
    
    if(broken>0)
        std::cout << "Voronoi cells of " << broken << " particles did not compute\n";
    
    return broken;      // 0 IF ALL CELLS COMPUTED, POSITIVE OTHERWISE
}



////////////////////////////////////////////////////
////
////   ADDS COMPUTED WVECTORS TO FILTER MAP
////
////////////////////////////////////////////////////

void calc_distribution(Filter &filter)
{
    // WVECTORS ARE ALREADY COMPUTED.  COPY THEM, SORT THEM, COMBINE
    // THEM, AND ADD THEM TO THE FILTER.
    std::vector< std::vector <int> > data_wvectors = all_wvectors;
    sort(data_wvectors.begin(), data_wvectors.end());

    int last = 0;
    int counter = 1;
    for(int c=1; c<number_of_particles; c++)
    {
        if(data_wvectors[c]!=data_wvectors[last])
        {
            filter.increment_or_add(data_wvectors[last],counter);
            last = c;
            counter=1;
        }
        else counter++;
    }
    filter.increment_or_add(data_wvectors[last],counter);
}







