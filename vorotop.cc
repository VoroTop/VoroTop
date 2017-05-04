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

////   File: vorotop.cc


#include <vector>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "filters.hh"
#include "variables.hh"
#include "functions.hh"



int main(int argc, char *argv[])
{
    ////////////////////////////////////////////////////
    ////
    ////   PARSE ALL COMMAND-LINE OPTIONS.  DETERMINE
    ////   FILE NAMES AND DIRECTORIES TO USE, AND
    ////   WHICH ANALYSES TO PERFORM.
    ////
    ////////////////////////////////////////////////////
    
    parse_arguments(argc, argv);
    
    
    ////////////////////////////////////////////////////
    ////
    ////   OPEN INPUT FILE
    ////
    ////////////////////////////////////////////////////
    
    std::ifstream data_file;
    data_file.open(name_of_data_file);
    if(!data_file.is_open())
    {
        help_message();
        std::cerr << "Unable to open input file " << name_of_data_file << "\n\n";
        exit(-1);
    }
    
    
    ////////////////////////////////////////////////////
    ////
    ////   PARSE DATA FILE TO DETERMINE FILE TYPE
    ////
    ////////////////////////////////////////////////////
    
    int file_type = parse_header(data_file);
    

    ////////////////////////////////////////////////////
    ////
    ////   CREATE SPACE FOR DATA
    ////
    ////////////////////////////////////////////////////
    
    int entries = attribute_labels.size();
    if(file_type==2 && no_velocity==0) entries+=3;
    
    // ALLOCATE MEMORY FOR INPUT PARTICLE DATA
    particle_data.resize(number_of_particles, std::vector<double>(entries));  // DATA FOR OUTPUT
    all_wvectors.resize (number_of_particles);
    xcoord.resize (number_of_particles);
    ycoord.resize (number_of_particles);
    zcoord.resize (number_of_particles);

    // ALLOCATE MEMORY FOR CLUSTER ANALYSIS
    if(c_switch)
    {
        neighbors_list.resize  (number_of_particles);
        neighbors_list_c.resize(number_of_particles);
        cluster_index.resize   (number_of_particles);
        cluster_sizes.resize   (number_of_particles);
        
        std::fill(cluster_sizes.begin(),cluster_sizes.end(), 0);
    }
    

    ////////////////////////////////////////////////////
    ////
    ////   VARIABLES AND HELP MESSAGES IF REQUESTED
    ////
    ////////////////////////////////////////////////////
    
    if(p_switch==1) print_variables();
    if(d_switch==0 && g_switch==0 && df_switch==0 &&
       w_switch==0 && c_switch==0 && f_switch==0)
    {
        help_message();
        exit(-1);
    }
    

    ////////////////////////////////////////////////////
    ////
    ////   CREATE VORO++ CONTAINER WITH NECESSARY DIMENSIONS;
    ////   THEN ADD PARTICLES TO THE CONTAINER
    ////
    ////////////////////////////////////////////////////
    
    voro::particle_order vo;
    voro::container_periodic con (supercell_edges[0][0],supercell_edges[1][0],supercell_edges[1][1],
                                  supercell_edges[2][0],supercell_edges[2][1],supercell_edges[2][2], n_x,n_y,n_z,8);
    
    if     (file_type==1) import_dump_file   (data_file, con, vo);
    else if(file_type==2) import_atomeye_file(data_file, con, vo);
    else
    {
        std::cout << "File format not recognized\n";
        exit(-1);
    }

    
    ////////////////////////////////////////////////////
    ////
    ////   LOAD FILTER IF REQUESTED
    ////
    ////////////////////////////////////////////////////

    Filter filter;
    if(f_switch) filter.loadFilter(filename_filter);

    
    ////////////////////////////////////////////////////
    ////
    ////   PROCESS DATA ACCORDING TO INPUT OPTIONS
    ////
    ////////////////////////////////////////////////////

    if(w_switch)
    {
        calc_all_wvectors(con,vo,1);
        print_wvectors(name_of_data_file);
    }
    
    else
    {
        calc_all_wvectors(con,vo,0);

        if(d_switch || df_switch)             calc_distribution(filter);
        if(g_switch)                          calc_gaussian_distribution(con,vo,filter);
        if(d_switch || df_switch || g_switch) filter.print_distribution(name_of_data_file);
        if(c_switch)                          cluster_analysis(filter);
        if(f_switch || df_switch)             create_cfg_file(filename_output,filter);
    }
    
    return 0;
}









