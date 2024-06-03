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

////   File: vorotop.cc


#include <random>
#include <ctime>
#include <cstdio>
#include <vector>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <chrono>

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
    
    data_file.open(filename_data, std::ifstream::in);
    if(!data_file.is_open())
    {
        help_message();
        std::cerr << "Unable to open input file " << filename_data << "\n\n";
        exit(-1);
    }
    
    
    ////////////////////////////////////////////////////
    ////
    ////   PARSE DATA FILE TO DETERMINE FILE TYPE
    ////
    ////////////////////////////////////////////////////
    
    parse_header(data_file);
    
    
    ////////////////////////////////////////////////////
    ////
    ////   CREATE SPACE FOR DATA
    ////
    ////////////////////////////////////////////////////
    
    vt_structure_types.resize          (number_of_particles);  // MEMORY FOR STRUCTURE TYPES
    vt_structure_types_resolved.resize (number_of_particles);  // MEMORY FOR RESOLVED STRUCTURE TYPES
    if(dimension==2 || u_switch || v_switch || c_switch)
    {
        neighbors_list_char.resize         (number_of_particles);  // MEMORY FOR LIST OF NEIGHBORS
        cell_neighbor_count.resize         (number_of_particles);
    }
    
    if(c_switch)                                        // MEMORY FOR CLUSTER ANALYSIS
    {
        cluster_index.resize   (number_of_particles);
        cluster_sizes.resize   (number_of_particles);
        std::fill(cluster_sizes.begin(),cluster_sizes.end(), 0);
    }
    
    if(particle_coloring_scheme==4)
        ring_index.resize(number_of_particles);
    
    particle_coordinates=new double[number_of_particles*dimension];
    particle_ids        =new    int[number_of_particles];
    
    
    if(r_switch)                                        // MEMORY FOR RESOLVED TYPES
    {
        vt_structure_types.resize(number_of_particles);
        volumes.resize           (number_of_particles);
    }
    
    ////////////////////////////////////////////////////
    ////
    ////   VARIABLES AND HELP MESSAGES IF REQUESTED
    ////
    ////////////////////////////////////////////////////
    
    if(d_switch==0 && g_switch==0 && vt_switch==0 && c_switch==0 &&
       f_switch==0 &&  u_switch==0 && v_switch==0 && e_switch==0)
    {
        help_message();
        exit(-1);
    }
    
    data_file.close();
    
    
    ////////////////////////////////////////////////////
    ////
    ////   CREATE VORO++ CONTAINER WITH NECESSARY DIMENSIONS;
    ////   THEN ADD PARTICLES TO THE CONTAINER
    ////
    ////////////////////////////////////////////////////
    
    // IMPORTS ALL COORDINATE DATA
    import_data(data_file);
    
    // CREATES CONTAINERS
    voro::container_2d con2d (origin[0],hi_bound[0],origin[1],hi_bound[1],n_x,n_y,true,true,4,threads);
    voro::container_3d con3d (supercell_edges[0][0],supercell_edges[1][0],supercell_edges[1][1],
                              supercell_edges[2][0],supercell_edges[2][1],supercell_edges[2][2],
                              n_x,n_y,n_z,true,true,true,8,threads);
    
    if(dimension==2)
    {
        con2d.add_parallel(particle_coordinates, number_of_particles, threads);
        con2d.put_reconcile_overflow();
    }

    if(dimension==3)
    {
        con3d.add_parallel(particle_coordinates, number_of_particles, threads);
        con3d.put_reconcile_overflow();
    }
    

    ////////////////////////////////////////////////////
    ////
    ////   LOAD FILTER IF REQUESTED
    ////
    ////////////////////////////////////////////////////
    
    Filter filter;
    if(f_switch) filter.load_filter(filename_filter);
    
    
    ////////////////////////////////////////////////////
    ////
    ////   PROCESS DATA ACCORDING TO INPUT OPTIONS
    ////
    ////////////////////////////////////////////////////
    
    // OUTPUTS VORONOI TOPOLOGY FOR EACH PARTICLE
    if(vt_switch)
    {
        if(dimension==2)
        {
            count_and_store_neighbors(con2d);
            print_topology_vectors2d(filename_data);
        }
        
        if(dimension==3)
        {
            print_topology_vectors(con3d,filename_data);
        }
    }
    
    // OUTPUTS DISTRIBUTION OF VORONOI TOPOLOGIES
    else if(d_switch)
    {
        if(dimension==2)
        {
            count_and_store_neighbors(con2d);
            calc_distribution(con2d, filter);
            filter.print_distribution(filename_data);
        }
        
        if(dimension==3)
        {
            calc_distribution(con3d, filter);
            filter.print_distribution(filename_data);
        }
    }
    
    
    // OUTPUTS DISTRIBUTION OF VORONOI TOPOLOGIES
    else if(g_switch)
    {
        if(dimension==2)
        {
            calc_gaussian_distribution(filter);
            filter.print_distribution(filename_data);
        }
        if(dimension==3)
        {
            calc_gaussian_distribution3d(filter);
            filter.print_distribution(filename_data);
        }
    }
    

    // OUTPUTS LAMMPS DUMP FILE, INCLUDING CLASSIFICATION USING FILTER
    else if(f_switch)
    {
        if(dimension==2)
        {
            count_and_store_neighbors(con2d);
            classify_particles_by_voronoi_topology(filter);
            output_system (filename_data,filter);
        }
        if(dimension==3)
        {
            classify_particles_by_voronoi_topology(con3d, filter);
            output_system (filename_data,filter);
        }
    }
    
    else if(u_switch || v_switch)  
    {
        if(dimension==2) count_and_store_neighbors(con2d);
        if(dimension==3) count_and_store_neighbors(con3d);
        pair_correlation_analysis ();
    }

    else if(e_switch && dimension==2)
    {
        count_and_store_neighbors(con2d);
        // COLORING SCHEMES 0,1,2,3 REQUIRE NO ANALYSIS;
        // COLORING SCHEME 4 REQUIRES LABELING PARTICLES ACCORDING TO RING NUMBER
        if(particle_coloring_scheme == 4) ring_coloring();
        output_eps(con2d,filename_data);
    }

    else
    {
        if(r_switch)             resolve_indeterminate_types(con3d,filter);
        if(c_switch)             cluster_analysis                 (filter);
    }
    
    delete particle_coordinates;
    delete particle_ids;
    
    return 0;
}




