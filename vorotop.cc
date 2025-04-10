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


#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "filters.hh"
#include "variables.hh"
#include "functions.hh"


void cleanup() {
    if (particle_coordinates) {
        delete[] particle_coordinates;
        particle_coordinates = nullptr;
    }
    if (particle_ids) {
        delete[] particle_ids;
        particle_ids = nullptr;
    }
}


void handle_error(const std::string& message) {
    cleanup();
    throw std::runtime_error(message);
}


int main(int argc, char *argv[])
{
    ////////////////////////////////////////////////////
    ////
    ////   PARSE ALL COMMAND-LINE OPTIONS.  DETERMINE
    ////   FILE NAMES AND DIRECTORIES TO USE, AND
    ////   WHICH ANALYSES TO PERFORM.
    ////
    ////////////////////////////////////////////////////
    
    parse_command_line_options(argc, argv);
    
    
    ////////////////////////////////////////////////////
    ////
    ////   OPEN INPUT FILE
    ////
    ////////////////////////////////////////////////////
    
    data_file.open(filename_data, std::ifstream::in);
    if(!data_file.is_open())
    {
        help_message();
        handle_error("Unable to open input file " + filename_data);
    }
    
    
    ////////////////////////////////////////////////////
    ////
    ////   PARSE DATA FILE TO DETERMINE FILE TYPE
    ////
    ////////////////////////////////////////////////////
    
    parse_header(data_file);
    data_file.close();


    ////////////////////////////////////////////////////
    ////
    ////   CREATE SPACE FOR DATA
    ////
    ////////////////////////////////////////////////////
    
    particle_ids        =new    int[number_of_particles];
    if (!particle_ids) {
        handle_error("Memory allocation for particle_ids failed.");
    }
    
    particle_coordinates=new double[number_of_particles*dimension];
    if (!particle_coordinates) {
        handle_error("Memory allocation for particle_coordinates failed.");
    }

    vt_structure_types.resize          (number_of_particles);  // MEMORY FOR STRUCTURE TYPES

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
    
    if(particle_coloring_scheme==4)                     // PARTICLE DISTANCE FROM CENTRAL PARTICLE
        ring_index.resize(number_of_particles);         // USED FOR COLORING PARTICLES IN EPS
        
    if(r_switch)                                        // MEMORY FOR RESOLVED TYPES
        vt_structure_types_resolved.resize(number_of_particles);
    
    ////////////////////////////////////////////////////
    ////
    ////   VARIABLES AND HELP MESSAGES IF REQUESTED
    ////
    ////////////////////////////////////////////////////
    
    if(d_switch==0 && g_switch==0 && vt_switch==0 && c_switch==0 &&
       f_switch==0 &&  u_switch==0 && v_switch==0 && e_switch==0)
    {
        help_message();
        handle_error("No valid options specified.");
    }
        
    
    // IMPORTS ALL COORDINATE DATA
    import_data();
    
    ////////////////////////////////////////////////////////
    ////  CREATE VORO++ CONTAINER WITH NECESSARY DIMENSIONS;
    ////  THEN ADD PARTICLES TO THE CONTAINER
    ////////////////////////////////////////////////////////
    
    // Declare pointers for the containers
    voro::container_2d* con2d = nullptr;
    voro::container_3d* con3d = nullptr;

    // Dynamically allocate the appropriate container based on the dimension
    if (dimension == 2) {
        con2d = new voro::container_2d(origin[0], hi_bound[0], origin[1], hi_bound[1], n_x, n_y, true, true, 4, threads);
        if (!con2d) {
            handle_error("Memory allocation for con2d failed.");
        }
    } else if (dimension == 3) {
        con3d = new voro::container_3d(supercell_edges[0][0], supercell_edges[1][0], supercell_edges[1][1],
                                    supercell_edges[2][0], supercell_edges[2][1], supercell_edges[2][2],
                                    n_x, n_y, n_z, true, true, true, 8, threads);
    }

    // Add all particles to the appropriate container
    if (dimension == 2) {
        for (int i = 0; i < number_of_particles; i++) {
            double* pp = particle_coordinates + 2 * i;
            con2d->put_parallel(i, *pp, pp[1]);
        }
        con2d->put_reconcile_overflow();
    } else if (dimension == 3) {
        for (int i = 0; i < number_of_particles; i++) {
            double* pp = particle_coordinates + 3 * i;
            con3d->put_parallel(i, *pp, pp[1], pp[2]);
        }
        con3d->put_reconcile_overflow();
    }


    ////////////////////////////////////////////////////
    ////
    ////   LOAD FILTER IF REQUESTED
    ////
    ////////////////////////////////////////////////////
    
    Filter filter;
    if(f_switch) filter.load_filter();

    
    ////////////////////////////////////////////////////
    ////
    ////   PROCESS DATA ACCORDING TO INPUT OPTIONS
    ////
    ////////////////////////////////////////////////////
    
    if(dimension==2) count_and_store_neighbors_2d(*con2d);

    // OUTPUTS VORONOI TOPOLOGY FOR EACH PARTICLE
    if(vt_switch)
    {
        if     (dimension==2) print_topology_vectors_2d(filename_data);
        else if(dimension==3) print_topology_vectors_3d(*con3d,filename_data);
    }
    
    // OUTPUTS DISTRIBUTION OF VORONOI TOPOLOGIES
    else if(d_switch)
    {
        if     (dimension==2) calc_distribution_2d(filter);        
        else if(dimension==3) calc_distribution_3d(*con3d, filter);
        filter.print_distribution(filename_data);
    }
        
    // OUTPUTS DISTRIBUTION OF VORONOI TOPOLOGIES
    else if(g_switch)
    {
        if     (dimension==2) calc_gaussian_distribution_2d(filter);
        else if(dimension==3) calc_gaussian_distribution_3d(filter);
        filter.print_distribution(filename_data);
    }
    
    // OUTPUTS LAMMPS DUMP FILE, INCLUDING CLASSIFICATION USING FILTER
    else if(f_switch)
    {
        if     (dimension==2) classify_particles_by_voronoi_topology_2d(filter);
        else if(dimension==3) classify_particles_by_voronoi_topology_3d(*con3d, filter);
        output_system(filename_data);
    }
    
    else if(u_switch || v_switch)  
    {
        if     (dimension==2) count_and_store_neighbors_2d(*con2d);
        else if(dimension==3) count_and_store_neighbors_3d(*con3d);
        pair_correlation_analysis();
    }
    
    else if(e_switch && dimension==2)
    {
        // COLORING SCHEMES 0,1,2,3 REQUIRE NO ANALYSIS;
        // COLORING SCHEME 4 REQUIRES LABELING PARTICLES ACCORDING TO RING NUMBER
        if(particle_coloring_scheme == 4) ring_coloring();
        output_eps(*con2d,filename_data);
    }
    
    else if(c_switch)
    {
        cluster_analysis();
    }
    
    // Clean up dynamically allocated containers
    if (con2d) {
        delete con2d;
        con2d = nullptr;
    }
    if (con3d) {
        delete con3d;
        con3d = nullptr;
    }
    
    cleanup();
    
    return 0;
}




