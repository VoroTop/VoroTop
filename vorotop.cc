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
    if (particle_types) {
        delete[] particle_types;
        particle_types = nullptr;
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
        std::cerr << "Error opening file: " << filename_data << std::endl;
        exit(0);
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
        exit(0);
    }

    particle_types      =new    int[number_of_particles];
    if (!particle_types) {
        handle_error("Memory allocation for particle_types failed.");
        exit(0);
    }

    particle_coordinates=new double[number_of_particles*dimension];
    if (!particle_coordinates) {
        handle_error("Memory allocation for particle_coordinates failed.");
        exit(0);
    }

    vt_structure_types.resize      (number_of_particles);  // MEMORY FOR STRUCTURE TYPES

    if(dimension==2 || u_switch || v_switch || c_switch)
    {
        list_of_neighbors.resize   (number_of_particles);  // MEMORY FOR LIST OF NEIGHBORS
        cell_neighbor_count.resize (number_of_particles);
    }
    
    if(c_switch)                                          // MEMORY FOR CLUSTER ANALYSIS
    {
        cluster_index.resize(number_of_particles);
        cluster_sizes.resize(number_of_particles);
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
    
    if(d_switch==0 && g_switch==0 && mf_switch==0 && vt_switch==0 && c_switch==0 &&
       l_switch==0 &&  u_switch==0 && v_switch==0 && e_switch==0)
    {
        if(f_switch)
        {
            std::cout << "Filter file specified, but no analysis requested.\n";
            std::cout << "By default will output LAMMPS dump file as if chosen using -l option.\n";
            l_switch = 1;
        }
        else
        {
            help_message();
            std::cout << "No valid options specified.\n";
        }
    }
        
    
    // IMPORTS ALL COORDINATE DATA
    import_data();
    
    ////////////////////////////////////////////////////////
    ////  CREATE VORO++ CONTAINER WITH NECESSARY DIMENSIONS;
    ////  THEN ADD PARTICLES TO THE CONTAINER
    ////////////////////////////////////////////////////////
    
    // Declare pointers for the containers
    voro::container_2d*        con2d  = nullptr;
    voro::container_3d*        con3d  = nullptr;
    voro::container_triclinic* contri = nullptr;

    // Dynamically allocate the appropriate container based on dimension and geometry
    if (dimension == 2)
    {
        con2d = new voro::container_2d(xlo,xhi,ylo,yhi,n_x,n_y,true,true,4,threads);
        if (!con2d) {
            handle_error("Memory allocation for con2d failed.");
        }
    }
    else if (triclinic_crystal_system)
    {
        double bx = xhi - xlo, by_val = yhi - ylo, bz_val = zhi - zlo;
        contri = new voro::container_triclinic(bx,xy,by_val,xz,yz,bz_val,n_x,n_y,n_z,8,threads);
        if (!contri) {
            handle_error("Memory allocation for contri failed.");
        }
    }
    else
    {
        con3d = new voro::container_3d(xlo,xhi,ylo,yhi,zlo,zhi,n_x,n_y,n_z,true,true,true,8,threads);
        if (!con3d) {
            handle_error("Memory allocation for con3d failed.");
        }
    }

    // Add all particles to the appropriate container
    if (dimension == 2) {
        for (int i = 0; i < number_of_particles; i++) {
            double* pp = particle_coordinates + 2 * i;
            con2d->put_parallel(i, *pp, pp[1]);
        }
        con2d->put_reconcile_overflow();
    } else if (contri) {
        for (int i = 0; i < number_of_particles; i++) {
            double* pp = particle_coordinates + 3 * i;
            contri->put_parallel(i, *pp, pp[1], pp[2]);
        }
        contri->put_reconcile_overflow();
    } else if (con3d) {
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
    
    // HELPER LAMBDA: DISPATCH A 3D OPERATION TO THE CORRECT CONTAINER TYPE
    // ACCEPTS A CALLABLE THAT TAKES A CONTAINER BY REFERENCE
    auto dispatch_3d = [&](auto&& func) {
        if(contri) func(*contri);
        else       func(*con3d);
    };

    // THESE CALCULATIONS ARE REQUIRED FOR ALL OPTIONS IN 2-DIMENSIONAL SYSTEMS
    if(dimension==2) count_and_store_neighbors_2d(*con2d);

    // OUTPUTS VORONOI TOPOLOGY FOR EACH PARTICLE
    if(vt_switch)
    {
        if     (dimension==2) print_topology_vectors_2d(filename_data);
        else if(dimension==3) dispatch_3d([&](auto& con){ print_topology_vectors_3d(con, filename_data); });
    }

    // GENERATES FILTER FILE FROM CRYSTAL STRUCTURE
    else if(mf_switch)
    {
        if(perturbation_samples > 0)  // PERTURBATION MODE VIA CELL REPLICATION
        {
            if     (dimension==2) calc_gaussian_distribution_2d(filter);
            else if(dimension==3) calc_gaussian_distribution_3d(filter);
        }
        else  // EXACT MODE
        {
            if     (dimension==2) calc_distribution_2d(filter);
            else if(dimension==3) dispatch_3d([&](auto& con){ calc_distribution_3d(con, filter); });
        }
        filter.print_filter(filename_data);
    }

    // OUTPUTS DISTRIBUTION OF VORONOI TOPOLOGIES
    else if(d_switch)
    {
        if     (dimension==2) calc_distribution_2d(filter);
        else if(dimension==3) dispatch_3d([&](auto& con){ calc_distribution_3d(con, filter); });
        filter.print_distribution(filename_data);
    }

    // OUTPUTS DISTRIBUTION OF VORONOI TOPOLOGIES
    else if(g_switch)
    {
        if     (dimension==2) calc_gaussian_distribution_2d(filter);
        else if(dimension==3) calc_gaussian_distribution_3d(filter);
        filter.print_distribution(filename_data);
    }

    // OUTPUTS PAIR CORRELATION ANALYSIS
    else if(u_switch || v_switch)
    {
        if(dimension==3) dispatch_3d([](auto& con){ count_and_store_neighbors_3d(con); });
        pair_correlation_analysis();
    }

    // OUTPUTS CLUSTER ANALYSIS
    else if(c_switch && !l_switch && !e_switch)
    {
        if     (dimension==2) classify_particles_by_voronoi_topology_2d(filter);
        else if(dimension==3) dispatch_3d([&](auto& con){ classify_particles_by_voronoi_topology_3d(con, filter); });
        if(dimension==3) dispatch_3d([](auto& con){ count_and_store_neighbors_3d(con); });
        cluster_analysis();
    }

    // OUTPUTS LAMMPS DUMP FILE, INCLUDING CLASSIFICATION USING FILTER
    if(l_switch)
    {
        if     (dimension==2) classify_particles_by_voronoi_topology_2d(filter);
        else if(dimension==3) dispatch_3d([&](auto& con){ classify_particles_by_voronoi_topology_3d(con, filter); });
        if(c_switch)
        {
            if(dimension==3) dispatch_3d([](auto& con){ count_and_store_neighbors_3d(con); });
            cluster_analysis();
        }
        output_lammps_dump(filename_data);
    }

    else if(e_switch)
    {
        // COLORING SCHEMES 3 REQUIRES CLASSIFYING PARTICLES USING FILTER
        if(particle_coloring_scheme == 3) classify_particles_by_voronoi_topology_2d(filter);

        // COLORING SCHEME 4 REQUIRES LABELING PARTICLES ACCORDING TO RING NUMBER
        if(particle_coloring_scheme == 4) ring_coloring();

        if(particle_coloring_scheme == 5 || particle_coloring_scheme == 6 || particle_coloring_scheme == 7 || particle_coloring_scheme == 8)
        {
            classify_particles_by_voronoi_topology_2d(filter);
            cluster_analysis();
        }

        output_eps(*con2d,filename_data);
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
    if (contri) {
        delete contri;
        contri = nullptr;
    }
    
    cleanup();
    
    return 0;
}




