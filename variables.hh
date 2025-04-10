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

////   File: variables.hh


#pragma once
#include <string>

extern int  dimension;    // DIMENSION OF SYSTEM, 2 OR 3

extern bool c_switch;    // PERFORM CLUSTER ANALYSIS
extern bool d_switch;    // COMPUTE DISTRIBUTION OF P-VECTORS OR W-VECTORS
extern bool e_switch;    // OUTPUT EPS FILE
extern bool f_switch;    // LOAD FILTER FILE
extern bool g_switch;    // COMPUTE DISTRIBUTION BASED ON PERTURBATIONS OF SYSTEM
extern bool n_switch;    // DO NOT DRAW VORONOI CELLS; ONLY DRAW PARTICLES
extern bool u_switch;    // COMPUTE UNNORMALIZED VORONOI PAIR CORRELATION FUNCTION
extern bool v_switch;    // COMPUTE VORONOI PAIR CORRELATION FUNCTION
extern bool vt_switch;   // OUTPUT VORONOI TOPOLOGY FOR EACH PARTICLE, EITHER AS P- OR W-VECTOR
extern bool r_switch;    // RESOLVE INDETERMINATE TYPES

extern double hi_bound[3];
extern double origin[3];
extern double supercell_edges[3][3];
extern double xlo_bound, xhi_bound, xy;
extern double ylo_bound, yhi_bound, xz;
extern double zlo_bound, zhi_bound, yz;
extern double perturbation_size;
extern bool   scaled_coordinates;

extern int threads;
extern int number_of_particles;
extern int n_x,n_y,n_z;

extern int clustering_default;
extern int clustering_default_switch;
extern int max_radius;
extern int particle_coloring_scheme;
extern int particles_in_eps;
extern int perturbation_samples;

extern int header_lines;
extern int particle_attributes;

extern int index_id;
extern int index_x;
extern int index_y;
extern int index_z;

extern std::string   filename_data;
extern std::string   filename_filter;
extern std::ifstream data_file;

extern int*    particle_ids;
extern double* particle_coordinates;

extern std::vector<int> cluster_index;
extern std::vector<int> cluster_sizes;
extern std::vector<int> ring_index;

extern std::vector<int> cell_neighbor_count;
extern std::vector<std::vector <int> > neighbors_list_char;

extern std::vector<int> vt_structure_types;
extern std::vector<int> vt_structure_types_resolved;


