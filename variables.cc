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

////   File: variables.cc

#include <vector>
#include <fstream>
#include <sstream>

int  dimension;    // DIMENSION OF SYSTEM, 2 OR 3

bool c_switch;    // PERFORM CLUSTER ANALYSIS
bool d_switch;    // COMPUTE DISTRIBUTION OF P-VECTORS OR W-VECTORS
bool e_switch;    // OUTPUT EPS FILE
bool f_switch;    // LOAD FILTER FILE
bool g_switch;    // COMPUTE DISTRIBUTION BASED ON PERTURBATIONS OF SYSTEM
bool n_switch;    // DO NOT DRAW VORONOI CELLS; ONLY DRAW PARTICLES
bool u_switch;    // COMPUTE UNNORMALIZED VORONOI PAIR CORRELATION FUNCTION
bool v_switch;    // COMPUTE VORONOI PAIR CORRELATION FUNCTION
bool vt_switch;   // OUTPUT VORONOI TOPOLOGY FOR EACH PARTICLE, EITHER AS P- OR W-VECTOR
bool r_switch;    // RESOLVE INDETERMINATE TYPES

double hi_bound[3];
double origin[3];
double supercell_edges[3][3];
double xlo_bound, xhi_bound, xy;
double ylo_bound, yhi_bound, xz;
double zlo_bound, zhi_bound, yz;
double perturbation_size;
bool   scaled_coordinates;

int threads;
int number_of_particles;
int n_x,n_y,n_z;

int clustering_default;
int clustering_default_switch;
int max_radius;
int particle_coloring_scheme;
int particles_in_eps;
int perturbation_samples;
int resolve_trials;

int header_lines;
int particle_attributes;
int index_id;
int index_x;
int index_y;
int index_z;

std::string   filename_data;
std::string   filename_filter;
std::ifstream data_file;

int*    particle_ids;
double* particle_coordinates;

std::vector<double> areas;
std::vector<double> volumes;

std::vector<int> cluster_index;
std::vector<int> cluster_sizes;
std::vector<int> ring_index;

std::vector<int> cell_neighbor_count;
std::vector<std::vector <int> > neighbors_list_char;

std::vector<int> vt_structure_types;
std::vector<int> vt_structure_types_resolved;
