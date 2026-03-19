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
bool l_switch;    // OUTPUT LAMMPS DUMP FILE
bool n_switch;    // DO NOT DRAW VORONOI CELLS; ONLY DRAW PARTICLES
bool u_switch;    // COMPUTE UNNORMALIZED VORONOI PAIR CORRELATION FUNCTION
bool v_switch;    // COMPUTE VORONOI PAIR CORRELATION FUNCTION
bool vt_switch;   // OUTPUT VORONOI TOPOLOGY FOR EACH PARTICLE, EITHER AS P- OR W-VECTOR
bool r_switch;    // RESOLVE INDETERMINATE TYPES

double xlo, xhi, xy;
double ylo, yhi, xz;
double zlo, zhi, yz;
bool   scaled_coordinates;
bool   triclinic_crystal_system;

int threads;
int number_of_particles;
int n_x,n_y,n_z;

int max_radius;
int particle_coloring_scheme;
int particles_in_eps;
int perturbation_samples;
double perturbation_size;

int file_format;
int header_lines;
int particle_attributes;

int index_id;
int index_type;
int index_species;
int index_x;
int index_y;
int index_z;

std::vector<char> column_types;

std::vector<int> header_assigned_types;
double coordinate_scale = 1.0;

std::string   filename_data;
std::string   filename_filter;
std::ifstream data_file;

int*    particle_ids;
int*    particle_types;
double* particle_coordinates;

std::vector<int> cluster_index;
std::vector<int> cluster_sizes;
std::vector<int> ring_index;

std::vector<int> cell_neighbor_count;
std::vector<std::vector <int> > list_of_neighbors;

int filter_structure_types;
std::vector<int> vt_structure_types;
std::vector<int> vt_structure_types_resolved;
