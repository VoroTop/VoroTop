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

////   File: variables.cc


#include <vector>
#include <fstream>
#include <sstream>



int    number_of_particles;
int    file_type;      // TYPE OF INPUT FILE
                       //  1. LAMMPS DUMP
                       //  2. ATOMEYE

double origin[3];
double supercell_edges[3][3];

double hi_bound[3];
double xlo_bound, xhi_bound, xy;
double ylo_bound, yhi_bound, xz;
double zlo_bound, zhi_bound, yz;

bool scaled_coordinates;
bool no_velocity;

int header_lines;

std::ifstream data_file;

std::vector<std::string> attribute_labels;
std::vector<std::vector <int> > all_wvectors;
std::vector<std::vector <int> > neighbors_list;
std::vector<int> vt_structure_types;

std::vector<int> neighbors_list_c;
std::vector<int> cluster_index;
std::vector<int> cluster_sizes;
std::vector<int> resolved_types;
std::vector<double> volumes;

std::string name_of_data_file;
std::string filename_output;
std::string filename_filter;

int  xindex;
int  n_x,n_y,n_z;

bool f_switch;    // LOAD FILTER FILE
bool w_switch;    // PRINT W-VECTORS
bool d_switch;    // COMPUTE DISTRIBUTION OF W-VECTORS
bool g_switch;    // COMPUTE DISTRIBUTION BASED ON PERTURBATIONS OF SYSTEM
bool r_switch;    // RESOLVE INDETERMINATE TYPES
bool c_switch;    // PERFORM CLUSTER ANALYSIS
bool o_switch;    // SPECIFY OUTPUT FILE
bool od_switch;   // SPECIFY OUTPUT DIRECTORY

int     perturbation_samples;
double  perturbation_size;

int     clustering_default;
int     clustering_default_switch;

int     resolve_trials;




