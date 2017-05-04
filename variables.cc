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

////   File: variables.cc


#include <vector>
#include <sstream>



int     timestep;
int     number_of_particles;
int     file_type;      // TYPE OF INPUT FILE
                        //  1. LAMMPS DUMP
                        //  2. ATOMEYE

double origin[3];
double xy, xz, yz;
double supercell_edges[3][3];

double cfg_lscale;
double cfg_atomic_mass;
std::string cfg_chem_symbol;
std::vector <double> xcoord;
std::vector <double> ycoord;
std::vector <double> zcoord;

bool scaled_coordinates;
bool no_velocity;

std::vector<std::string> attribute_labels;
std::vector<std::vector <int> > all_wvectors;
std::vector<std::vector <int> > neighbors_list;
std::vector<std::vector <double> > particle_data;
std::vector<int> neighbors_list_c;
std::vector<int> cluster_index;
std::vector<int> cluster_sizes;

std::string name_of_data_file;
std::string filename_output;
std::string filename_filter;

int  xindex;

int  n_x,n_y,n_z;

bool w_switch =0;  // PRINT WVECTORS
bool d_switch =0;  // CREATE DISTRIBUTION
bool g_switch =0;  // CREATE EINSTEIN DISTRIBUTION
bool df_switch=0;  // CREATE CFG USING DISTRIBUTION FILTER
bool c_switch =0;  // CREATE CFG USING JUST CLUSTERING TO REMOVE ``SMALL'' DEFECTS
bool f_switch =0;  // CREATE CFG USING SPECIFIED FILTER FILE
bool o_switch =0;  // SPECIFY OUTPUT FILE
bool od_switch=0;  // SPECIFY OUTPUT DIRECTORY
bool p_switch =0;  // PRINT VARIABLES
bool h_switch =0;  // PRINT HELP MENU

int     perturbation_samples;
double  perturbation_size;





