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

////   File: vorotop.hh


#ifndef __VOROTOP_H_INCLUDED__
#define __VOROTOP_H_INCLUDED__   


#include "voro++.hh"
#include <string>


extern int  timestep;
extern int  number_of_particles;
extern int  file_type;
extern int  xindex;

extern double origin[3];
extern double xy, xz, yz;
extern double supercell_edges[3][3];
extern double cfg_lscale;
extern double cfg_atomic_mass;

extern std::string cfg_chem_symbol;

extern std::vector <double> xcoord;
extern std::vector <double> ycoord;
extern std::vector <double> zcoord;

extern bool scaled_coordinates;
extern bool no_velocity;

extern std::string name_of_data_file;
extern std::string filename_output;
extern std::string filename_filter;

extern std::vector<std::string> attribute_labels;

extern std::vector<std::vector <int> > all_wvectors;
extern std::vector<std::vector <int> > neighbors_list;
extern std::vector<std::vector <double> > particle_data;     // STORE PARTICLE DATA

extern std::vector<int> neighbors_list_c;
extern std::vector <int> cluster_index;
extern std::vector <int> cluster_sizes;


extern int n_x,n_y,n_z;

extern bool w_switch;    // PRINT W-VECTORS
extern bool d_switch;    // CREATE DISTRIBUTION OF W-VECTORS
extern bool g_switch;   // CREATE EINSTEIN DISTRIBUTION
extern bool f_switch;    // USE FILTER FILE
extern bool df_switch;   // CREATE CFG USING DISTRIBUTION FILTER
extern bool c_switch;    // CREATE CFG FOR BCC USING CLUSTER-BASED FILTERING
extern bool o_switch;    // SPECIFY OUTPUT FILE
extern bool od_switch;   // SPECIFY OUTPUT DIRECTORY
extern bool p_switch;    // PRINT VARIABLES
extern bool h_switch;    // PRINT HELP MENU

extern int     perturbation_samples;
extern double  perturbation_size;



#endif







