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

////   File: variables.hh


#ifndef __VOROTOP_H_INCLUDED__
#define __VOROTOP_H_INCLUDED__   


#include <string>


extern int  number_of_particles;
extern int  file_type;
extern int  xindex;

extern std::ifstream data_file;

extern double origin[3];
extern double supercell_edges[3][3];

extern bool scaled_coordinates;
extern bool no_velocity;

extern double hi_bound[3];
extern double xlo_bound, xhi_bound, xy;
extern double ylo_bound, yhi_bound, xz;
extern double zlo_bound, zhi_bound, yz;

extern int header_lines;

extern std::string name_of_data_file;
extern std::string filename_output;
extern std::string filename_filter;

extern std::vector<std::string> attribute_labels;

extern std::vector<std::vector <int> > all_wvectors;
extern std::vector<std::vector <int> > neighbors_list;
extern std::vector<int> vt_structure_types;

extern std::vector <int> neighbors_list_c;
extern std::vector <int> cluster_index;
extern std::vector <int> cluster_sizes;
extern std::vector <int> resolved_types;
extern std::vector <double> volumes;

extern int n_x,n_y,n_z;

extern bool f_switch;    // LOAD FILTER FILE
extern bool w_switch;    // PRINT W-VECTORS
extern bool d_switch;    // COMPUTE DISTRIBUTION OF W-VECTORS
extern bool g_switch;    // COMPUTE DISTRIBUTION BASED ON PERTURBATIONS OF SYSTEM
extern bool r_switch;    // RESOLVE INDETERMINATE TYPES
extern bool c_switch;    // PERFORM CLUSTER ANALYSIS
extern bool o_switch;    // SPECIFY OUTPUT FILE
extern bool od_switch;   // SPECIFY OUTPUT DIRECTORY

extern int    perturbation_samples;
extern double perturbation_size;

extern int    clustering_default;
extern int    clustering_default_switch;

extern int    resolve_trials;

#endif







