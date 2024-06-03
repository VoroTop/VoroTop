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

////   File: functions.hh


#ifndef __FUNCTIONS_H_INCLUDED__
#define __FUNCTIONS_H_INCLUDED__  

void help_message    (void);
void parse_arguments (int argc, char *argv[]);
void parse_header    (std::ifstream&);
void import_data     (std::ifstream&);

void count_and_store_neighbors(voro::container_2d& con);
void count_and_store_neighbors(voro::container_3d& con);
void calc_distribution (voro::container_2d& con, Filter &filter);
void calc_distribution (voro::container_3d& con, Filter &filter);

void calc_gaussian_distribution  (Filter& filter);
void calc_gaussian_distribution3d(Filter& filter);
void resolve_indeterminate_types(voro::container_3d& con, Filter &filter);
void cluster_analysis  (Filter& filter);

void print_topology_vectors2d(std::string filename);
int  print_topology_vectors(voro::container_3d& con, std::string filename);

void cluster_analysis         (void);
void defect_cluster_analysis  (void);
void pair_correlation_analysis(void);

void output_system    (std::string filename, Filter &filter);
void output_eps       (voro::container_2d& con, std::string filename);
void ring_coloring    ();

void classify_particles_by_voronoi_topology(Filter &filter);
int  classify_particles_by_voronoi_topology(voro::container_3d& con, Filter &filter);

int  compute_canonical_code(vector<int>& canonical_code, voro::voronoicell_3d& vcell);

#endif




