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


#pragma once

#include <vector>

struct VoronoiTopology {
    std::vector<int> canonical_code;
    int chirality;
    int symmetry_counter;
    int face_count;
    int max_face_edges;
    std::vector<int> face_edge_counts;  // number of faces with each number of edges
};

void help_message    (void);
void parse_command_line_options (int argc, char *argv[]);
void parse_header    (std::ifstream&);
void import_data     ();

void count_and_store_neighbors_2d(voro::container_2d& con);
void count_and_store_neighbors_3d(voro::container_3d& con);

void calc_distribution_2d (Filter &filter);
void calc_distribution_3d (voro::container_3d& con, Filter &filter);

void calc_gaussian_distribution_2d(Filter& filter);
void calc_gaussian_distribution_3d(Filter& filter);

void print_topology_vectors_2d(std::string filename);
void print_topology_vectors_3d(voro::container_3d& con, std::string filename);

void cluster_analysis         (void);
void cluster_analysis2        (void);

void defect_cluster_analysis  (void);   // EXPERIMENTAL
void pair_correlation_analysis(void);

void output_lammps_dump (std::string filename);
void output_eps         (voro::container_2d& con, std::string filename);
void ring_coloring      ();

void classify_particles_by_voronoi_topology_2d(Filter &filter);
int  classify_particles_by_voronoi_topology_3d(voro::container_3d& con, Filter &filter);

int  compute_canonical_code_2d(vector<int>& canonical_code, int pid);
VoronoiTopology compute_canonical_code_3d(voro::voronoicell_3d& vcell);

void validate_max_radius(int max_radius);




