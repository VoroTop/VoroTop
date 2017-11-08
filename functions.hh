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

////   File: functions.hh


#ifndef __FUNCTIONS_H_INCLUDED__
#define __FUNCTIONS_H_INCLUDED__  


void help_message       (void);
void parse_arguments    (int argc, char *argv[]);

int  parse_header       (std::ifstream&);
void import_dump_file   (std::ifstream&, voro::container_periodic &con, voro::particle_order &vo);
void import_atomeye_file(std::ifstream&, voro::container_periodic &con, voro::particle_order &vo);

int  calc_all_wvectors (voro::container_periodic& con, voro::particle_order& vo, bool extended);
void calc_structure_types(Filter &filter);
void calc_distribution (Filter& filter);
void cluster_analysis  (Filter& filter);
void calc_gaussian_distribution (voro::container_periodic& con, voro::particle_order& vo, Filter& filter);
void resolve_indeterminate_types(voro::container_periodic& con, voro::particle_order& vo, Filter &filter);

int  print_wvectors    (std::string filename);
void output_system     (std::string filename, Filter &filter);

#endif







