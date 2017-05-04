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

////   File: functions.hh


#ifndef __FUNCTIONS_H_INCLUDED__
#define __FUNCTIONS_H_INCLUDED__  


int  parse_header       (std::ifstream&);
void import_dump_file   (std::ifstream&, voro::particle_order &vo, voro::container_periodic &con);
void import_atomeye_file(std::ifstream&, voro::particle_order &vo, voro::container_periodic &con);

void help_message   (void);
void print_variables(void);
void parse_arguments(int argc, char *argv[]);

int  calc_all_wvectors (voro::container_periodic& con, voro::particle_order& vo, bool extended);
int  print_wvectors    (std::string filename);

void calc_distribution (Filter& filter);
void calc_gaussian_distribution(voro::container_periodic& con, voro::particle_order& vo, Filter& filter);

void create_cfg_file   (std::string filename, Filter& filter);
void cluster_analysis  (Filter& filter);


#endif







