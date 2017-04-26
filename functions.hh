////////////////////////////////////////////////////////
////                                                ////
////   ******************************************   ////
////   *                 VTop++                 *   ////
////   *     Complete Voronoi Cell Topology     *   ////
////   *   Analysis and Visualization Toolkit   *   ////
////   *             (Version 0.3)              *   ////
////   *                                        *   ////
////   *           Emanuel A. Lazar             *   ////
////   *      University of Pennsylvania        *   ////
////   *           December 5, 2014             *   ////
////   ******************************************   ////
////                                                ////
////////////////////////////////////////////////////////

////   File: functions.hh


#ifndef __FUNCTIONS_H_INCLUDED__
#define __FUNCTIONS_H_INCLUDED__  


unsigned int countWordsInString(std::string const& str);

void help_message   (void);
void print_variables(void);
void parse_arguments(int argc, char *argv[]);

int  calc_all_wvectors (container_periodic& con, particle_order& vo, bool extended);
int  print_wvectors    (std::string filename);

void calc_distribution (Filter& filter);
void calc_gaussian_distribution(voro::container_periodic& con, voro::particle_order& vo, Filter& filter);

void create_cfg_file   (std::string filename, Filter& filter);
void cluster_analysis  (Filter& filter);

void create_cfg_file_with_perturbations(std::string filename, container_periodic& con, particle_order& vo, Filter& filter);


#endif







