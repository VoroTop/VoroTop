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

////   File: import.hh


#ifndef __IMPORT_H_INCLUDED__
#define __IMPORT_H_INCLUDED__  


#include "voro++.hh"

using namespace voro;


int  parse_header       (std::ifstream&);
void import_dump_file   (std::ifstream&, particle_order &vo, container_periodic &con);
void import_atomeye_file(std::ifstream&, particle_order &vo, container_periodic &con);

#endif







