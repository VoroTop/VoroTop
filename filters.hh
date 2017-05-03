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

////   File: filters.hh


#ifndef __FILTERS_H_INCLUDED__
#define __FILTERS_H_INCLUDED__   


#include "voro++.hh"
#include <map>
#include <vector>
#include <algorithm>


struct FilterEntry
{
    int type;
    int count;
    bool source; // 0 filter file, 1 data
    
    FilterEntry(int t) : type(t) {};
    FilterEntry(int t, int c, bool s) : type(t), count(c), source(s) {};
};



class Filter {
public:
    std::map<std::vector<int>,FilterEntry> entries;
    std::vector<std::string> structure_types;
    int max_file_filter_type;
    int max_filter_type;
    int file_filter_types;
    
    
public:
    Filter() : max_file_filter_type(0), max_filter_type(0) {}
    
    void loadFilter(std::string);
    void print_distribution(std::string);

    void increment_or_add(std::vector<int> wvector) {increment_or_add(wvector,0);}
    void increment_or_add(std::vector<int>, int n);
    
    int  wvector_type(std::vector<int>);

    int  get_max_type(void)    { return max_filter_type; }
    int  get_max_ff_type(void) { return max_file_filter_type; }
};



unsigned int countWordsInString(std::string const& str);

std::vector<int> calc_wvector (voro::voronoicell_base &c);
std::vector<int> calc_wvector (voro::voronoicell_base &c, bool extended);



#endif







