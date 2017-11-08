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

////   File: filters.hh


#ifndef __FILTERS_H_INCLUDED__
#define __FILTERS_H_INCLUDED__   


#include "voro++.hh"
#include <map>
#include <vector>
#include <algorithm>


using WeinbergVector = std::vector<int>;

struct FilterEntry
{
    int type;
    int count;
    int source; // 0 filter file, 1 data
    
    FilterEntry(int t)               : type(t), count(0), source(0) {};
    FilterEntry(int t, int c, int s) : type(t), count(c), source(s) {};
};


class Filter {
    std::map<WeinbergVector,FilterEntry> entries;
    std::vector<std::string> structure_types;
    
    std::vector<bool>               indeterminate;    // TRACK WHETHER TYPE IS INDETERMINATE
public:
    std::vector<std::pair<int,int> > resolved_types;   // PRIMARY AND SECONDARY RESOLVED TYPES

    int max_file_filter_type;
    int max_filter_type;
    int file_filter_types;
    
    
public:
    Filter() : max_file_filter_type(0), max_filter_type(0) {}
    
    void loadFilter(std::string);
    void print_distribution(std::string);
    void relabel_data_types();
    
    int  wvector_type    (WeinbergVector);
    void increment_or_add(WeinbergVector, int n);

    int  get_max_type   (void) { return max_filter_type; }
    int  get_max_ff_type(void) { return max_file_filter_type; }
    
    std::map<WeinbergVector,FilterEntry>::iterator find(WeinbergVector w) { return entries.find(w);   }
    std::map<WeinbergVector,FilterEntry>::iterator end (void)             { return entries.end ( );   }
    int get_file_filter_types(void)                                       { return file_filter_types; }
    bool is_indeterminate(int c)                                          { return indeterminate[c];  }
};


WeinbergVector calc_wvector (voro::voronoicell_base &c);
WeinbergVector calc_wvector (voro::voronoicell_base &c, bool extended);


#endif







