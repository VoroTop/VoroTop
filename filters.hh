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
#include <vector>


struct FilterEntry
{
    int type;
    int count;
    std::vector<int> wvector;
    
    FilterEntry() :                                          count(0)             {}
    FilterEntry(int c, std::vector<int> w) :                 count(c), wvector(w) {};
    FilterEntry(int t, int c, std::vector<int> w) : type(t), count(c), wvector(w) {};
};



inline bool compByCount(FilterEntry a, FilterEntry b)
{
    return a.count > b.count;
}

inline bool compByWvector(FilterEntry a, FilterEntry b)
{
    return a.wvector < b.wvector;
}

inline bool compByType(FilterEntry a, FilterEntry b)
{
    return a.type < b.type;
}

inline bool compByTypeThenCount(FilterEntry a, FilterEntry b)
{
    if(a.type == b.type) return a.count > b.count;
    else return a.type < b.type;
}



class Filter {
    std::vector<FilterEntry> entries;
    std::vector<std::string> structure_types;
    int max_file_filter_type;
    int max_filter_type;
    int file_filter_types;
    

public:
    Filter() : max_filter_type(0), max_file_filter_type(0) {}
    
    void loadFilter(std::string);  
    void print_distribution(std::string);
    void increment_or_add(std::vector<int> wvector) {increment_or_add(wvector,0);}
    void increment_or_add(std::vector<int>, int n);
    
    int  wvector_type (std::vector<int>);
    int  wvector_index(std::vector<int>);
    int  size(void) { return entries.size(); }
    int  get_max_type(void)    { return max_filter_type; }
    int  get_max_ff_type(void) { return max_file_filter_type; }
    
    int  get_entry_type (int n);
    int  get_entry_count(int n);
    void sort_by_wvector(void) { sort(entries.begin(), entries.end(), compByWvector); }
    void sort_by_count(void) { sort(entries.begin(), entries.end(), compByCount); }
    void sort_by_type(void) { sort(entries.begin(), entries.end(), compByType); }
    void sort_by_type_then_count(void) {sort(entries.begin(), entries.end(), compByTypeThenCount);}
    void sort_dist_by_count(void);
    void sort_for_clustering(void);
};


bool compByCount  (FilterEntry a, FilterEntry b);
bool compByWvector(FilterEntry a, FilterEntry b);


std::vector<int> calc_wvector (voro::voronoicell_base &c);
std::vector<int> calc_wvector (voro::voronoicell_base &c, bool extended);



#endif







