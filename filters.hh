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

////   File: filters.hh


// A FILTER SPECIFICIES TOPOLOGICAL TYPES ASSOCIATED TO PARTICULAR
// LOCAL STRUCTURE TYPES.  EACH TOPOLOGICAL TYPE IS RECORDED WITH
// A P-VECTOR IN TWO DIMENSIONS OR A WEINBERG CODE IN THREE.  IF
// WE KNOW WHICH TOPOLOGICAL TYPES ARE ASSOCIATED WITH A GIVEN BULK
// CRYSTAL, THEN PARTICLES ASSOCIATED TO DEFECTS CAN BE IDENTIFIED
// AS THOSE WITH OTHER TOPOLOGICAL TYPES.

// VoroTop ALLOWS A USER TO SPECIFY A PARTICULAR FILTER TO STUDY
// DATA.  FOR EXAMPLE, WHEN STUDYING AN FCC POLYCRYSTAL, ONE MIGHT
// USE A FILTER SPECIFYING ALL TOPOLOGICAL TYPES ASSOCIATED WITH
// FINITE-TEMPERATURE FCC AND HCP CRYSTALS.  THIS WILL ENABLE THE
// IDENTIFICATION OF FCC STRUCTURE ASSOCIATED WITH THE BULK
// AND HCP STRUCTURE ASSOCIATED WITH DISLOCATIONS, TWIN PLANES,
// AND STACKING FAULTS.

// FILTER FILES FOR COMMON STRUCTURES CAN BE FOUND ONLINE, AT
//    https://www.vorotop.org

// FILTER FILES ARE ORGANIZED AS FOLLOWS:
//
//   LINES BEGINNING WITH # ARE IGNORED AND ARE USED FOR COMMENTS.
//
//   LINES BEGINNING WITH * SPECIFICY STRUCTURE TYPES IN THE FITLER.
//   EACH SUCH LINE, AFTER THE *, INCLUDES AN INTEGER AND THEN A PLAIN-
//   TEXT DESCRIPTION OF THE TYPE.  STRUCTURE TYPES ARE ALWAYS LISTED
//   IN INCREASING ORDER, AND BEGIN COUNTING FROM 1.
//
//   REMAINING LINES RECORD WEINBERG CODES AND ASSOCIATED STRUCTURE
//   TYPES.  EACH LINE BEGINS WITH AN INTEGER AND IS FOLLOWED BY A
//   SEQUENCE OF INTEGERS (A WEINBERG CODE) REPRESENTING A TOPOLOGICAL
//   TYPE.  THE FORMAT OF THE WEINBERG CODES IS THAT FOUND IN THE WORK
//   OF WEINBERG; SEE DOCUMENTATION FOR REFERENCES TO APPROPRIATE WORK.


#ifndef __FILTERS_H_INCLUDED__
#define __FILTERS_H_INCLUDED__   


#include "voro++.hh"
#include <map>
#include <vector>
#include <algorithm>


using TopologyVector = std::vector<int>;

// EACH FILTER ENTRY CORRESPONDS TO A TOPOLOGICAL TYPE, GIVEN EITHER BY
// A P-VECTOR OR A WEINBERG-VECTOR.  WE SAVE ITS ASSOCIATED STRUCTURE
// TYPE OR TYPES, ITS COUNTS (ACCOUNTING FOR DIFFERENCES IN CHIRALITY),
// AND ALSO THE SOURCE OF THE TYPE, WHETHER FROM A FILTER DATA FILE OR
// OR DATA FROM A PARTICLE SYSTEM.
struct FilterEntry
{
    std::vector<int> structure_types;  // TRACKS ASSOCIATED STRUCTURE TYPES. IF > 1, THEN TOPOLOGY IS INDETERMINATE

    int total;      // total count, sum of counts[3]
    int counts[3];  // 0 for left-handed, 1 for chiral, and 2 for right-handed
    int source;     // 0 filter file, 1 data
    
    FilterEntry(int type)                           : source(0) { counts[0]=0; counts[1]=0; counts[2]=0; structure_types.push_back(type); total=0;};
    FilterEntry(int type, int chirality, int count) : source(1) { counts[0]=0; counts[1]=0; counts[2]=0; total=count;
        counts[chirality+1]=count;
        structure_types.push_back(0);
    };
    FilterEntry(int type, int left, int non, int right) {
        structure_types.push_back(0);
        counts[0]=left;
        counts[1]=non;
        counts[2]=right;
        total=left+non+right;
    }
};


// EACH FILTER IS A LIST OF VORONOI TOPOLOGIES, AND SOME
class Filter {
    std::map<TopologyVector,FilterEntry> entries;

public:
    int max_file_filter_type;
    int max_filter_type;
    void copy_filter(Filter to_copy);
    
public:
    Filter() : max_file_filter_type(0), max_filter_type(0) {}
    
    void load_filter(std::string);               // FILE NAME OF FILTER
    void print_distribution(std::string);       // PRINT DISTRIBUTION TO FILE
    void relabel_data_types();                  // FIX: DON'T REMEMBER WHAT THIS IS FOR
    
    int  vt_structure_type(TopologyVector);     // RETURNS ASSOCIATED STRUCTURE TYPE
    void increment_or_add (TopologyVector);     // ADDS VORONOI TOPOLOGY TO FILTER
    void increment_or_add (TopologyVector, int chirality, int count);

    int  get_max_ff_type(void) { return max_file_filter_type; }
    
    std::map<TopologyVector,FilterEntry>::iterator find(TopologyVector w) { return entries.find(w);   }
    std::map<TopologyVector,FilterEntry>::iterator end (void)             { return entries.end ( );   }
    int count_indeterminate_types();
};

#endif




