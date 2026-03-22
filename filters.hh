////////////////////////////////////////////////////////
////                                                ////
////   ******************************************   ////
////   *                                        *   ////
////   *     VoroTop: Voronoi Cell Topology     *   ////
////   *   Visualization and Analysis Toolkit   *   ////
////   *             (Version 1.1)              *   ////
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


#pragma once

#include "voro++.hh"
#include <map>
#include <vector>
#include <algorithm>

using VoronoiTopologyVector = std::vector<int>;

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
    
    FilterEntry(int type)                           { counts[0]=0; counts[1]=0; counts[2]=0; structure_types.push_back(type); total=0;};
    FilterEntry(int chirality, int count)           { counts[0]=0; counts[1]=0; counts[2]=0; total=count;
        counts[chirality+1]=count;
        structure_types.push_back(0);
    };
    FilterEntry(int left_handed, int non_handed, int right_handed) {
        // Initialize structure_types with 0 to indicate an indeterminate or default structure type.
        structure_types.push_back(0);
        counts[0]= left_handed;
        counts[1]=  non_handed;
        counts[2]=right_handed;
        total=left_handed+non_handed+right_handed;
    }
};


// EACH FILTER IS A LIST OF VORONOI TOPOLOGIES, AND SOME
class Filter {
    std::map<VoronoiTopologyVector,FilterEntry> entries;

public:
    int max_filter_type_from_file;
    void copy_filter(const Filter& to_copy);
    int count_indeterminate_types();
    
    Filter() : max_filter_type_from_file(0) {}

    // LOADS A FILTER FROM A FILE
    // The method takes a string representing the file name of the filter
    // and loads the filter data into the entries map. 
    void load_filter();               // FILE NAME OF FILTER

    // Prints the distribution to a file. The parameter specifies the file path where the output will be written.
    void print_distribution(std::string file_path);       // PRINT DISTRIBUTION TO FILE
    void print_filter(std::string file_path);              // PRINT FILTER TO FILE

    // RETURNS ASSOCIATED STRUCTURE TYPE:
    // The method returns an integer representing the structure type
    // associated with the given Voronoi topology vector. If no match
    // is found, it returns 0 to indicate an indeterminate type.
    int  vt_structure_type(VoronoiTopologyVector);

    void increment_or_add (VoronoiTopologyVector);     // ADDS VORONOI TOPOLOGY TO FILTER
    void increment_or_add (VoronoiTopologyVector& topology_vector, int chirality, int count);

    std::map<VoronoiTopologyVector,FilterEntry>::iterator find(VoronoiTopologyVector w) { return entries.find(w);   }
    std::map<VoronoiTopologyVector,FilterEntry>::iterator end (void)                    { return entries.end ( );   }
};



