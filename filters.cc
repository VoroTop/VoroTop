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

////   File: filters.cc


// A FILTER IS A SPECIFICATION OF TOPOLOGICAL TYPES, RECORDED
// WITH WEINBERG CODES, ASSOCIATED TO PARTICULAR LOCAL STRUCTURE
// TYPES.  IF WE KNOW WHICH TOPOLOGICAL TYPES ARE ASSOCIATED WITH
// A GIVEN BULK CRYSTAL, FOR EXAMPLE, THEN PARTICLES ASSOCIATED
// TO DEFECTS CAN BE IDENTIFIED AS THOSE WITH OTHER TOPOLOGICAL
// TYPES.

// VoroTop ALLOWS A USER TO SPECIFICY A PARTICULAR FILTER TO STUDY
// DATA.  FOR EXAMPLE, WHEN STUDYING AN FCC POLYCRYSTAL, ONE MIGHT
// USE A FILTER SPECIFYING ALL TOPOLOGICAL TYPES ASSOCIATED WITH
// FINITE-TEMPERATURE FCC AND HCP CRYSTALS.  THIS WILL THEN ALLOW
// FOR THE IDENTIFICATION OF FCC STRUCTURE (ASSOCIATED WITH THE
// BULK CRYSTAL), HCP STRUCTURE (E.G., ASSOCIATED WITH DISLOCATIONS,
// TWIN PLANES, STACKING FAULTS), AND OTHER DEFECTS.

// FILTER FILES FOR SEVERAL COMMON STRUCTURES CAN BE FOUND ONLINE,
// AT https://www.vorotop.org

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
//   SEQUENCE OF INTEGERS (WEINBERG CODE) REPRESENTING A TOPOLOGICAL
//   TYPE.  THE FORMAT OF THE WEINBERG CODES IS THAT FOUND IN THE WORK
//   OF WEINBERG; SEE DOCUMENTATION FOR REFERENCES TO APPROPRIATE WORK.


#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "filters.hh"
#include "variables.hh"



////////////////////////////////////////////////////
////
////   IMPORTS FILTER FROM FILE
////
////////////////////////////////////////////////////

////   THREE TYPES OF LINES:
////   1) COMMENTS, BEGIN WITH #
////   2) SPECIFICATION OF STRUCTURE TYPES, BEGIN WITH *
////   3) FILTER, STRUCTURE TYPES AND WEINBERG CODES

void Filter::loadFilter(std::string filename)
{
    std::ifstream filter_file(filename_filter);
    if(!filter_file.is_open())
    {
        std::cerr << "Failed to open input file " << filename_filter << "\n\n";
        exit(-1);
    }
    
    structure_types.push_back(std::string());   // STRUCTURE TYPES WILL BE INDEXED BY TYPE, BEGINNING FROM 1.
    
    std::string line;
    max_file_filter_type = 0;
    int line_counter = 0;
    
    while ( getline(filter_file, line) )
    {
        std::istringstream is(line);
        
        line_counter++;
        if(line[0] == '#') continue;    // IGNORE COMMENTED LINES
        else if(line[0] == '*')         // ADD STRUCTURE TYPES
        {
            int index;
            std::string name;
            
            is.ignore();  // IGNORES *
            is >> index;
            is >> name;
            
            if(index != ++max_file_filter_type)
            {
                std::cout << "Invalid filter provided.  Structure types in filter must \n";
                std::cout << "be ordered consecutively beginning from 1.\n";
                exit(-1);
            }
            
            structure_types.push_back(name);
            continue;
        }
        
        else                            // ADD WEINBERG CODES
        {
            int type;
            is >> type;
            
            if(type < 1 || type > max_file_filter_type)
            {
                std::cout << "Invalid filter provided.  Structure type " << type << " on line " << line_counter << " does\n";
                std::cout << "not match any of those specified in filter header.\n";
                exit(-1);
            }
            
            std::vector<int> new_wvector;
            is.ignore(256,'(');
            int next_label;
            while (is >> next_label)
            {
                new_wvector.push_back(next_label);
                if (is.peek() == ',') is.ignore();
                if (is.peek() == ')') break;  // IGNORES DATA AFTER FINAL ')'
            }
            
            entries.insert({std::move(new_wvector), FilterEntry(type)});
        }
    }
    
    file_filter_types = entries.size();
    max_filter_type   = max_file_filter_type;
}



bool compare_by_count_wvector(std::map<std::vector<int>,FilterEntry>::iterator a, std::map<std::vector<int>,FilterEntry>::iterator b)
{
    if(a->second.count == b->second.count)
        return a->first < b->first;
    else return a->second.count > b->second.count;
}



void Filter::increment_or_add(std::vector<int> wvector, int count)
{
    auto it = entries.find(wvector);
    if (it != entries.end()) it->second.count += count;
    else entries.insert({std::move(wvector), FilterEntry(0, count, 1)});
}



////////////////////////////////////////////////////
////
////   COMPUTES AND PRINTS DISTRIBUTION OF TYPES
////
////////////////////////////////////////////////////

void Filter::print_distribution(std::string filename)
{
    std::string distribution_name(filename);
    distribution_name.append(".distribution");
    std::ofstream distribution_file;
    distribution_file.open(distribution_name.c_str());
    
    // OUTPUT INFORMATION ABOUT SOURCE OF DISTRIBUTION
                    distribution_file << "#\tDISTRIBUTION CREATED FROM ";
    if(g_switch==1) distribution_file << "PERTURBATIONS OF ";
                    distribution_file << filename;
    if(g_switch==1) distribution_file << ", USING " << perturbation_samples << " PERTURBATIONS WITH MAGNITUDE " << perturbation_size;
                    distribution_file << '\n';
    
    if(f_switch)    distribution_file << "#\tTypes 1 through " << max_file_filter_type << " obtained from filter in file " << filename_filter << '\n';
    
    if(max_filter_type > max_file_filter_type && (d_switch || df_switch))
    {
       if(f_switch) distribution_file << "#\tTypes " << max_file_filter_type+1 << " through " << max_filter_type << " obtained from other types in data\n";
       else         distribution_file << "#\tTypes " << max_file_filter_type+1 << " through " << max_filter_type << " obtained from types in data\n";
    }
    
                    distribution_file << "#\tColumns indicate: type, Weinberg vector, and count\n";
    

    // SORT BY COUNT, THEN BY WVECTOR
    std::vector<std::map<std::vector<int>,FilterEntry>::iterator> sorted_entries;
    for (auto it = entries.begin(); it != entries.end(); ++it)
        sorted_entries.push_back(it);
    std::sort(sorted_entries.begin(), sorted_entries.end(), compare_by_count_wvector);

    // OUTPUT THE SORTED DISTRIBUTION
    for(auto it = sorted_entries.begin(); it != sorted_entries.end(); ++it) if((*it)->second.count > 0)
    {
        distribution_file << (*it)->second.type << '\t';
        distribution_file << '(';
        for(int i=0; i<(*it)->first.size()-1; i++)
            distribution_file << (*it)->first[i] << ',';
        distribution_file << (*it)->first.back() << ')' << '\t';
        distribution_file << (*it)->second.count << '\n';
    }
}



////////////////////////////////////////////////////
////
////   FINDS FILTER NUMBER ASSOCIATED WITH A CELL
////
////////////////////////////////////////////////////

int Filter::wvector_type(std::vector<int> wvector)
{
    std::map<std::vector<int>,FilterEntry>::iterator it = entries.find(wvector);
    if (it != entries.end()) return it->second.type;
    else                     return 0;
}



////////////////////////////////////////////////////
////
////   RELABELS TYPES OBTAINED FROM DATA IN
////   DECREASING ORDER OF FREQUENCY
////
////////////////////////////////////////////////////

void Filter::relabel_data_types(void)
{
    std::vector<std::map<std::vector<int>,FilterEntry>::iterator> sorted_entries;
    for (auto it = entries.begin(); it != entries.end(); ++it)
        sorted_entries.push_back(it);
    std::sort(sorted_entries.begin(), sorted_entries.end(), compare_by_count_wvector);
    
    for(auto it = sorted_entries.begin(); it != sorted_entries.end(); ++it)
        if((*it)->second.source == 1)
            (*it)->second.type = ++max_filter_type;
}




