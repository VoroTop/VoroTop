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
//   LINES BEGINNING WITH # ARE IGNORED AND ARE USED FOR
//   COMMENTS.
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
#include "vorotop.hh"



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
            
            FilterEntry new_entry;
            new_entry.type = type;

            is.ignore(256,'(');
            int next_label;
            while (is >> next_label)
            {
                new_entry.wvector.push_back(next_label);
                if (is.peek() == ',') is.ignore();
                if (is.peek() == ')') break;  // IGNORES DATA AFTER FINAL ')'
            }
            
            entries.push_back(new_entry);
        }
    }
    
    file_filter_types = entries.size();
    max_filter_type   = max_file_filter_type;
}



void Filter::increment_or_add(std::vector<int> wvector, int n)
{
    // IMPLEMENTS A BINARY SEARCH
    int first = 0;
    int last = size()-1;
    int middle = (first+last)/2;
    
    while( first <= last )
    {
        if ( entries[middle].wvector < wvector )
            first = middle + 1;
        else if ( entries[middle].wvector == wvector )
        {
            entries[middle].count+=n;
            break;
        }
        else
            last = middle - 1;
        
        middle = (first + last)/2;
    }
    
    if ( first > last )
        entries.insert(entries.begin()+first, FilterEntry(++max_filter_type,n,wvector));
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
    
    if(f_switch)
        distribution_file << "#\tTypes 1 through " << max_file_filter_type << " obtained from filter in file " << filename_filter << '\n';
    if(max_filter_type > max_file_filter_type && (d_switch || df_switch))
    {
        if(f_switch) distribution_file << "#\tTypes " << max_file_filter_type+1 << " through " << max_filter_type << " obtained from other types in data\n";
        else         distribution_file << "#\tTypes " << max_file_filter_type+1 << " through " << max_filter_type << " obtained from types in data\n";
    }
    distribution_file << "#\tColumns indicate: type, Weinberg vector, and count\n";
    
    // OUTPUT THE DISTRIBUTION
    for (std::vector<FilterEntry>::iterator it = entries.begin() ; it != entries.end(); ++it) if(it->count>0)
    {
        distribution_file << it->type << '\t';
        distribution_file << '(';
        for(int i=0; i<it->wvector.size()-1; i++)
            distribution_file << it->wvector[i] << ',';
        distribution_file << it->wvector.back() << ')' << '\t';
        distribution_file << it->count << '\n';
    }
}



////////////////////////////////////////////////////
////
////   FINDS FILTER NUMBER ASSOCIATED WITH A CELL
////
////////////////////////////////////////////////////

int Filter::wvector_type(std::vector<int> wvector)
{
    if(wvector.empty()) return -1;
    
    // IMPLEMENTS A BINARY SEARCH
    int first  = 0;
    int last   = file_filter_types-1;//entries.size()-1;
    int middle = (first+last)/2;
    
    while( first <= last )
    {
        if ( entries[middle].wvector < wvector )
            first = middle + 1;
        else if ( entries[middle].wvector == wvector )
        {
            return entries[middle].type;
            break;
        }
        else
            last = middle - 1;
        
        middle = (first + last)/2;
    }
    
    if ( first > last )
        return 0;
    
    return 0;
}



////////////////////////////////////////////////////
////
////   FINDS FILTER NUMBER ASSOCIATED WITH A CELL
////
////////////////////////////////////////////////////

int Filter::wvector_index(std::vector<int> wvector)
{
    if(wvector.empty()) return -1;
    
    // IMPLEMENTS A BINARY SEARCH
    int first = 0;
    int last = entries.size()-1;
    int middle = (first+last)/2;
    
    while( first <= last )
    {
        if ( entries[middle].wvector < wvector )
            first = middle + 1;
        else if ( entries[middle].wvector == wvector )
        {
            return middle;
            break;
        }
        else
            last = middle - 1;
        
        middle = (first + last)/2;
    }
    
    if ( first > last )
        return -1;
    
    return -1;
}



void Filter::sort_dist_by_count(void)
{
    sort_by_type();
    
    sort(entries.begin(), entries.begin()+file_filter_types, compByCount);
    sort(entries.begin()+file_filter_types,   entries.end(), compByCount);
    
    int label = max_file_filter_type;
    for(int c=file_filter_types; c<entries.size(); c++)
        entries[c].type = ++label;
}



// FOR CLUSTERING, WE NEED TO PLACE ALL TYPES FROM THE FILTER FILE
void Filter::sort_for_clustering(void)
{
    sort_by_type();
    sort(entries.begin(), entries.begin()+file_filter_types, compByWvector);
}



int  Filter::get_entry_type (int n)
{
    if(n>0) return entries[n].type;
    else    return 0;
}



int  Filter::get_entry_count(int n)
{
    if(n>0) return entries[n].count;
    else    return 0;
}





