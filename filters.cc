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

////   File: filters.cc


#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <stdexcept>

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

void Filter::load_filter()
{
    std::ifstream filter_file(filename_filter);
    if(!filter_file.is_open())
    {
        throw std::runtime_error("Failed to open input file " + filename_filter);
    }
    
    std::string line;
    max_filter_type_from_file = 0;
    int line_counter          = 0;
    
    while (getline(filter_file, line))
    {
        std::istringstream is(line);
        
        line_counter++;
        if(line[0] == '#') continue;    // IGNORE COMMENTED LINES
        else if(line[0] == '*')         // ADD STRUCTURE TYPES
        {
            int index;
            std::string name;
            
            is.ignore();    // IGNORES *

            is >> index;    // READ STRUCTURE TYPE INDEX
            if(index != ++max_filter_type_from_file)
            {
                throw std::runtime_error("Invalid filter provided. Structure types in filter must be ordered consecutively beginning from 1.");
            }
            
            // READ STRUCTURE TYPE NAME UP TO NEXT TAB
            is >> std::ws;  // IGNORE WHITESPACE
            char ch;
            while((ch=is.get())!='\t' && !is.eof())
                name += ch;
            
            // IF NUMBERS FOLLOW, THEN THIS STRUCTURE TYPE IS INDETERMINATE, AND
            // CAN BE RESOLVED USING THE -r OPTION.
            is >> std::ws;           // IGNORE WHITESPACE
            if(isdigit(is.peek()))
            {
                int resolved_primary   = 0;
                int resolved_secondary = 0;
                is >> resolved_primary;
                is >> resolved_secondary;
            }
            continue;
        }
        
        else                         // READ WEINBERG CODES
        {
            int structure_type;
            is >> structure_type;
            
            if(structure_type < 1 || structure_type > max_filter_type_from_file)
            {
                throw std::runtime_error("Invalid filter provided. Structure type " + std::to_string(structure_type) + " on line " + std::to_string(line_counter) + " does not match any of those specified in filter header.");
            }
            
            std::vector<int> new_vector;
            is.ignore(256,'(');
            int next_label;
            while (is >> next_label)
            {
                new_vector.push_back(next_label);
                if (is.peek() == ',') is.ignore();
                if (is.peek() == ')') break;  // IGNORES DATA AFTER FINAL ')'
            }
            
            auto it = entries.find(new_vector);
            if (it != entries.end()) it->second.structure_types.push_back(structure_type);
            else entries.insert({std::move(new_vector), FilterEntry(structure_type)});
        }
    }

    if(r_switch && (count_indeterminate_types()==0))
        std::cout << "Warning: -r switch chosen, but no indeterminate structure types found in filter\n";
}



bool compare_by_count_vector(std::map<std::vector<int>,FilterEntry>::iterator a, std::map<std::vector<int>,FilterEntry>::iterator b)
{
    if(a->second.total == b->second.total)
        return a->first < b->first;
    else return a->second.total > b->second.total;
}


void Filter::increment_or_add(std::vector<int>& topology_vector, int chirality, int count)
{
    auto it = entries.find(topology_vector);
    if (it != entries.end()) { it->second.counts[chirality+1] += count; it->second.total += count; }
    else entries.insert({std::move(topology_vector), FilterEntry(chirality, count)});
}



////////////////////////////////////////////////////
////
////   FINDS FILTER NUMBER ASSOCIATED WITH A CELL
////
////////////////////////////////////////////////////

int Filter::vt_structure_type(std::vector<int> topology_vector)
{
    auto it = entries.find(topology_vector);
    if (it != entries.end()) return it->second.structure_types[0]; 
    else                     return 0;
}



/////////////////////////////////////////////////////
////
////   COMPUTES AND PRINTS CHIRAL DISTRIBUTION OF TYPES
////
/////////////////////////////////////////////////////

void Filter::print_distribution(std::string filename)
{
    std::string distribution_name(filename);
    distribution_name += ".distribution";
    std::ofstream distribution_file;
    distribution_file.open(distribution_name.c_str());
    if (!distribution_file.is_open())
    {
        throw std::runtime_error("Failed to open distribution file: " + distribution_name);
    }
    
    // OUTPUT INFORMATION ABOUT SOURCE OF DISTRIBUTION
    distribution_file << "#\tDISTRIBUTION CREATED FROM ";
    if(g_switch==1) distribution_file << "PERTURBATIONS OF ";
    distribution_file << filename;
    if(g_switch==1) distribution_file << ", USING " << perturbation_samples << " PERTURBATIONS WITH MAGNITUDE " << perturbation_size;
    distribution_file << '\n';

    distribution_file << "#\tColumns indicate: Voronoi topology vector, total count, left-handed, non-chiral, and right-handed types\n";
    
    // SORT BY COUNT, THEN BY VORONOI TOPOLOGY VECTOR
    std::vector<std::map<std::vector<int>,FilterEntry>::iterator> sorted_entries;
    for (auto it = entries.begin(); it != entries.end(); ++it)
        sorted_entries.push_back(it);
    std::sort(sorted_entries.begin(), sorted_entries.end(), compare_by_count_vector);
    
    // OUTPUT THE SORTED DISTRIBUTION
    for(auto it = sorted_entries.begin(); it != sorted_entries.end(); ++it) 
    {
        distribution_file << '(';
        const auto& topology_vector = (*it)->first;
        for(auto value = topology_vector.begin(); value != topology_vector.end() - 1; ++value)
            distribution_file << *value << ',';
        distribution_file << topology_vector.back() << ')' << '\t';
        distribution_file << (*it)->second.total     << '\t';
        distribution_file << (*it)->second.counts[0] << '\t';
        distribution_file << (*it)->second.counts[1] << '\t';
        distribution_file << (*it)->second.counts[2] << '\n';
    }
}



int Filter::count_indeterminate_types()
{
    int count = 0;
    for(auto it = entries.begin(); it != entries.end(); ++it) if(it->second.structure_types.size() > 1) count++;
    return count;
}


// WE WANT TO COPY DATA FROM to_copy FILTER TO THIS FILTER
void Filter::copy_filter(const Filter& to_copy)
{
    for (const auto& entry : to_copy.entries)
    {
        auto it2 =  entries.find(entry.first);
        if  (it2 != entries.end())
        {
            it2->second.counts[0] += entry.second.counts[0];
            it2->second.counts[1] += entry.second.counts[1];
            it2->second.counts[2] += entry.second.counts[2];
            it2->second.total     += entry.second.total;
        }
        else
        {
            entries.emplace(entry.first, FilterEntry(entry.second.counts[0], entry.second.counts[1], entry.second.counts[2]));
        }
    }
}


