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

////   File: output.cc


#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "filters.hh"
#include "variables.hh"


void output_system(std::string filename, Filter &filter)
{
    std::string output_file_name(filename);
    
    if(file_type==1) output_file_name.append(".dump");
    if(file_type==2) output_file_name.append(".cfg");
    std::ofstream output_file(output_file_name.c_str(),std::ofstream::out);
    
    std::string full_line;
    data_file.seekg (0, data_file.beg);
    
    
    ////////////////////////////////////////////////////
    ////
    ////    PRINT HEADER INFO TO FILE
    ////
    ////////////////////////////////////////////////////

    if(file_type==1)    // LAMMPS DUMP FILE
    {
        for(int c=0; c<header_lines; c++)
        {
            getline(data_file, full_line);
            output_file << full_line;
            if(full_line.find("ITEM: ATOMS") == 0)
            {
                output_file << "\tvt ";
                
                if (r_switch) output_file << "resolved_type ";
                if (c_switch)
                {
                    output_file << "cluster_index ";
                    output_file << "cluster_size ";
                }
                output_file << '\n';
            }
            else
                output_file << '\n';
        }
    }
    
    if(file_type==2)    // ATOMEYE CFG FILE
    {
        int entry_count   = 0;
        int extra_entries = 1;          // FOR VORONOI TOPOLOGY
        
        if(r_switch) extra_entries+=1;  // FOR RESOLVED TYPE
        if(c_switch) extra_entries+=2;  // FOR CLUSTER ID AND COUNT
        
        for(int c=0; c<header_lines; c++)
        {
            getline(data_file, full_line);
            
            if (full_line.find("entry_count") == 0)
            {
                std::stringstream s(full_line.substr(14));
                s >> entry_count;
                
                output_file << "entry_count = " << entry_count+extra_entries << "\n";
            }
            
            else output_file << full_line << '\n';
        }
        
        int auxiliary_counter = entry_count-3;
        
        // VORONOI ANALYSIS ATTRIBUTES
        output_file << "auxiliary[" << auxiliary_counter++ << "] = vt\n";
        if (r_switch) output_file << "auxiliary[" << auxiliary_counter++ << "] = resolved_type\n";
        if (c_switch)
        {
            output_file << "auxiliary[" << auxiliary_counter++ << "] = cluster_index\n";
            output_file << "auxiliary[" << auxiliary_counter++ << "] = cluster_size\n";
        }
        
        getline(data_file, full_line); output_file << full_line << '\n';
        getline(data_file, full_line); output_file << full_line << '\n';
    }
    
    
    ////////////////////////////////////////////////////
    ////
    ////    PRINT PARTICLE DATA TO FILE
    ////
    ////////////////////////////////////////////////////
    
    for(int c=0; c<number_of_particles; c++)
    {
        // OUTPUT INITIAL DATA
        getline(data_file, full_line);
        
        // THESE ARE LINES THAT HAVE EITHER ONLY THE ATOMIC MASS
        // OR ELSE ONLY THE CHEMICAL SYMBOL ON THEM
        if(file_type==2)
        if(std::count(full_line.begin(), full_line.end(), '.') < 2)
        {
            output_file << full_line << '\n';
            c--;
            continue;
        }

        output_file << full_line << '\t';
        
        // OUTPUT VORONOI TOPOLOGY
        output_file << vt_structure_types[c] << '\t';
        
        // OUTPUT RESOLVED TYPE AND CLUSTER INFORMATION, AS SPECIFIED
        if(r_switch) output_file << resolved_types[c] << '\t';
        if(c_switch)
        {
            output_file << cluster_index[c] << '\t';
            output_file << cluster_sizes[c] << '\t';
        }
        
        output_file << '\n';
    }
}



////////////////////////////////////////////////////
////
////   PRINT WEINBERG VECTOR FOR EACH PARTICLE
////
////////////////////////////////////////////////////

int print_wvectors(std::string filename)
{
    std::string wvectors_name(filename);
    wvectors_name.append(".wvectors");
    std::ofstream wvector_file(wvectors_name.c_str(), std::ofstream::out);
    
    for(int c=0; c<number_of_particles; c++)
    {
        int  length = all_wvectors[c].back();                          // LENGTH OF WVECTOR
        int plength = all_wvectors[c].end()[-2];//back()-1;            // LENGTH OF PVECTOR
        
        wvector_file << all_wvectors[c][length+plength-3] << '\t';     // NUMBER OF FACES
        
        wvector_file << '(';                                           // P VECTOR
        for(int d=length; d<length+plength-4; d++)
            wvector_file << all_wvectors[c][d] << ',';
        wvector_file << all_wvectors[c][length+plength-4] << ')' << '\t';
        
        wvector_file << '(';                                            // WEINBERG VECTOR
        for(int d=0; d<length; d++)
            wvector_file << all_wvectors[c][d] << ',';
        wvector_file << 1 << ')' << '\t';
        
        wvector_file << all_wvectors[c][length+plength-2] << '\t';     // SYMMETRIES
        wvector_file << all_wvectors[c][length+plength-1] << '\t';     // CHIRALITY
        //wvector_file << all_wvectors[c][length+plength]   << '\t';     // STABLE
        
        wvector_file << '\n';
    }
    
    return 0;
}



