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

////   File: import.cc


#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <stdio.h>
#include <stdexcept>

#include "filters.hh"
#include "variables.hh"


// TAKEN FROM https://stackoverflow.com/questions/41045236
// TO FACILITATE IGNORING DATA IN FILE
template <typename CharT>
std::basic_istream<CharT>& ignore_m(std::basic_istream<CharT>& in){
    std::string ignoredValue;
    return in >> ignoredValue;
}


////////////////////////////////////////////////////
////
////   PARSE HEADER OF INPUT FILE.  CURRENTLY RECOGNIZES
////   LAMMPS DUMP AND ATOMEYE EXTENDED CFG FORMATS.
////
////////////////////////////////////////////////////

const int LAMMPS_HEADER_LINES = 2;
const int BOX_BOUNDS_LINES = 4;

void parse_header(std::ifstream& fp)
{
    ////////////////////////////////////////////////////
    ////
    ////    INITIALIZE VARIABLES
    ////
    ////////////////////////////////////////////////////
    
    origin[0]=0.;
    origin[1]=0.;
    origin[2]=0.;
    for(int c=0; c<3; c++)
        for(int d=0; d<3; d++)
            supercell_edges[c][d]=0;
    
    
    ////////////////////////////////////////////////////
    ////
    ////    DETERMINE FILE TYPE
    ////
    ////////////////////////////////////////////////////
    
    if (!fp.is_open()) {
        throw std::runtime_error("Error opening file");
    }

    std::string full_line;
    getline(fp, full_line);
    
    if (full_line.find("ITEM: TIMESTEP") != std::string::npos)
    {
        // LAMMPS DUMP FILE FORMAT
        header_lines = LAMMPS_HEADER_LINES;
    }
    
    else
    {
        throw std::runtime_error("Unrecognized data file format");
    }
    
    
    ////////////////////////////////////////////////////
    ////
    ////    PARSE FILE
    ////
    ////////////////////////////////////////////////////
    
    bool done=0;
    while(done==0)
    {
        getline(fp, full_line);
        std::istringstream iss(full_line);
        
        if (full_line.find("ITEM: NUMBER OF ATOMS") != std::string::npos)
        {
            fp >> number_of_particles;
            header_lines+=2;
        }
        
        if (full_line.find("ITEM: DIMENSION") != std::string::npos)
        {
            fp >> dimension;
            header_lines+=2;
        }
        
        else if (full_line.find("ITEM: BOX BOUNDS pp pp pp") != std::string::npos)          // PERIODIC IN ALL THREE DIRECTIONS
        {
            for(int c=0; c<3; c++)
            {
                fp >> origin[c];
                fp >> hi_bound[c];
                supercell_edges[c][c] = hi_bound[c]-origin[c];
            }
            header_lines += BOX_BOUNDS_LINES;
        }
        
        else if (full_line.find("ITEM: BOX BOUNDS xy xz yz pp pp pp") != std::string::npos) // TRICLINIC CRYSTAL SYSTEM
        {
            fp >> xlo_bound >> xhi_bound >> xy;
            fp >> ylo_bound >> yhi_bound >> xz;
            fp >> zlo_bound >> zhi_bound >> yz;
            
            origin  [2] = zlo_bound;
            hi_bound[2] = zhi_bound;
            
            origin  [1] = ylo_bound - fmin(0.,yz);
            hi_bound[1] = yhi_bound - fmax(0.,yz);
            
            origin  [0] = xlo_bound - fmin(fmin(0.0,xy),fmin(xz,xy+xz));
            hi_bound[0] = xhi_bound - fmax(fmax(0.0,xy),fmax(xz,xy+xz));
            
            supercell_edges[1][0] = xy;
            supercell_edges[2][0] = xz;
            supercell_edges[2][1] = yz;
            
            for(int c=0; c<3; c++)
                supercell_edges[c][c] = hi_bound[c]-origin[c];
            
            header_lines+=4;
        }
        
        else if (full_line.find("ITEM: ATOMS ") != std::string::npos)
        {
            // DETERMINE ATOM ATTRIBUTES, INDEX OF COORDINATES, AND SCALING
            std::string entry;
            iss >> entry;   // READS IN "ITEM:"
            iss >> entry;   // READS IN "ATOMS"
            
            std::vector<std::string> attribute_labels;
            while(iss >> entry)                             // READ ATTRIBUTE LABELS, I.E., 'id', 'type', 'x', 'y', and 'z'
                attribute_labels.push_back(entry);
            particle_attributes = attribute_labels.size();
            
            index_id=-1;
            index_x=-1;
            index_y=-1;
            index_z=-1;
            scaled_coordinates=0;
            
            for(long unsigned int c=0; c<attribute_labels.size(); c++)
            {
                if(attribute_labels[c]=="id") index_id = c;

                if(attribute_labels[c]=="x") index_x = c;
                if(attribute_labels[c]=="y") index_y = c;
                if(attribute_labels[c]=="z") index_z = c;

                if(attribute_labels[c]=="xs") { index_x = c; scaled_coordinates=1; }
                if(attribute_labels[c]=="ys") { index_y = c; scaled_coordinates=1; }
                if(attribute_labels[c]=="zs") { index_z = c; scaled_coordinates=1; }
            }
            
            if(index_z==-1) dimension = 2;
            if(index_x==-1 || index_y==-1)
            {
                throw std::runtime_error("Insufficient xyz coordinates included in dump file.");
            }
            
            header_lines+=1;
            done=1;
        }
    }
    
    if (dimension != 2 && dimension != 3) {
        throw std::runtime_error("Unsupported dimension: " + std::to_string(dimension));
    }
    
    // EACH BLOCK SHOULD BE ROUGLY A SQUARE OR CUBE AND CONTAIN 
    // ROUGHLY THE SAME NUMBER OF PARTICLES
    if(dimension==2)
    {
        int total_blocks = number_of_particles/4.+1;
        double area = supercell_edges[0][0]*supercell_edges[1][1];  
        double area_per_block = area/(double)total_blocks;          
        double length_per_cube = pow(area_per_block, 1./2.);
        n_x = (int)supercell_edges[0][0]/length_per_cube+1;
        n_y = (int)supercell_edges[1][1]/length_per_cube+1;
    }

    if(dimension==3)
    {
        int total_blocks = number_of_particles/4.+1;
        double volume = supercell_edges[0][0]*supercell_edges[1][1]*supercell_edges[2][2]; 
        double volume_per_block = volume/(double)total_blocks;
        double length_per_cube = pow(volume_per_block, 1./3.);
        n_x = (int)supercell_edges[0][0]/length_per_cube+1;
        n_y = (int)supercell_edges[1][1]/length_per_cube+1;
        n_z = (int)supercell_edges[2][2]/length_per_cube+1;
    }
}


// THIS USES A C-STYLE FILE* FOR READING DATA FROM THE INPUT FILE.
void import_data()
{
    FILE *in_file = fopen(filename_data.c_str(), "r");
    if (!in_file) {
        throw std::runtime_error("Error opening file");
    }

    for(int i=0; i<header_lines; i++)
        fscanf(in_file, "%*[^\n]\n");

    for(int c=0; c<number_of_particles; c++)
    {
        int id;
        double x,y,z=0;
        double junk;

        for(int d=0; d<particle_attributes; d++)
        {
            int result;
            if     (d==index_id) result = fscanf(in_file,"%d",&id);
            else if(d==index_x)  result = fscanf(in_file,"%lg",&x);
            else if(d==index_y)  result = fscanf(in_file,"%lg",&y);
            else if(d==index_z)  result = fscanf(in_file,"%lg",&z);
            else                 result = fscanf(in_file,"%lg",&junk);

            if (result != 1) {
                throw std::runtime_error("Error reading attribute data from file");
            }
        }

        if(scaled_coordinates==0)
        {
            x -= origin[0];
            y -= origin[1];
            z -= origin[2];
        }
        
        if(scaled_coordinates==1)
        {
            double newx = supercell_edges[0][0]*x + supercell_edges[1][0]*y + supercell_edges[2][0]*z;
            double newy = supercell_edges[0][1]*x + supercell_edges[1][1]*y + supercell_edges[2][0]*z;
            double newz = supercell_edges[0][2]*x + supercell_edges[1][2]*y + supercell_edges[2][0]*z;

            x=newx;
            y=newy;
            z=newz;
        }

        if(dimension==2)
        {
            particle_coordinates[dimension*c]  =x;
            particle_coordinates[dimension*c+1]=y;
        }
        
        if(dimension==3)
        {
            particle_coordinates[dimension*c]  =x;
            particle_coordinates[dimension*c+1]=y;
            particle_coordinates[dimension*c+2]=z;
        }
        
        if(index_id!=-1) particle_ids[c]=id;    // IF PROVIDED, STORE GIVEN PARTICLE IDS
        else particle_ids[c]=c+1;               // IF ID NOT PROVIDED, LABEL EACH PARTICLE BY ITS
                                                // GIVEN ORDER, BEGINNING FROM 1. INTERNALLY, PARTICLES
                                                // WILL USE IDS 0, 1, ..., number_of_particles-1.
    }
    
    fclose(in_file);
}



