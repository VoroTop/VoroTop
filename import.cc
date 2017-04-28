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

////   File: import.cc


#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "import.hh"
#include "voro++.hh"
#include "filters.hh"
#include "vorotop.hh"
#include "functions.hh"

using namespace voro;



////////////////////////////////////////////////////
////
////   PARSE HEADER OF INPUT FILE.  CURRENTLY RECOGNIZES
////   LAMMPS DUMP AND ATOMEYE EXTENDED CFG FORMATS.
////
////////////////////////////////////////////////////

int parse_header(std::ifstream& fp)
{
    ////////////////////////////////////////////////////
    ////
    ////    INITIALIZE VARIABLES
    ////
    ////////////////////////////////////////////////////
    
    file_type=0;    // 1 LAMMPS DUMP FORMAT
                    // 2 AtomEye cfg file ("Extended" format)
    
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
    
    std::string full_line;
    getline(fp, full_line);
    
    if (full_line.find("ITEM: TIMESTEP") != std::string::npos)
    {
        file_type = 1;    // LAMMPS DUMP
        fp >> timestep;   // RECORD timestep
    }
    
    else if (full_line.find("Number of particles") != std::string::npos)
    {
        file_type = 2;    // ATOMEYE CFG
        std::stringstream s(full_line.substr(21));
        s >> number_of_particles;   // RECORD number_of_particles
        scaled_coordinates=1;   // ATOMEYE CFG IS ALWAYS SCALED COORDINATES
    }
    
    else
    {
        std::cout << "Unrecognized data file format\n";
        exit(0);
    }

    
    ////////////////////////////////////////////////////
    ////
    ////    PARSE FILE
    ////
    ////////////////////////////////////////////////////
    
    if(file_type == 1)  // LAMMPS DUMP FILE
    {
        bool done=0;
        double hi_bound[3];
        
        while(done==0)
        {
            getline(fp, full_line);
            std::istringstream iss(full_line);
            
            if (full_line.find("ITEM: NUMBER OF ATOMS") != std::string::npos)
            {
                fp >> number_of_particles;
            }
            
            else if (full_line.find("ITEM: BOX BOUNDS pp pp pp") != std::string::npos)          // PERIODIC IN ALL THREE DIRECTIONS
            {
                for(int c=0; c<3; c++)
                {
                    fp >> origin[c];
                    fp >> hi_bound[c];
                    supercell_edges[c][c] = hi_bound[c]-origin[c];
                }
            }
            
            else if (full_line.find("ITEM: BOX BOUNDS xy xz yz pp pp pp") != std::string::npos) // TRICLINIC CRYSTAL SYSTEM
            {
                double xlo_bound, xhi_bound, xy;
                double ylo_bound, yhi_bound, xz;
                double zlo_bound, zhi_bound, yz;
                
                fp >> xlo_bound; fp >> xhi_bound; fp >> xy;
                fp >> ylo_bound; fp >> yhi_bound; fp >> xz;
                fp >> zlo_bound; fp >> zhi_bound; fp >> yz;
                
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
            }
            
            else if (full_line.find("ITEM: ATOMS ") != std::string::npos)
            {
                // DETERMINE ATOM ATTRIBUTES, INDEX OF COORDINATES, AND SCALING
                std::string entry;
                iss >> entry;   // READS IN "ITEM:"
                iss >> entry;   // READS IN "ATOMS"
                
                while(iss >> entry)                             // READ ATTRIBUTE LABELS, I.E., 'id', 'type', 'x', 'y', and 'z'
                    attribute_labels.push_back(entry);
                
                xindex=-1;
                for(int c=0; c<attribute_labels.size(); c++)
                {
                    if(attribute_labels[c]=="x" || attribute_labels[c]=="xs")
                    {
                        if(attribute_labels[c]=="x") scaled_coordinates=0;
                        else                         scaled_coordinates=1;

                        attribute_labels.erase(attribute_labels.begin()+c, attribute_labels.begin()+c+3);
                        xindex=c;
                    }
                }
                
                if(xindex==-1)
                {
                    std::cout << "Insufficient xyz coordinates included in dump file.\n";
                    exit(-1);
                }
                
                done=1;
            }
        }
        
        cfg_lscale = 1.;
        cfg_atomic_mass = 1.;   // WE SHOULD FIX, TO ALLOW FOR INCLUSION OF MASS
        cfg_chem_symbol = "X";  // WE SHOULD FIX, TO ALLOW FOR INCLUSION OF ELEMENT
    }
    
    else if(file_type == 2) // AtomEye "Extended" format
    {
        int entries=0;
        bool done  =0;
        while(done==0)
        {
            // BY DEFAULT THIS IGNORE BLANK LINES, OR LINES BEGINNING WITH # (OR ANY LINE NOT CONTAINING ANY OF THE BELOW)
            getline(fp, full_line);
            
            // READS IN SUPERCELL EDGES
            if (full_line.find("H0") != std::string::npos)
            {
                int i,j;
                double value;
                std::string stuff;
                
                std::stringstream s(full_line.substr(3));
                s >> i; s.ignore();     // IGNORES ','
                s >> j; s.ignore();     // IGNORES ')'
                s >> stuff;             // IGNORES " = "
                s >> value;
                
                supercell_edges[i-1][j-1]=value;
            }
            
            if (full_line.find("A = ") != std::string::npos)
            {
                std::stringstream s(full_line.substr(4));
                s >> cfg_lscale;
            }
            else if (full_line.find(".NO_VELOCITY.") != std::string::npos)
            {
                no_velocity = 1;
            }
            else if (full_line.find("entry_count = ") != std::string::npos)
            {
                std::stringstream s(full_line.substr(14));
                s >> entries;// cfg_entry_count;
                done = 1;   // ONCE WE HIT entry_count, WE STOP READING IN auxiliary NAMES
            }
            else
            {
                // THERE ARE OTHER LINES THAT WE CURRENTLY IGNORE
            }
        }

        if(no_velocity==1) entries -= 3;
        else               entries -= 6;
        
        // entries DETERMINES HOW MANY attribute_labels TO READ IN
        for(int i=0; i<entries; i++)
        {
            getline(fp, full_line);
            std::stringstream s(full_line.substr(15));  // THIS PLACES US AFTER auxiliary[X]
            
            std::string entry;
            s >> entry;
            attribute_labels.push_back(entry);
        }
    }
    
    int total_blocks = number_of_particles/8;
    double volume = supercell_edges[0][0]*supercell_edges[1][1]*supercell_edges[2][2];  // APPROXIMATION
    double volume_per_block = volume/(double)total_blocks;
    double length_per_cube = pow(volume_per_block, 1./3.);
    
    n_x = int(supercell_edges[0][0]/length_per_cube+1);
    n_y = int(supercell_edges[1][1]/length_per_cube+1);
    n_z = int(supercell_edges[2][2]/length_per_cube+1);
    
    return file_type;
}



////////////////////////////////////////////////////
////
////   IMPORTS DATA POINTS FROM DUMP FILE
////
////////////////////////////////////////////////////

void import_dump_file(std::ifstream& fp, particle_order &vo, container_periodic &con)
{
    double x,y,z;
    
    // STORES ALL PARTICLE DATA TO particle_data[][]
    for(int c=0; c<number_of_particles; c++)
    {
        for(int d=0; d<xindex; d++)
            fp >> particle_data[c][d];
        
        fp >> xcoord[c]; x = xcoord[c];
        fp >> ycoord[c]; y = ycoord[c];
        fp >> zcoord[c]; z = zcoord[c];
        
        for(int d=xindex+3; d<attribute_labels.size()+3; d++)
            fp >> particle_data[c][d-3];
        
        // ADJUST COORDINATES SO SYSTEM UNSCALED AND AT ORIGIN
        if(scaled_coordinates==0)
        {
            x -= origin[0];
            y -= origin[1];
            z -= origin[2];
        }
        
        if(scaled_coordinates==1)
        {
            double newx = supercell_edges[0][0]*x + supercell_edges[1][0]*y + supercell_edges[2][0]*z;
            double newy = supercell_edges[0][1]*x + supercell_edges[1][1]*y + supercell_edges[2][1]*z;
            double newz = supercell_edges[0][2]*x + supercell_edges[1][2]*y + supercell_edges[2][2]*z;
            
            x=newx;
            y=newy;
            z=newz;
        }
        
        con.put(vo,c,x,y,z);
    }
}



void import_atomeye_file(std::ifstream& fp, particle_order &vo, container_periodic &con)
{
    double x,y,z;
    std::string full_line;

    for(int c=0; c<number_of_particles; c++)
    {
        getline(fp, full_line);
        
        if(countWordsInString(full_line)==1)                // ATOMIC MASS AND CHEM SYMBOL
        {
            cfg_atomic_mass = ::atof(full_line.c_str());
            fp >> cfg_chem_symbol; getline(fp, full_line);
            c--;  continue;
        }

        std::stringstream ss(full_line);

        ss >> xcoord[c];
        ss >> ycoord[c];
        ss >> zcoord[c];
        
        for(int d=0; d<attribute_labels.size(); d++)
            ss >> particle_data[c][d];
        
        x = xcoord[c];
        y = ycoord[c];
        z = zcoord[c];
        
        if(x<0) x+=1.; if(x>1) x-=1.;
        if(y<0) y+=1.; if(y>1) y-=1.;
        if(z<0) z+=1.; if(z>1) z-=1.;
        
        double newx = supercell_edges[0][0]*x + supercell_edges[1][0]*y + supercell_edges[2][0]*z;
        double newy = supercell_edges[0][1]*x + supercell_edges[1][1]*y + supercell_edges[2][1]*z;
        double newz = supercell_edges[0][2]*x + supercell_edges[1][2]*y + supercell_edges[2][2]*z;
        
        con.put(vo,c,newx,newy,newz);
    }
}





