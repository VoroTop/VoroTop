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
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "filters.hh"
#include "variables.hh"


////////////////////////////////////////////////////
////
////   PARSE HEADER OF LAMMPS DUMP FILE.
////
////   LAMMPS dump files contain ITEM: sections in
////   a fixed order.  Optional ITEM: UNITS and
////   ITEM: TIME sections may precede ITEM: TIMESTEP.
////
////////////////////////////////////////////////////

void parse_header(std::ifstream& fp)
{
    if (!fp.is_open())
        throw std::runtime_error("Error opening file");

    std::string line;
    header_lines = 0;

    bool found_atoms = false;

    while(!found_atoms && getline(fp, line))
    {
        header_lines++;

        // SKIP OPTIONAL ITEM: UNITS (appears once at start of file)
        if (line.find("ITEM: UNITS") == 0)
        {
            getline(fp, line);
            header_lines++;
            continue;
        }

        // SKIP OPTIONAL ITEM: TIME (distinct from ITEM: TIMESTEP)
        if (line.find("ITEM: TIME") == 0 && line.find("ITEM: TIMESTEP") != 0)
        {
            getline(fp, line);
            header_lines++;
            continue;
        }

        // SKIP ITEM: TIMESTEP AND ITS VALUE
        if (line.find("ITEM: TIMESTEP") == 0)
        {
            getline(fp, line);
            header_lines++;
            continue;
        }

        // ITEM: NUMBER OF ATOMS
        if (line.find("ITEM: NUMBER OF ATOMS") == 0)
        {
            getline(fp, line);
            header_lines++;
            number_of_particles = std::stoi(line);
            continue;
        }

        // ITEM: BOX BOUNDS
        if (line.find("ITEM: BOX BOUNDS") == 0)
        {
            // DETECT GENERAL TRICLINIC (abc origin ...)
            if (line.find("abc origin") != std::string::npos)
            {
                std::cerr << "VoroTop does not support general triclinic boxes." << std::endl;
                std::cerr << "Please convert to an orthogonal box." << std::endl;
                exit(1);
            }

            // DETECT RESTRICTED TRICLINIC (xy xz yz ...)
            if (line.find("xy xz yz") != std::string::npos)
            {
                triclinic_crystal_system = 1;
                std::cerr << "VoroTop does not currently support triclinic systems." << std::endl;
                std::cerr << "Please convert to an orthogonal box." << std::endl;
                exit(1);
            }

            // ORTHOGONAL BOX: READ xlo xhi, ylo yhi, zlo zhi
            triclinic_crystal_system = 0;

            getline(fp, line); header_lines++;
            std::istringstream(line) >> xlo >> xhi;

            getline(fp, line); header_lines++;
            std::istringstream(line) >> ylo >> yhi;

            getline(fp, line); header_lines++;
            std::istringstream(line) >> zlo >> zhi;

            continue;
        }

        // ITEM: ATOMS — DETERMINE COLUMN LAYOUT
        if (line.find("ITEM: ATOMS") == 0)
        {
            std::istringstream iss(line);
            std::string entry;
            iss >> entry;   // "ITEM:"
            iss >> entry;   // "ATOMS"

            std::vector<std::string> attribute_labels;
            while(iss >> entry)
                attribute_labels.push_back(entry);
            particle_attributes = attribute_labels.size();

            index_id = -1;
            index_x  = -1;
            index_y  = -1;
            index_z  = -1;
            scaled_coordinates = 0;

            for(int c = 0; c < (int)attribute_labels.size(); c++)
            {
                if(attribute_labels[c] == "id") index_id = c;

                // ABSOLUTE COORDINATES (WRAPPED OR UNWRAPPED)
                if(attribute_labels[c] == "x"  || attribute_labels[c] == "xu")  index_x = c;
                if(attribute_labels[c] == "y"  || attribute_labels[c] == "yu")  index_y = c;
                if(attribute_labels[c] == "z"  || attribute_labels[c] == "zu")  index_z = c;

                // SCALED/FRACTIONAL COORDINATES (WRAPPED OR UNWRAPPED)
                if(attribute_labels[c] == "xs" || attribute_labels[c] == "xsu") { index_x = c; scaled_coordinates = 1; }
                if(attribute_labels[c] == "ys" || attribute_labels[c] == "ysu") { index_y = c; scaled_coordinates = 1; }
                if(attribute_labels[c] == "zs" || attribute_labels[c] == "zsu") { index_z = c; scaled_coordinates = 1; }
            }

            if(index_x == -1 || index_y == -1)
                throw std::runtime_error("Insufficient coordinate columns in dump file.");

            if(index_z == -1) dimension = 2;

            found_atoms = true;
        }
    }

    if(!found_atoms)
        throw std::runtime_error("ITEM: ATOMS not found in dump file.");

    if(index_id == -1)
        std::cerr << "Warning: No particle IDs found in file. Using default IDs." << std::endl;

    if(dimension != 2 && dimension != 3)
        throw std::runtime_error("Unsupported dimension: " + std::to_string(dimension));

    // COMPUTE VORO++ GRID BLOCK SIZES.
    // EACH BLOCK SHOULD BE ROUGHLY A SQUARE OR CUBE AND CONTAIN
    // ROUGHLY THE SAME NUMBER OF PARTICLES.
    if(dimension == 2)
    {
        int total_blocks = number_of_particles / 4. + 1;
        double area = (xhi - xlo) * (yhi - ylo);
        double length_per_block = pow(area / (double)total_blocks, 1./2.);
        n_x = (int)((xhi - xlo) / length_per_block) + 1;
        n_y = (int)((yhi - ylo) / length_per_block) + 1;
    }

    if(dimension == 3)
    {
        int total_blocks = number_of_particles / 4. + 1;
        double volume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo);
        double length_per_block = pow(volume / (double)total_blocks, 1./3.);
        n_x = (int)((xhi - xlo) / length_per_block) + 1;
        n_y = (int)((yhi - ylo) / length_per_block) + 1;
        n_z = (int)((zhi - zlo) / length_per_block) + 1;
    }

    if(particles_in_eps == 0) particles_in_eps = number_of_particles;
    if(particles_in_eps > number_of_particles)
    {
        std::cout << "Warning: particles_in_eps exceeds number of particles. "
                  << "Drawing all " << number_of_particles << "." << std::endl;
        particles_in_eps = number_of_particles;
    }
}


////////////////////////////////////////////////////
////
////   IMPORT PARTICLE DATA FROM FILE.
////   HANDLES ABSOLUTE, SCALED, AND UNWRAPPED
////   COORDINATE TYPES.
////
////////////////////////////////////////////////////

void import_data()
{
    FILE *in_file = fopen(filename_data.c_str(), "r");
    if (!in_file)
        throw std::runtime_error("Error opening file: " + filename_data);

    for(int i = 0; i < header_lines; i++)
        fscanf(in_file, "%*[^\n]\n");

    double Lx = xhi - xlo;
    double Ly = yhi - ylo;
    double Lz = zhi - zlo;

    for(int c = 0; c < number_of_particles; c++)
    {
        int id = c + 1;
        double x = 0, y = 0, z = 0;
        double junk;

        for(int d = 0; d < particle_attributes; d++)
        {
            int result;
            if     (d == index_id) result = fscanf(in_file, "%d",  &id);
            else if(d == index_x)  result = fscanf(in_file, "%lg", &x);
            else if(d == index_y)  result = fscanf(in_file, "%lg", &y);
            else if(d == index_z)  result = fscanf(in_file, "%lg", &z);
            else                   result = fscanf(in_file, "%lg", &junk);

            if (result != 1)
                throw std::runtime_error("Error reading particle data from file");
        }

        // CONVERT SCALED/FRACTIONAL COORDINATES TO ABSOLUTE
        if(scaled_coordinates)
        {
            x = xlo + x * Lx;
            y = ylo + y * Ly;
            if(dimension == 3) z = zlo + z * Lz;
        }

        // WRAP COORDINATES INTO PERIODIC BOX.
        // USES FMOD TO HANDLE UNWRAPPED COORDINATES THAT MAY BE
        // MULTIPLE BOX LENGTHS OUTSIDE THE BOX.
        x = xlo + fmod(x - xlo, Lx);
        if(x < xlo) x += Lx;

        y = ylo + fmod(y - ylo, Ly);
        if(y < ylo) y += Ly;

        if(dimension == 3)
        {
            z = zlo + fmod(z - zlo, Lz);
            if(z < zlo) z += Lz;
        }

        // STORE COORDINATES
        particle_coordinates[dimension * c]     = x;
        particle_coordinates[dimension * c + 1] = y;
        if(dimension == 3)
            particle_coordinates[dimension * c + 2] = z;

        // STORE PARTICLE ID
        if(index_id != -1) particle_ids[c] = id;
        else               particle_ids[c] = c + 1;
    }

    fclose(in_file);
}
