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
#include <cctype>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "filters.hh"
#include "variables.hh"


////////////////////////////////////////////////////
////
////   PARSE LAMMPS DUMP FILE HEADER.
////
////   LAMMPS dump files contain ITEM: sections in
////   a fixed order.  Optional ITEM: UNITS and
////   ITEM: TIME sections may precede ITEM: TIMESTEP.
////
////////////////////////////////////////////////////

static void parse_lammps_dump(std::ifstream& fp)
{
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

            index_id   = -1;
            index_type = -1;
            index_x    = -1;
            index_y    = -1;
            index_z    = -1;
            scaled_coordinates = 0;

            for(int c = 0; c < (int)attribute_labels.size(); c++)
            {
                if(attribute_labels[c] == "id")   index_id   = c;
                if(attribute_labels[c] == "type") index_type = c;

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
}


////////////////////////////////////////////////////
////
////   PARSE LAMMPS DATA FILE HEADER.
////
////   Data files begin with a comment line, followed
////   by header keywords ("N atoms", "xlo xhi", etc.)
////   and then named sections ("Atoms", "Masses", etc.)
////   in any order.
////
////   The column layout of the Atoms section depends
////   on the atom_style.  If a comment is present
////   (e.g. "Atoms # full"), it is used to determine
////   the layout; otherwise, the layout is inferred
////   from the number of columns.
////
////////////////////////////////////////////////////

static void parse_lammps_data(std::ifstream& fp)
{
    std::string line;
    header_lines = 0;
    scaled_coordinates = 0;
    triclinic_crystal_system = 0;

    // LINE 1: TITLE/COMMENT (SKIP)
    getline(fp, line);
    header_lines++;

    // PARSE HEADER KEYWORDS.
    // Each header line has format: value(s) keyword (e.g. "100 atoms").
    // The header ends when a section keyword (starting with a letter) is found.
    while(getline(fp, line))
    {
        header_lines++;

        // SKIP BLANK LINES
        size_t start = line.find_first_not_of(" \t\r");
        if(start == std::string::npos) continue;

        // A LINE STARTING WITH A LETTER IS A SECTION KEYWORD — HEADER IS DONE
        if(std::isalpha(line[start])) break;

        // STRIP TRAILING COMMENTS
        size_t comment = line.find('#');
        std::string content = (comment != std::string::npos) ? line.substr(0, comment) : line;

        // PARSE HEADER KEYWORDS
        if(content.find("atoms") != std::string::npos && content.find("atom types") == std::string::npos)
            std::istringstream(content) >> number_of_particles;
        else if(content.find("xlo xhi") != std::string::npos)
            std::istringstream(content) >> xlo >> xhi;
        else if(content.find("ylo yhi") != std::string::npos)
            std::istringstream(content) >> ylo >> yhi;
        else if(content.find("zlo zhi") != std::string::npos)
            std::istringstream(content) >> zlo >> zhi;
        else if(content.find("xy xz yz") != std::string::npos)
        {
            triclinic_crystal_system = 1;
            std::cerr << "VoroTop does not currently support triclinic systems." << std::endl;
            std::cerr << "Please convert to an orthogonal box." << std::endl;
            exit(1);
        }
    }

    // `line` NOW CONTAINS THE FIRST SECTION KEYWORD.
    // FIND THE "Atoms" SECTION, SKIPPING ANY OTHER SECTIONS.
    bool found_atoms = false;

    while(!found_atoms)
    {
        size_t start = line.find_first_not_of(" \t\r");
        if(start != std::string::npos && line.substr(start, 5) == "Atoms"
           && (start + 5 >= line.size() || !std::isalpha(line[start + 5])))
        {
            found_atoms = true;

            // DETERMINE ATOM_STYLE FROM OPTIONAL COMMENT (e.g. "Atoms # full")
            // style: 0=atomic, 1=charge, 2=molecular/bond/angle, 3=full, -1=auto-detect
            int style = -1;
            size_t hash = line.find('#');
            if(hash != std::string::npos)
            {
                std::string style_str = line.substr(hash + 1);
                size_t s = style_str.find_first_not_of(" \t");
                size_t e = style_str.find_last_not_of(" \t\r\n");
                if(s != std::string::npos && e != std::string::npos)
                    style_str = style_str.substr(s, e - s + 1);

                if     (style_str == "atomic")                                        style = 0;
                else if(style_str == "charge")                                        style = 1;
                else if(style_str == "molecular" || style_str == "bond" ||
                        style_str == "angle")                                         style = 2;
                else if(style_str == "full")                                          style = 3;
                else
                    std::cerr << "Warning: unrecognized atom_style '" << style_str
                              << "', will auto-detect from column count." << std::endl;
            }

            // SKIP MANDATORY BLANK LINE AFTER SECTION KEYWORD
            getline(fp, line);
            header_lines++;

            // COUNT COLUMNS ON FIRST ATOM LINE FOR AUTO-DETECTION
            std::streampos atom_data_pos = fp.tellg();
            getline(fp, line);
            int num_columns = 0;
            {
                std::istringstream iss(line);
                std::string token;
                while(iss >> token) num_columns++;
            }
            fp.seekg(atom_data_pos);  // REWIND TO BEFORE FIRST ATOM LINE

            // SET COLUMN INDICES BASED ON ATOM_STYLE
            index_id = 0;

            if(style >= 0)
            {
                switch(style)
                {
                    case 0: index_type = 1; index_x = 2; index_y = 3; index_z = 4; break;  // atomic:     id type x y z
                    case 1: index_type = 1; index_x = 3; index_y = 4; index_z = 5; break;  // charge:     id type q x y z
                    case 2: index_type = 2; index_x = 3; index_y = 4; index_z = 5; break;  // molecular:  id mol type x y z
                    case 3: index_type = 2; index_x = 4; index_y = 5; index_z = 6; break;  // full:       id mol type q x y z
                }
            }
            else
            {
                // AUTO-DETECT FROM COLUMN COUNT.
                // Without image flags: 5=atomic, 6=charge/molecular, 7=full
                // With image flags (+3): 8=atomic, 9=charge/molecular, 10=full
                int base = (num_columns > 7) ? num_columns - 3 : num_columns;
                switch(base)
                {
                    case 5:  index_type = 1; index_x = 2; index_y = 3; index_z = 4; break;
                    case 6:  index_type = 1; index_x = 3; index_y = 4; index_z = 5; break;
                    case 7:  index_type = 2; index_x = 4; index_y = 5; index_z = 6; break;
                    default:
                        std::cerr << "Warning: cannot determine atom_style from "
                                  << num_columns << " columns. Assuming atomic." << std::endl;
                        index_x = 2; index_y = 3; index_z = 4;
                        break;
                }
            }

            particle_attributes = num_columns;
        }
        else
        {
            // SKIP THIS SECTION — READ LINES UNTIL THE NEXT SECTION KEYWORD
            bool found_next = false;
            while(getline(fp, line))
            {
                header_lines++;
                size_t s = line.find_first_not_of(" \t\r");
                if(s == std::string::npos) continue;
                if(std::isalpha(line[s])) { found_next = true; break; }
            }
            if(!found_next)
                throw std::runtime_error("Atoms section not found in LAMMPS data file.");
        }
    }
}


////////////////////////////////////////////////////
////
////   PARSE HEADER: DETECT FILE FORMAT AND DISPATCH
////   TO THE APPROPRIATE PARSER.
////
////   LAMMPS dump files are identified by "ITEM:" on
////   the first non-blank line.  All other files are
////   assumed to be LAMMPS data files.
////
////////////////////////////////////////////////////

void parse_header(std::ifstream& fp)
{
    if (!fp.is_open())
        throw std::runtime_error("Error opening file");

    // READ FIRST LINE TO DETECT FILE FORMAT
    std::string first_line;
    getline(fp, first_line);
    fp.seekg(0);  // REWIND

    if (first_line.find("ITEM:") != std::string::npos)
    {
        file_format = 0;  // LAMMPS DUMP
        parse_lammps_dump(fp);
    }
    else
    {
        file_format = 1;  // LAMMPS DATA
        parse_lammps_data(fp);
    }

    // COMMON VALIDATION
    if(number_of_particles <= 0)
        throw std::runtime_error("No particles found in file.");

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
////   COORDINATE TYPES.  WORKS FOR BOTH LAMMPS
////   DUMP AND DATA FILE FORMATS.
////
////////////////////////////////////////////////////

void import_data()
{
    FILE *in_file = fopen(filename_data.c_str(), "r");
    if (!in_file)
        throw std::runtime_error("Error opening file: " + filename_data);

    char skip_buf[1024];
    for(int i = 0; i < header_lines; i++)
        fgets(skip_buf, sizeof(skip_buf), in_file);

    double Lx = xhi - xlo;
    double Ly = yhi - ylo;
    double Lz = zhi - zlo;

    for(int c = 0; c < number_of_particles; c++)
    {
        int id = c + 1;
        int type = 1;
        double x = 0, y = 0, z = 0;
        double junk;

        for(int d = 0; d < particle_attributes; d++)
        {
            int result;
            if     (d == index_id)   result = fscanf(in_file, "%d",  &id);
            else if(d == index_type) result = fscanf(in_file, "%d",  &type);
            else if(d == index_x)    result = fscanf(in_file, "%lg", &x);
            else if(d == index_y)    result = fscanf(in_file, "%lg", &y);
            else if(d == index_z)    result = fscanf(in_file, "%lg", &z);
            else                     result = fscanf(in_file, "%lg", &junk);

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

        // STORE PARTICLE ID AND TYPE
        if(index_id != -1) particle_ids[c] = id;
        else               particle_ids[c] = c + 1;
        particle_types[c] = type;
    }

    fclose(in_file);
}
