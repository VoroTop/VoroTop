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
#include <map>

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
////   PARSE EXTENDED XYZ FILE HEADER.
////
////   Line 1: atom count.
////   Line 2: key=value pairs including Lattice=
////   (9 numbers: lattice vectors a, b, c) and
////   Properties= (colon-separated triplets defining
////   column names, types, and widths).
////
////////////////////////////////////////////////////

static void parse_extended_xyz(std::ifstream& fp)
{
    std::string line;
    header_lines = 2;
    scaled_coordinates = 0;
    triclinic_crystal_system = 0;

    // LINE 1: ATOM COUNT
    getline(fp, line);
    number_of_particles = std::stoi(line);

    // LINE 2: KEY=VALUE PAIRS
    getline(fp, line);

    // EXTRACT Lattice="..." FROM LINE 2.
    // THE 9 NUMBERS ARE THREE LATTICE VECTORS a, b, c IN SEQUENCE.
    size_t lat_pos = line.find("Lattice=");
    if(lat_pos == std::string::npos)
        lat_pos = line.find("lattice=");

    if(lat_pos != std::string::npos)
    {
        // FIND THE QUOTED VALUE
        size_t q1 = line.find('"', lat_pos);
        if(q1 == std::string::npos)
            throw std::runtime_error("Malformed Lattice= in Extended XYZ (missing opening quote).");
        size_t q2 = line.find('"', q1 + 1);
        if(q2 == std::string::npos)
            throw std::runtime_error("Malformed Lattice= in Extended XYZ (missing closing quote).");

        std::string lattice_str = line.substr(q1 + 1, q2 - q1 - 1);
        double v[9];
        std::istringstream lss(lattice_str);
        for(int i = 0; i < 9; i++)
        {
            if(!(lss >> v[i]))
                throw std::runtime_error("Lattice= must contain 9 numbers.");
        }

        // LATTICE VECTORS: a = (v[0],v[1],v[2]), b = (v[3],v[4],v[5]), c = (v[6],v[7],v[8])
        // CHECK IF ORTHOGONAL (ALL OFF-DIAGONAL COMPONENTS ARE ZERO)
        if(v[1] != 0 || v[2] != 0 || v[3] != 0 || v[5] != 0 || v[6] != 0 || v[7] != 0)
        {
            triclinic_crystal_system = 1;
            std::cerr << "VoroTop does not currently support triclinic systems." << std::endl;
            std::cerr << "Please convert to an orthogonal box." << std::endl;
            exit(1);
        }

        // ORTHOGONAL BOX WITH ORIGIN AT (0,0,0)
        xlo = 0;  xhi = v[0];
        ylo = 0;  yhi = v[4];
        zlo = 0;  zhi = v[8];
    }
    else
    {
        throw std::runtime_error("Extended XYZ file must contain Lattice= specification.");
    }

    // EXTRACT Properties="..." FROM LINE 2.
    // FORMAT: name:type:ncols:name:type:ncols:...
    // TYPE CODES: S=STRING, R=REAL, I=INTEGER, L=LOGICAL
    size_t prop_pos = line.find("Properties=");
    if(prop_pos == std::string::npos)
        prop_pos = line.find("properties=");

    std::vector<std::string> prop_names;
    std::vector<char>        prop_types;
    std::vector<int>         prop_ncols;

    if(prop_pos != std::string::npos)
    {
        // FIND THE VALUE (MAY OR MAY NOT BE QUOTED)
        size_t val_start = prop_pos + 11;  // LENGTH OF "Properties="
        std::string prop_str;

        if(val_start < line.size() && line[val_start] == '"')
        {
            size_t q2 = line.find('"', val_start + 1);
            if(q2 == std::string::npos)
                throw std::runtime_error("Malformed Properties= (missing closing quote).");
            prop_str = line.substr(val_start + 1, q2 - val_start - 1);
        }
        else
        {
            // UNQUOTED: VALUE RUNS UNTIL NEXT WHITESPACE
            size_t end = line.find_first_of(" \t", val_start);
            prop_str = (end != std::string::npos) ? line.substr(val_start, end - val_start) : line.substr(val_start);
        }

        // PARSE COLON-SEPARATED TRIPLETS
        std::istringstream pss(prop_str);
        std::string token;
        std::vector<std::string> tokens;
        while(getline(pss, token, ':'))
            tokens.push_back(token);

        if(tokens.size() % 3 != 0)
            throw std::runtime_error("Properties= must contain triplets of name:type:ncols.");

        for(size_t i = 0; i < tokens.size(); i += 3)
        {
            prop_names.push_back(tokens[i]);
            char type_char = std::toupper(tokens[i+1][0]);
            prop_types.push_back(type_char);
            prop_ncols.push_back(std::stoi(tokens[i+2]));
        }
    }
    else
    {
        // NO Properties= FOUND: ASSUME PLAIN XYZ (SPECIES X Y Z)
        prop_names = {"species", "pos"};
        prop_types = {'S', 'R'};
        prop_ncols = {1, 3};
    }

    // MAP PROPERTIES TO COLUMN INDICES AND BUILD column_types VECTOR
    index_id      = -1;
    index_type    = -1;
    index_species = -1;
    index_x       = -1;
    index_y       = -1;
    index_z       = -1;

    int col = 0;
    for(size_t i = 0; i < prop_names.size(); i++)
    {
        std::string name = prop_names[i];

        // CONVERT NAME TO LOWERCASE FOR MATCHING
        std::string lname = name;
        for(char& ch : lname) ch = std::tolower(ch);

        if(lname == "species" || lname == "element" || lname == "symbol")
        {
            index_species = col;
        }
        else if(lname == "z" && prop_types[i] == 'I' && prop_ncols[i] == 1)
        {
            // ATOMIC NUMBER — TREAT AS TYPE
            index_type = col;
        }
        else if(lname == "pos" || lname == "positions")
        {
            index_x = col;
            index_y = col + 1;
            if(prop_ncols[i] >= 3)
                index_z = col + 2;
        }
        else if(lname == "id")
        {
            index_id = col;
        }

        // ADD COLUMN TYPES FOR EACH COLUMN IN THIS PROPERTY
        for(int j = 0; j < prop_ncols[i]; j++)
            column_types.push_back(prop_types[i]);

        col += prop_ncols[i];
    }

    particle_attributes = col;

    if(index_x == -1 || index_y == -1)
        throw std::runtime_error("Extended XYZ file must contain position data (pos:R:3).");

    if(index_z == -1) dimension = 2;
}


////////////////////////////////////////////////////
////
////   PARSE HEADER: DETECT FILE FORMAT AND DISPATCH
////   TO THE APPROPRIATE PARSER.
////
////   LAMMPS dump: identified by "ITEM:" on line 1.
////   Extended XYZ: line 1 is an integer, line 2
////   contains "Lattice=" or "Properties=".
////   Otherwise: assumed to be LAMMPS data file.
////
////////////////////////////////////////////////////

void parse_header(std::ifstream& fp)
{
    if (!fp.is_open())
        throw std::runtime_error("Error opening file");

    // READ FIRST TWO LINES TO DETECT FILE FORMAT
    std::string first_line, second_line;
    getline(fp, first_line);
    getline(fp, second_line);
    fp.seekg(0);  // REWIND

    if (first_line.find("ITEM:") != std::string::npos)
    {
        file_format = 0;  // LAMMPS DUMP
        parse_lammps_dump(fp);
    }
    else if (second_line.find("Lattice=")   != std::string::npos ||
             second_line.find("lattice=")   != std::string::npos ||
             second_line.find("Properties=") != std::string::npos ||
             second_line.find("properties=") != std::string::npos)
    {
        file_format = 2;  // EXTENDED XYZ
        parse_extended_xyz(fp);
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
////   COORDINATE TYPES.  WORKS FOR LAMMPS DUMP,
////   LAMMPS DATA, AND EXTENDED XYZ FILE FORMATS.
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

    // SPECIES-TO-TYPE MAPPING FOR EXTENDED XYZ FILES.
    // EACH UNIQUE SPECIES STRING IS ASSIGNED A SEQUENTIAL INTEGER TYPE.
    std::map<std::string, int> species_map;
    int next_type = 1;

    for(int c = 0; c < number_of_particles; c++)
    {
        int id = c + 1;
        int type = 1;
        double x = 0, y = 0, z = 0;
        double junk;
        char str_buf[256];

        for(int d = 0; d < particle_attributes; d++)
        {
            int result;
            if     (d == index_id)      result = fscanf(in_file, "%d",  &id);
            else if(d == index_type)    result = fscanf(in_file, "%d",  &type);
            else if(d == index_species)
            {
                result = fscanf(in_file, "%255s", str_buf);
                if(result == 1)
                {
                    std::string species(str_buf);
                    if(species_map.find(species) == species_map.end())
                        species_map[species] = next_type++;
                    type = species_map[species];
                }
            }
            else if(d == index_x)       result = fscanf(in_file, "%lg", &x);
            else if(d == index_y)       result = fscanf(in_file, "%lg", &y);
            else if(d == index_z)       result = fscanf(in_file, "%lg", &z);
            else if(!column_types.empty() && d < (int)column_types.size() && column_types[d] == 'S')
                                        result = fscanf(in_file, "%255s", str_buf);
            else                        result = fscanf(in_file, "%lg", &junk);

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
