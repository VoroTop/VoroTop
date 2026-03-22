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
#include <algorithm>

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
                std::cerr << "Please convert to a restricted triclinic box." << std::endl;
                exit(1);
            }

            // DETECT RESTRICTED TRICLINIC (xy xz yz ...)
            bool is_triclinic = (line.find("xy xz yz") != std::string::npos);

            if(is_triclinic)
            {
                triclinic_crystal_system = 1;

                // RESTRICTED TRICLINIC: EACH LINE HAS lo_bound hi_bound tilt
                double xlo_bound, xhi_bound;
                double ylo_bound, yhi_bound;
                double zlo_bound, zhi_bound;

                getline(fp, line); header_lines++;
                std::istringstream(line) >> xlo_bound >> xhi_bound >> xy;

                getline(fp, line); header_lines++;
                std::istringstream(line) >> ylo_bound >> yhi_bound >> xz;

                getline(fp, line); header_lines++;
                std::istringstream(line) >> zlo_bound >> zhi_bound >> yz;

                // CONVERT BOUNDING BOX TO ACTUAL BOX EDGES
                xlo = xlo_bound - std::min({0.0, xy, xz, xy+xz});
                xhi = xhi_bound - std::max({0.0, xy, xz, xy+xz});
                ylo = ylo_bound - std::min(0.0, yz);
                yhi = yhi_bound - std::max(0.0, yz);
                zlo = zlo_bound;
                zhi = zhi_bound;
            }
            else
            {
                // ORTHOGONAL BOX: READ xlo xhi, ylo yhi, zlo zhi
                triclinic_crystal_system = 0;
                xy = xz = yz = 0;

                getline(fp, line); header_lines++;
                std::istringstream(line) >> xlo >> xhi;

                getline(fp, line); header_lines++;
                std::istringstream(line) >> ylo >> yhi;

                getline(fp, line); header_lines++;
                std::istringstream(line) >> zlo >> zhi;
            }

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
            std::istringstream(content) >> xy >> xz >> yz;
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
        // RESTRICTED TRICLINIC REQUIRES LOWER-TRIANGULAR FORM:
        //   a = (bx, 0, 0), b = (bxy, by, 0), c = (bxz, byz, bz)
        // I.E., UPPER TRIANGLE (v[1], v[2], v[5]) MUST BE ZERO.
        if(v[1] != 0 || v[2] != 0 || v[5] != 0)
        {
            std::cerr << "Extended XYZ lattice vectors are not in restricted triclinic form." << std::endl;
            std::cerr << "The first vector must be along x, the second in the xy plane." << std::endl;
            std::cerr << "Please reorient the lattice vectors." << std::endl;
            exit(1);
        }

        if(v[3] != 0 || v[6] != 0 || v[7] != 0)
        {
            // TRICLINIC: LOWER TRIANGLE HAS NON-ZERO TILT FACTORS
            triclinic_crystal_system = 1;
            xy = v[3];   // x-COMPONENT OF SECOND LATTICE VECTOR
            xz = v[6];   // x-COMPONENT OF THIRD LATTICE VECTOR
            yz = v[7];   // y-COMPONENT OF THIRD LATTICE VECTOR
        }
        else
        {
            triclinic_crystal_system = 0;
            xy = xz = yz = 0;
        }

        // BOX WITH ORIGIN AT (0,0,0)
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
////   CIF HELPER FUNCTIONS.
////
////   strip_cif_uncertainty: Remove parenthesized
////   uncertainty from CIF numbers, e.g. "3.52(5)"
////   becomes "3.52".
////
////   tokenize_cif_line: Split a CIF data line on
////   whitespace, respecting single- and double-quoted
////   strings as single tokens.
////
////////////////////////////////////////////////////

static std::string strip_cif_uncertainty(const std::string& s)
{
    size_t p = s.find('(');
    if(p != std::string::npos)
        return s.substr(0, p);
    return s;
}

static std::vector<std::string> tokenize_cif_line(const std::string& line)
{
    std::vector<std::string> tokens;
    size_t i = 0;
    while(i < line.size())
    {
        // SKIP WHITESPACE
        while(i < line.size() && (line[i] == ' ' || line[i] == '\t'))
            i++;
        if(i >= line.size()) break;

        // QUOTED STRING
        if(line[i] == '\'' || line[i] == '"')
        {
            char quote = line[i];
            i++;
            size_t start = i;
            while(i < line.size() && line[i] != quote)
                i++;
            tokens.push_back(line.substr(start, i - start));
            if(i < line.size()) i++;  // SKIP CLOSING QUOTE
        }
        else
        {
            // UNQUOTED TOKEN
            size_t start = i;
            while(i < line.size() && line[i] != ' ' && line[i] != '\t')
                i++;
            tokens.push_back(line.substr(start, i - start));
        }
    }
    return tokens;
}


////////////////////////////////////////////////////
////
////   CIF SYMMETRY OPERATION PARSER.
////
////   Each symmetry operation is a string like
////   "-x+1/2, y, -z+1/4".  This is parsed into
////   an affine transformation:
////
////     x' = rot[0]*x + rot[1]*y + rot[2]*z + trans
////
////   The rotation coefficients are always -1, 0, or 1.
////   The translation is a rational fraction (0, 1/6,
////   1/4, 1/3, 1/2, 2/3, 3/4, 5/6).
////
////////////////////////////////////////////////////

struct SymOp {
    double rot[3][3];   // ROTATION/REFLECTION MATRIX
    double trans[3];    // TRANSLATION VECTOR
};

static void parse_symop_component(const std::string& expr,
                                  double& cx, double& cy, double& cz, double& constant)
{
    cx = cy = cz = constant = 0.0;
    int sign = 1;
    size_t i = 0;

    while(i < expr.size())
    {
        char ch = expr[i];

        // SKIP WHITESPACE
        if(ch == ' ' || ch == '\t') { i++; continue; }

        // SIGN
        if(ch == '+') { sign = 1;  i++; continue; }
        if(ch == '-') { sign = -1; i++; continue; }

        // VARIABLE: x, y, z
        if(ch == 'x' || ch == 'X') { cx += sign; sign = 1; i++; continue; }
        if(ch == 'y' || ch == 'Y') { cy += sign; sign = 1; i++; continue; }
        if(ch == 'z' || ch == 'Z') { cz += sign; sign = 1; i++; continue; }

        // NUMBER: COULD BE INTEGER, DECIMAL, OR FRACTION
        if(std::isdigit(ch) || ch == '.')
        {
            size_t start = i;
            while(i < expr.size() && (std::isdigit(expr[i]) || expr[i] == '.'))
                i++;

            double num = std::stod(expr.substr(start, i - start));

            // SKIP WHITESPACE BEFORE POSSIBLE '/'
            while(i < expr.size() && (expr[i] == ' ' || expr[i] == '\t'))
                i++;

            // CHECK FOR FRACTION
            if(i < expr.size() && expr[i] == '/')
            {
                i++;
                while(i < expr.size() && (expr[i] == ' ' || expr[i] == '\t'))
                    i++;
                size_t dstart = i;
                while(i < expr.size() && (std::isdigit(expr[i]) || expr[i] == '.'))
                    i++;
                double denom = std::stod(expr.substr(dstart, i - dstart));
                num /= denom;
            }

            // CHECK IF FOLLOWED BY A VARIABLE (E.G. "2x" — RARE BUT VALID)
            while(i < expr.size() && (expr[i] == ' ' || expr[i] == '\t'))
                i++;
            if(i < expr.size() && (expr[i] == 'x' || expr[i] == 'X'))
                { cx += sign * num; sign = 1; i++; }
            else if(i < expr.size() && (expr[i] == 'y' || expr[i] == 'Y'))
                { cy += sign * num; sign = 1; i++; }
            else if(i < expr.size() && (expr[i] == 'z' || expr[i] == 'Z'))
                { cz += sign * num; sign = 1; i++; }
            else
                { constant += sign * num; sign = 1; }

            continue;
        }

        // SKIP UNRECOGNIZED CHARACTERS
        i++;
    }
}

static SymOp parse_symop(const std::string& op_str)
{
    SymOp op = {};

    // SPLIT ON COMMAS TO GET THREE COMPONENTS
    std::vector<std::string> components;
    std::istringstream css(op_str);
    std::string comp;
    while(getline(css, comp, ','))
        components.push_back(comp);

    if(components.size() != 3)
        throw std::runtime_error("CIF: symmetry operation must have 3 components: " + op_str);

    for(int i = 0; i < 3; i++)
    {
        double cx, cy, cz, constant;
        parse_symop_component(components[i], cx, cy, cz, constant);
        op.rot[i][0] = cx;
        op.rot[i][1] = cy;
        op.rot[i][2] = cz;
        op.trans[i]  = constant;
    }

    return op;
}


////////////////////////////////////////////////////
////
////   CIF DUPLICATE DETECTION.
////
////   After applying symmetry operations, atoms on
////   special positions or cell boundaries may produce
////   duplicate positions.  Two fractional coordinates
////   are considered identical if they differ by less
////   than a tolerance, accounting for periodic wrapping
////   near 0/1.
////
////////////////////////////////////////////////////

static double frac_wrap(double f)
{
    f = fmod(f, 1.0);
    if(f < 0) f += 1.0;
    if(f >= 1.0 - 1e-6) f = 0.0;  // SNAP VALUES NEAR 1.0 TO 0.0
    return f;
}

static bool is_duplicate(double fx, double fy, double fz,
                         const std::vector<double>& coords, int n_atoms,
                         double tol)
{
    for(int i = 0; i < n_atoms; i++)
    {
        double dx = fabs(fx - coords[3*i]);
        double dy = fabs(fy - coords[3*i+1]);
        double dz = fabs(fz - coords[3*i+2]);

        // ACCOUNT FOR PERIODIC WRAPPING
        if(dx > 0.5) dx = 1.0 - dx;
        if(dy > 0.5) dy = 1.0 - dy;
        if(dz > 0.5) dz = 1.0 - dz;

        if(dx < tol && dy < tol && dz < tol)
            return true;
    }
    return false;
}


////////////////////////////////////////////////////
////
////   PARSE CIF (CRYSTALLOGRAPHIC INFORMATION FILE).
////
////   CIF files contain key-value pairs and loop_
////   sections.  The parser extracts cell parameters,
////   symmetry operations, and atom site positions.
////   Symmetry operations are applied to the asymmetric
////   unit to generate the full unit cell.
////
////   Only orthogonal cells (alpha=beta=gamma=90) are
////   currently supported.
////
////////////////////////////////////////////////////

static void parse_cif(std::ifstream& fp)
{
    scaled_coordinates = 0;
    triclinic_crystal_system = 0;
    data_already_imported = true;

    // READ ENTIRE FILE INTO LINES
    std::vector<std::string> lines;
    std::string line;
    while(getline(fp, line))
        lines.push_back(line);

    // PHASE 1: EXTRACT CELL PARAMETERS
    double cell_a = 0, cell_b = 0, cell_c = 0;
    double cell_alpha = 90, cell_beta = 90, cell_gamma = 90;

    for(size_t i = 0; i < lines.size(); i++)
    {
        std::vector<std::string> tok = tokenize_cif_line(lines[i]);
        if(tok.empty()) continue;

        if(tok[0] == "_cell_length_a" && tok.size() >= 2)
            cell_a = std::stod(strip_cif_uncertainty(tok[1]));
        else if(tok[0] == "_cell_length_b" && tok.size() >= 2)
            cell_b = std::stod(strip_cif_uncertainty(tok[1]));
        else if(tok[0] == "_cell_length_c" && tok.size() >= 2)
            cell_c = std::stod(strip_cif_uncertainty(tok[1]));
        else if(tok[0] == "_cell_angle_alpha" && tok.size() >= 2)
            cell_alpha = std::stod(strip_cif_uncertainty(tok[1]));
        else if(tok[0] == "_cell_angle_beta" && tok.size() >= 2)
            cell_beta = std::stod(strip_cif_uncertainty(tok[1]));
        else if(tok[0] == "_cell_angle_gamma" && tok.size() >= 2)
            cell_gamma = std::stod(strip_cif_uncertainty(tok[1]));
    }

    if(cell_a <= 0 || cell_b <= 0 || cell_c <= 0)
        throw std::runtime_error("CIF: missing or invalid cell dimensions.");

    // CONVERT CELL PARAMETERS (a, b, c, alpha, beta, gamma) TO LATTICE VECTORS
    // IN RESTRICTED TRICLINIC (LOWER-TRIANGULAR) FORM:
    //   a = (bx, 0, 0)
    //   b = (bxy, by, 0)
    //   c = (bxz, byz, bz)
    double deg2rad = M_PI / 180.0;
    double cos_alpha = cos(cell_alpha * deg2rad);
    double cos_beta  = cos(cell_beta  * deg2rad);
    double cos_gamma = cos(cell_gamma * deg2rad);
    double sin_gamma = sin(cell_gamma * deg2rad);

    double bx  = cell_a;
    double bxy = cell_b * cos_gamma;
    double by  = cell_b * sin_gamma;
    double bxz = cell_c * cos_beta;
    double byz = cell_c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
    double bz  = sqrt(cell_c * cell_c - bxz * bxz - byz * byz);

    if(fabs(bxy) > 0.01 || fabs(bxz) > 0.01 || fabs(byz) > 0.01)
    {
        triclinic_crystal_system = 1;
        xy = bxy;
        xz = bxz;
        yz = byz;
    }
    else
    {
        triclinic_crystal_system = 0;
        xy = xz = yz = 0;
    }

    xlo = 0;  xhi = bx;
    ylo = 0;  yhi = by;
    zlo = 0;  zhi = bz;

    // PHASE 2: EXTRACT SPACE GROUP CENTERING AND SYMMETRY OPERATIONS

    // EXTRACT SPACE GROUP NAME TO DETERMINE CENTERING TYPE
    char centering = 'P';  // DEFAULT: PRIMITIVE
    for(size_t i = 0; i < lines.size(); i++)
    {
        std::vector<std::string> tok = tokenize_cif_line(lines[i]);
        if(tok.size() >= 2 &&
           (tok[0] == "_symmetry_space_group_name_H-M" ||
            tok[0] == "_space_group_name_H-M_alt"))
        {
            // EXTRACT FIRST NON-WHITESPACE CHARACTER OF THE SPACE GROUP SYMBOL
            for(char ch : tok[1])
            {
                if(ch != ' ' && ch != '\t')
                {
                    centering = std::toupper(ch);
                    break;
                }
            }
            break;
        }
    }

    std::vector<SymOp> symops;

    for(size_t i = 0; i < lines.size(); i++)
    {
        std::string trimmed = lines[i];
        size_t start = trimmed.find_first_not_of(" \t\r");
        if(start == std::string::npos) continue;
        trimmed = trimmed.substr(start);

        if(trimmed != "loop_") continue;

        // READ LOOP HEADER TAGS
        std::vector<std::string> tags;
        size_t j = i + 1;
        while(j < lines.size())
        {
            std::string t = lines[j];
            size_t s = t.find_first_not_of(" \t\r");
            if(s == std::string::npos) { j++; continue; }
            if(t[s] == '_')
            {
                std::vector<std::string> tk = tokenize_cif_line(t);
                if(!tk.empty()) tags.push_back(tk[0]);
                j++;
            }
            else break;
        }

        // CHECK IF THIS LOOP CONTAINS SYMMETRY OPERATIONS
        int symop_col = -1;
        for(int c = 0; c < (int)tags.size(); c++)
        {
            if(tags[c] == "_symmetry_equiv_pos_as_xyz" ||
               tags[c] == "_space_group_symop_operation_xyz")
                symop_col = c;
        }
        if(symop_col == -1) continue;

        // READ SYMMETRY OPERATION DATA LINES
        while(j < lines.size())
        {
            std::string dl = lines[j];
            size_t s = dl.find_first_not_of(" \t\r");
            if(s == std::string::npos) { j++; continue; }
            if(dl[s] == '_' || dl[s] == '#') break;

            // CHECK FOR loop_ OR data_ (CASE-SENSITIVE CIF KEYWORDS)
            std::string keyword = dl.substr(s, 5);
            if(keyword == "loop_" || keyword == "data_") break;

            std::vector<std::string> vals = tokenize_cif_line(dl);
            if((int)vals.size() <= symop_col) break;

            symops.push_back(parse_symop(vals[symop_col]));
            j++;
        }
        break;  // FOUND THE SYMMETRY LOOP, NO NEED TO CONTINUE
    }

    // DEFAULT: IDENTITY OPERATION IF NO SYMMETRY FOUND
    if(symops.empty())
    {
        SymOp identity = {};
        identity.rot[0][0] = identity.rot[1][1] = identity.rot[2][2] = 1.0;
        symops.push_back(identity);
    }

    // EXPAND SYMMETRY OPERATIONS WITH LATTICE CENTERING TRANSLATIONS.
    // SOME CIF FILES LIST ONLY THE POINT GROUP OPERATIONS WITHOUT INCLUDING
    // THE CENTERING TRANSLATIONS.  WE DETECT THE CENTERING TYPE FROM THE
    // SPACE GROUP NAME AND ADD MISSING TRANSLATIONS.
    //
    // CENTERING TYPES AND THEIR TRANSLATIONS:
    //   I (BODY):   (1/2, 1/2, 1/2)
    //   F (FACE):   (1/2, 1/2, 0), (1/2, 0, 1/2), (0, 1/2, 1/2)
    //   A (A-FACE): (0, 1/2, 1/2)
    //   B (B-FACE): (1/2, 0, 1/2)
    //   C (C-FACE): (1/2, 1/2, 0)
    //   R (RHOMBOHEDRAL): (2/3, 1/3, 1/3), (1/3, 2/3, 2/3)
    //   P (PRIMITIVE): NO ADDITIONAL TRANSLATIONS

    struct CenteringTranslation { double t[3]; };
    std::vector<CenteringTranslation> centering_translations;

    if(centering == 'I')
        centering_translations = {{{0.5, 0.5, 0.5}}};
    else if(centering == 'F')
        centering_translations = {{{0.5, 0.5, 0.0}}, {{0.5, 0.0, 0.5}}, {{0.0, 0.5, 0.5}}};
    else if(centering == 'A')
        centering_translations = {{{0.0, 0.5, 0.5}}};
    else if(centering == 'B')
        centering_translations = {{{0.5, 0.0, 0.5}}};
    else if(centering == 'C')
        centering_translations = {{{0.5, 0.5, 0.0}}};
    else if(centering == 'R')
        centering_translations = {{{2.0/3, 1.0/3, 1.0/3}}, {{1.0/3, 2.0/3, 2.0/3}}};

    if(!centering_translations.empty())
    {
        // CHECK WHETHER THE CENTERING TRANSLATIONS ARE ALREADY PRESENT IN THE
        // SYMMETRY OPERATIONS.  WE TEST THE FIRST CENTERING TRANSLATION: IF AN
        // OPERATION WITH THE IDENTITY ROTATION AND THAT TRANSLATION EXISTS, THE
        // CIF ALREADY INCLUDES THE EXPANDED SET.
        bool already_expanded = false;
        double ct0 = centering_translations[0].t[0];
        double ct1 = centering_translations[0].t[1];
        double ct2 = centering_translations[0].t[2];

        for(const auto& op : symops)
        {
            // CHECK IF THIS IS IDENTITY ROTATION + CENTERING TRANSLATION
            if(fabs(op.rot[0][0] - 1.0) < 1e-6 && fabs(op.rot[1][1] - 1.0) < 1e-6 &&
               fabs(op.rot[2][2] - 1.0) < 1e-6 &&
               fabs(op.rot[0][1]) < 1e-6 && fabs(op.rot[0][2]) < 1e-6 &&
               fabs(op.rot[1][0]) < 1e-6 && fabs(op.rot[1][2]) < 1e-6 &&
               fabs(op.rot[2][0]) < 1e-6 && fabs(op.rot[2][1]) < 1e-6 &&
               fabs(fmod(op.trans[0] - ct0 + 10.0, 1.0)) < 1e-6 &&
               fabs(fmod(op.trans[1] - ct1 + 10.0, 1.0)) < 1e-6 &&
               fabs(fmod(op.trans[2] - ct2 + 10.0, 1.0)) < 1e-6)
            {
                already_expanded = true;
                break;
            }
        }

        if(!already_expanded)
        {
            // DUPLICATE EACH EXISTING OPERATION WITH EACH CENTERING TRANSLATION
            size_t orig_size = symops.size();
            for(const auto& ct : centering_translations)
            {
                for(size_t s = 0; s < orig_size; s++)
                {
                    SymOp new_op = symops[s];
                    new_op.trans[0] += ct.t[0];
                    new_op.trans[1] += ct.t[1];
                    new_op.trans[2] += ct.t[2];
                    symops.push_back(new_op);
                }
            }
        }
    }

    // PHASE 3: EXTRACT ATOM SITE POSITIONS
    struct AsymAtom {
        std::string species;
        double fx, fy, fz;
    };
    std::vector<AsymAtom> asym_atoms;

    for(size_t i = 0; i < lines.size(); i++)
    {
        std::string trimmed = lines[i];
        size_t start = trimmed.find_first_not_of(" \t\r");
        if(start == std::string::npos) continue;
        trimmed = trimmed.substr(start);

        if(trimmed != "loop_") continue;

        // READ LOOP HEADER TAGS
        std::vector<std::string> tags;
        size_t j = i + 1;
        while(j < lines.size())
        {
            std::string t = lines[j];
            size_t s = t.find_first_not_of(" \t\r");
            if(s == std::string::npos) { j++; continue; }
            if(t[s] == '_')
            {
                std::vector<std::string> tk = tokenize_cif_line(t);
                if(!tk.empty()) tags.push_back(tk[0]);
                j++;
            }
            else break;
        }

        // CHECK IF THIS LOOP CONTAINS ATOM SITES
        int col_fx = -1, col_fy = -1, col_fz = -1;
        int col_type_symbol = -1, col_label = -1;
        for(int c = 0; c < (int)tags.size(); c++)
        {
            if(tags[c] == "_atom_site_fract_x")    col_fx = c;
            if(tags[c] == "_atom_site_fract_y")    col_fy = c;
            if(tags[c] == "_atom_site_fract_z")    col_fz = c;
            if(tags[c] == "_atom_site_type_symbol") col_type_symbol = c;
            if(tags[c] == "_atom_site_label")       col_label = c;
        }
        if(col_fx == -1 || col_fy == -1 || col_fz == -1) continue;

        // DETERMINE WHICH COLUMN TO USE FOR SPECIES
        int col_species = (col_type_symbol != -1) ? col_type_symbol : col_label;

        // READ ATOM SITE DATA LINES
        while(j < lines.size())
        {
            std::string dl = lines[j];
            size_t s = dl.find_first_not_of(" \t\r");
            if(s == std::string::npos) { j++; continue; }
            if(dl[s] == '_' || dl[s] == '#') break;

            std::string keyword = dl.substr(s, 5);
            if(keyword == "loop_" || keyword == "data_") break;

            std::vector<std::string> vals = tokenize_cif_line(dl);
            int max_col = std::max({col_fx, col_fy, col_fz});
            if((int)vals.size() <= max_col) break;

            AsymAtom atom;
            atom.fx = std::stod(strip_cif_uncertainty(vals[col_fx]));
            atom.fy = std::stod(strip_cif_uncertainty(vals[col_fy]));
            atom.fz = std::stod(strip_cif_uncertainty(vals[col_fz]));

            // EXTRACT SPECIES NAME
            if(col_species != -1 && col_species < (int)vals.size())
            {
                atom.species = vals[col_species];
                // IF USING LABEL (e.g. "Fe1"), STRIP TRAILING DIGITS
                if(col_species == col_label)
                {
                    size_t end = atom.species.find_first_of("0123456789");
                    if(end != std::string::npos)
                        atom.species = atom.species.substr(0, end);
                }
            }
            else
            {
                atom.species = "X";
            }

            asym_atoms.push_back(atom);
            j++;
        }
        break;  // FOUND THE ATOM SITE LOOP
    }

    if(asym_atoms.empty())
        throw std::runtime_error("CIF: no atom sites found.");

    // PHASE 4: APPLY SYMMETRY OPERATIONS AND REMOVE DUPLICATES
    std::vector<double> frac_coords;  // FRACTIONAL COORDS OF UNIQUE ATOMS (3 PER ATOM)
    std::vector<std::string> atom_species;
    int n_unique = 0;

    for(const auto& atom : asym_atoms)
    {
        for(const auto& op : symops)
        {
            double fx = op.rot[0][0]*atom.fx + op.rot[0][1]*atom.fy + op.rot[0][2]*atom.fz + op.trans[0];
            double fy = op.rot[1][0]*atom.fx + op.rot[1][1]*atom.fy + op.rot[1][2]*atom.fz + op.trans[1];
            double fz = op.rot[2][0]*atom.fx + op.rot[2][1]*atom.fy + op.rot[2][2]*atom.fz + op.trans[2];

            fx = frac_wrap(fx);
            fy = frac_wrap(fy);
            fz = frac_wrap(fz);

            if(!is_duplicate(fx, fy, fz, frac_coords, n_unique, 0.01))
            {
                frac_coords.push_back(fx);
                frac_coords.push_back(fy);
                frac_coords.push_back(fz);
                atom_species.push_back(atom.species);
                n_unique++;
            }
        }
    }

    // ASSIGN SPECIES TO INTEGER TYPES
    std::map<std::string, int> species_map;
    int next_type = 1;

    number_of_particles = n_unique;

    // CONVERT FRACTIONAL TO CARTESIAN COORDINATES.
    // FOR LOWER-TRIANGULAR LATTICE MATRIX:
    //   x = fx*bx + fy*bxy + fz*bxz
    //   y =         fy*by  + fz*byz
    //   z =                  fz*bz
    double Lx = xhi - xlo;  // bx
    double Ly = yhi - ylo;  // by
    double Lz = zhi - zlo;  // bz

    cif_coordinates.resize(3 * number_of_particles);
    cif_types.resize(number_of_particles);

    for(int c = 0; c < number_of_particles; c++)
    {
        double fx = frac_coords[3*c];
        double fy = frac_coords[3*c + 1];
        double fz = frac_coords[3*c + 2];

        cif_coordinates[3*c]     = fx * Lx + fy * xy + fz * xz;
        cif_coordinates[3*c + 1] = fy * Ly + fz * yz;
        cif_coordinates[3*c + 2] = fz * Lz;

        if(species_map.find(atom_species[c]) == species_map.end())
            species_map[atom_species[c]] = next_type++;
        cif_types[c] = species_map[atom_species[c]];
    }
}


////////////////////////////////////////////////////
////
////   PARSE VASP POSCAR/CONTCAR FILE.
////
////   Line 1: comment.
////   Line 2: universal scale factor.
////   Lines 3-5: lattice vectors (3x3).
////   Line 6: species names (VASP 5+) or atom counts.
////   Line 7: atom counts (if line 6 was species names).
////   Optional: "Selective dynamics".
////   Next line: "Direct" or "Cartesian".
////
////////////////////////////////////////////////////

static void parse_poscar(std::ifstream& fp)
{
    std::string line;
    scaled_coordinates = 0;
    triclinic_crystal_system = 0;
    header_lines = 0;

    // LINE 1: COMMENT (SKIP)
    getline(fp, line);
    header_lines++;

    // LINE 2: SCALE FACTOR
    getline(fp, line);
    header_lines++;
    double scale_factor = std::stod(line);

    // LINES 3-5: LATTICE VECTORS
    double v[9];
    for(int i = 0; i < 3; i++)
    {
        getline(fp, line);
        header_lines++;
        std::istringstream iss(line);
        if(!(iss >> v[3*i] >> v[3*i+1] >> v[3*i+2]))
            throw std::runtime_error("POSCAR: could not read lattice vector on line " + std::to_string(header_lines) + ".");
    }

    // HANDLE NEGATIVE SCALE FACTOR (VALUE IS DESIRED CELL VOLUME)
    if(scale_factor < 0)
    {
        double vol = v[0]*(v[4]*v[8] - v[5]*v[7])
                   - v[1]*(v[3]*v[8] - v[5]*v[6])
                   + v[2]*(v[3]*v[7] - v[4]*v[6]);
        if(vol < 0) vol = -vol;
        scale_factor = pow(fabs(scale_factor) / vol, 1.0/3.0);
    }

    // SCALE LATTICE VECTORS
    for(int i = 0; i < 9; i++)
        v[i] *= scale_factor;

    // RESTRICTED TRICLINIC REQUIRES LOWER-TRIANGULAR FORM:
    //   a = (bx, 0, 0), b = (bxy, by, 0), c = (bxz, byz, bz)
    if(v[1] != 0 || v[2] != 0 || v[5] != 0)
    {
        std::cerr << "POSCAR lattice vectors are not in restricted triclinic form." << std::endl;
        std::cerr << "The first vector must be along x, the second in the xy plane." << std::endl;
        std::cerr << "Please reorient the lattice vectors." << std::endl;
        exit(1);
    }

    if(v[3] != 0 || v[6] != 0 || v[7] != 0)
    {
        triclinic_crystal_system = 1;
        xy = v[3];
        xz = v[6];
        yz = v[7];
    }
    else
    {
        triclinic_crystal_system = 0;
        xy = xz = yz = 0;
    }

    xlo = 0;  xhi = v[0];
    ylo = 0;  yhi = v[4];
    zlo = 0;  zhi = v[8];

    // LINE 6: SPECIES NAMES (VASP 5+) OR ATOM COUNTS.
    // IF THE FIRST TOKEN IS NOT A NUMBER, THIS LINE CONTAINS SPECIES NAMES
    // AND THE NEXT LINE CONTAINS ATOM COUNTS.
    getline(fp, line);
    header_lines++;

    std::vector<int> atom_counts;
    std::istringstream iss6(line);
    std::string token;
    iss6 >> token;

    bool is_number = !token.empty() && (std::isdigit(token[0]) || token[0] == '-' || token[0] == '+');
    if(is_number)
    {
        // LINE 6 CONTAINS ATOM COUNTS DIRECTLY (PRE-VASP5 FORMAT)
        atom_counts.push_back(std::stoi(token));
        int count;
        while(iss6 >> count)
            atom_counts.push_back(count);
    }
    else
    {
        // LINE 6 CONTAINS SPECIES NAMES; LINE 7 HAS ATOM COUNTS
        getline(fp, line);
        header_lines++;
        std::istringstream iss7(line);
        int count;
        while(iss7 >> count)
            atom_counts.push_back(count);
    }

    if(atom_counts.empty())
        throw std::runtime_error("POSCAR: no atom counts found.");

    // COMPUTE TOTAL PARTICLES AND BUILD TYPE ARRAY FROM HEADER
    number_of_particles = 0;
    for(int c : atom_counts)
        number_of_particles += c;

    header_assigned_types.resize(number_of_particles);
    int idx = 0;
    for(int t = 0; t < (int)atom_counts.size(); t++)
        for(int j = 0; j < atom_counts[t]; j++)
            header_assigned_types[idx++] = t + 1;

    // CHECK FOR OPTIONAL "Selective dynamics" LINE.
    // THIS LINE STARTS WITH 'S' OR 's'; COORDINATE TYPE LINES START WITH
    // 'D'/'d' (DIRECT) OR 'C'/'c'/'K'/'k' (CARTESIAN).
    getline(fp, line);
    header_lines++;

    bool selective_dynamics = false;
    size_t start = line.find_first_not_of(" \t\r");
    if(start != std::string::npos && std::tolower(line[start]) == 's')
    {
        selective_dynamics = true;
        getline(fp, line);
        header_lines++;
        start = line.find_first_not_of(" \t\r");
    }

    // COORDINATE TYPE: "Direct" (FRACTIONAL) OR "Cartesian"/"Kartesian"
    if(start == std::string::npos)
        throw std::runtime_error("POSCAR: missing coordinate type line.");

    char first = std::tolower(line[start]);
    if(first == 'd')
    {
        scaled_coordinates = 1;
        coordinate_scale = 1.0;
    }
    else if(first == 'c' || first == 'k')
    {
        scaled_coordinates = 0;
        coordinate_scale = scale_factor;
    }
    else
    {
        throw std::runtime_error("POSCAR: expected 'Direct' or 'Cartesian', got: " + line);
    }

    // SET COLUMN INDICES (POSCAR ATOM LINES: x y z [T T T])
    index_id      = -1;
    index_type    = -1;
    index_species = -1;
    index_x       = 0;
    index_y       = 1;
    index_z       = 2;

    if(selective_dynamics)
    {
        particle_attributes = 6;
        column_types = {'R', 'R', 'R', 'S', 'S', 'S'};
    }
    else
    {
        particle_attributes = 3;
    }
}


////////////////////////////////////////////////////
////
////   PARSE ATOMEYE EXTENDED CFG FILE.
////
////   The AtomEye Extended CFG format begins with
////   header lines specifying the number of particles,
////   lattice vectors (H0), a scale factor (A), and
////   auxiliary property definitions.  Coordinates are
////   always fractional.
////
////   Because species markers (mass + symbol lines)
////   are interspersed with coordinate data, this
////   parser reads all data during header parsing and
////   stores results via data_already_imported.
////
////////////////////////////////////////////////////

static void parse_atomeye_cfg(std::ifstream& fp)
{
    std::string line;
    scaled_coordinates = 1;
    triclinic_crystal_system = 0;
    data_already_imported = true;

    // LATTICE VECTORS (3x3 MATRIX, ROW-MAJOR)
    double H[3][3] = {};
    double scale = 1.0;
    int entry_count = 0;
    bool no_velocity = false;

    // PARSE HEADER LINES UNTIL entry_count IS FOUND
    while(getline(fp, line))
    {
        // SKIP BLANK LINES AND COMMENTS
        size_t start = line.find_first_not_of(" \t\r");
        if(start == std::string::npos) continue;
        if(line[start] == '#') continue;

        if(line.find("Number of particles") != std::string::npos)
        {
            size_t eq = line.find('=');
            if(eq != std::string::npos)
                number_of_particles = std::stoi(line.substr(eq + 1));
        }
        else if(line.find("A =") != std::string::npos)
        {
            size_t eq = line.find('=');
            if(eq != std::string::npos)
                scale = std::stod(line.substr(eq + 1));
        }
        else if(line.find("H0(") != std::string::npos)
        {
            // FORMAT: H0(i,j) = value [A]
            // INDICES ARE 1-BASED
            size_t p1 = line.find('(');
            size_t p2 = line.find(')');
            if(p1 != std::string::npos && p2 != std::string::npos)
            {
                std::string indices = line.substr(p1 + 1, p2 - p1 - 1);
                int i, j;
                char comma;
                std::istringstream iss(indices);
                iss >> i >> comma >> j;

                size_t eq = line.find('=');
                if(eq != std::string::npos)
                {
                    double value;
                    std::istringstream vss(line.substr(eq + 1));
                    vss >> value;
                    H[i - 1][j - 1] = value;
                }
            }
        }
        else if(line.find(".NO_VELOCITY.") != std::string::npos)
        {
            no_velocity = true;
        }
        else if(line.find("entry_count") != std::string::npos)
        {
            size_t eq = line.find('=');
            if(eq != std::string::npos)
                entry_count = std::stoi(line.substr(eq + 1));
            break;  // STOP READING HEADER AFTER entry_count
        }
    }

    // APPLY SCALE FACTOR TO LATTICE VECTORS
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            H[i][j] *= scale;

    // THE ATOMEYE CONVENTION STORES LATTICE VECTORS AS COLUMNS OF H0
    // (JU LI, MODELLING SIMUL. MATER. SCI. ENG. 11, 2003), BUT SOME
    // TOOLS (E.G. ASE) WRITE THEM AS ROWS.  FOR RESTRICTED TRICLINIC,
    // THE COLUMN CONVENTION GIVES AN UPPER-TRIANGULAR MATRIX AND THE
    // ROW CONVENTION GIVES A LOWER-TRIANGULAR MATRIX.  BOTH ARE
    // ACCEPTED SINCE THEY PRODUCE IDENTICAL CARTESIAN COORDINATES.

    bool upper_zero = (H[0][1] == 0 && H[0][2] == 0 && H[1][2] == 0);
    bool lower_zero = (H[1][0] == 0 && H[2][0] == 0 && H[2][1] == 0);

    if(!upper_zero && !lower_zero)
    {
        std::cerr << "AtomEye CFG: H0 matrix is not in restricted triclinic form." << std::endl;
        std::cerr << "H0 must be upper-triangular (column convention) or" << std::endl;
        std::cerr << "lower-triangular (row convention)." << std::endl;
        exit(1);
    }

    if(upper_zero && lower_zero)
    {
        // DIAGONAL MATRIX: ORTHOGONAL CELL
        triclinic_crystal_system = 0;
        xy = xz = yz = 0;
    }
    else if(upper_zero)
    {
        // LOWER-TRIANGULAR: ROW CONVENTION (ROWS = LATTICE VECTORS)
        triclinic_crystal_system = 1;
        xy = H[1][0];
        xz = H[2][0];
        yz = H[2][1];
    }
    else
    {
        // UPPER-TRIANGULAR: COLUMN CONVENTION (COLUMNS = LATTICE VECTORS)
        triclinic_crystal_system = 1;
        xy = H[0][1];
        xz = H[0][2];
        yz = H[1][2];
    }

    xlo = 0;  xhi = H[0][0];
    ylo = 0;  yhi = H[1][1];
    zlo = 0;  zhi = H[2][2];

    // NUMBER OF AUXILIARY PROPERTIES PER ATOM (AFTER x y z AND OPTIONAL vx vy vz)
    int aux_count;
    if(no_velocity) aux_count = entry_count - 3;
    else            aux_count = entry_count - 6;

    // SKIP AUXILIARY PROPERTY NAMES
    for(int i = 0; i < aux_count; i++)
        getline(fp, line);

    // READ PARTICLE DATA.
    // EACH SPECIES BLOCK BEGINS WITH A LINE: <mass> <symbol>
    // FOLLOWED BY COORDINATE LINES FOR ATOMS OF THAT SPECIES.
    // COORDINATE LINES: x y z [vx vy vz] [aux1 aux2 ...]
    cif_coordinates.resize(3 * number_of_particles);
    cif_types.resize(number_of_particles);

    std::map<std::string, int> species_map;
    int next_type = 1;
    int current_type = 1;
    int atoms_read = 0;

    // NUMBER OF EXTRA VALUES TO SKIP PER COORDINATE LINE (VELOCITIES + AUXILIARIES)
    int skip_count;
    if(no_velocity) skip_count = aux_count;
    else            skip_count = 3 + aux_count;

    double Lx = xhi - xlo;
    double Ly = yhi - ylo;
    double Lz = zhi - zlo;

    while(atoms_read < number_of_particles && getline(fp, line))
    {
        // SKIP BLANK LINES
        size_t start = line.find_first_not_of(" \t\r");
        if(start == std::string::npos) continue;

        // DETERMINE IF THIS IS A SPECIES LINE OR A COORDINATE LINE.
        // A SPECIES LINE HAS: <mass (number)> <symbol (alphabetic)>
        // A COORDINATE LINE HAS: <x (number)> <y (number)> ...
        std::istringstream iss(line);
        double first_val;
        iss >> first_val;

        std::string second_token;
        iss >> second_token;

        if(!second_token.empty() && std::isalpha(second_token[0]))
        {
            // THIS IS A SPECIES LINE: mass symbol
            if(species_map.find(second_token) == species_map.end())
                species_map[second_token] = next_type++;
            current_type = species_map[second_token];
        }
        else
        {
            // THIS IS A COORDINATE LINE
            double fx = first_val;
            double fy = second_token.empty() ? 0 : std::stod(second_token);
            double fz;
            iss >> fz;

            // SKIP REMAINING VALUES ON LINE (VELOCITIES AND AUXILIARIES)
            double junk;
            for(int i = 0; i < skip_count; i++)
                iss >> junk;

            // WRAP FRACTIONAL COORDINATES INTO [0, 1)
            fx = fmod(fx, 1.0); if(fx < 0) fx += 1.0;
            fy = fmod(fy, 1.0); if(fy < 0) fy += 1.0;
            fz = fmod(fz, 1.0); if(fz < 0) fz += 1.0;

            // CONVERT FRACTIONAL TO CARTESIAN
            double x, y, z;
            if(triclinic_crystal_system)
            {
                x = fx * Lx + fy * xy + fz * xz;
                y =           fy * Ly + fz * yz;
                z =                     fz * Lz;
            }
            else
            {
                x = xlo + fx * Lx;
                y = ylo + fy * Ly;
                z = zlo + fz * Lz;
            }

            cif_coordinates[3 * atoms_read]     = x;
            cif_coordinates[3 * atoms_read + 1] = y;
            cif_coordinates[3 * atoms_read + 2] = z;
            cif_types[atoms_read] = current_type;
            atoms_read++;
        }
    }

    if(atoms_read != number_of_particles)
        throw std::runtime_error("AtomEye CFG: expected " + std::to_string(number_of_particles)
                                 + " atoms but read " + std::to_string(atoms_read) + ".");
}


////////////////////////////////////////////////////
////
////   PARSE HEADER: DETECT FILE FORMAT AND DISPATCH
////   TO THE APPROPRIATE PARSER.
////
////   LAMMPS dump: identified by "ITEM:" on line 1.
////   Extended XYZ: line 2 contains "Lattice=" or
////   "Properties=".
////   AtomEye CFG: line 1 contains "Number of
////   particles".
////   CIF: a line starts with "data_" or contains
////   "_cell_length_a".
////   POSCAR: line 2 is a single number and lines 3-5
////   each contain exactly 3 numbers.
////   Otherwise: assumed to be LAMMPS data file.
////
////////////////////////////////////////////////////

void parse_header(std::ifstream& fp)
{
    if (!fp.is_open())
        throw std::runtime_error("Error opening file");

    // READ FIRST FIVE LINES TO DETECT FILE FORMAT.
    // CLEAR EOF BIT BEFORE REWINDING (FILE MAY HAVE FEWER THAN 5 LINES).
    std::string detect_lines[5];
    for(int i = 0; i < 5; i++)
        getline(fp, detect_lines[i]);
    fp.clear();
    fp.seekg(0);  // REWIND

    if (detect_lines[0].find("ITEM:") != std::string::npos)
    {
        file_format = 0;  // LAMMPS DUMP
        parse_lammps_dump(fp);
    }
    else if (detect_lines[1].find("Lattice=")   != std::string::npos ||
             detect_lines[1].find("lattice=")   != std::string::npos ||
             detect_lines[1].find("Properties=") != std::string::npos ||
             detect_lines[1].find("properties=") != std::string::npos)
    {
        file_format = 2;  // EXTENDED XYZ
        parse_extended_xyz(fp);
    }
    else if (detect_lines[0].find("Number of particles") != std::string::npos)
    {
        file_format = 5;  // ATOMEYE CFG
        parse_atomeye_cfg(fp);
    }
    else if ([&]() -> bool {
        // CIF: ANY OF THE FIRST 5 LINES STARTS WITH "data_" OR CONTAINS "_cell_length_a"
        for(int i = 0; i < 5; i++)
        {
            size_t s = detect_lines[i].find_first_not_of(" \t\r");
            if(s != std::string::npos && detect_lines[i].substr(s, 5) == "data_")
                return true;
            if(detect_lines[i].find("_cell_length_a") != std::string::npos)
                return true;
        }
        return false;
    }())
    {
        file_format = 4;  // CIF
        parse_cif(fp);
    }
    else if ([&]() -> bool {
        // POSCAR: LINE 2 IS A SINGLE NUMBER, LINES 3-5 EACH HAVE EXACTLY 3 NUMBERS
        double val;
        std::istringstream iss2(detect_lines[1]);
        std::string extra;
        if(!(iss2 >> val) || (iss2 >> extra)) return false;

        for(int i = 2; i <= 4; i++)
        {
            double a, b, c;
            std::istringstream iss(detect_lines[i]);
            if(!(iss >> a >> b >> c) || (iss >> extra)) return false;
        }
        return true;
    }())
    {
        file_format = 3;  // POSCAR
        parse_poscar(fp);
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
////   LAMMPS DATA, EXTENDED XYZ, POSCAR, CIF, AND
////   ATOMEYE CFG.
////
////////////////////////////////////////////////////

void import_data()
{
    // CIF DATA WAS ALREADY FULLY PARSED DURING parse_header()
    if(data_already_imported)
    {
        for(int c = 0; c < number_of_particles; c++)
        {
            particle_coordinates[3*c]     = cif_coordinates[3*c];
            particle_coordinates[3*c + 1] = cif_coordinates[3*c + 1];
            particle_coordinates[3*c + 2] = cif_coordinates[3*c + 2];
            particle_ids[c]   = c + 1;
            particle_types[c] = cif_types[c];
        }
        cif_coordinates.clear();
        cif_types.clear();
        return;
    }

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

        // APPLY COORDINATE SCALE (FOR POSCAR CARTESIAN MODE)
        x *= coordinate_scale;
        y *= coordinate_scale;
        if(dimension == 3) z *= coordinate_scale;

        // USE HEADER-ASSIGNED TYPES (FOR POSCAR)
        if(index_type == -1 && !header_assigned_types.empty())
            type = header_assigned_types[c];

        if(triclinic_crystal_system)
        {
            // TRICLINIC: CONVERT AND WRAP USING THE FULL LATTICE MATRIX.
            // THE LOWER-TRIANGULAR MATRIX IS:
            //   | Lx  0   0  |
            //   | xy  Ly  0  |
            //   | xz  yz  Lz |

            // STEP 1: CONVERT SCALED/FRACTIONAL TO CARTESIAN, OR SHIFT TO ORIGIN
            if(scaled_coordinates)
            {
                double fx = x, fy = y, fz = z;
                x = fx * Lx + fy * xy + fz * xz;
                y =           fy * Ly + fz * yz;
                z =                     fz * Lz;
            }
            else
            {
                x -= xlo;
                y -= ylo;
                z -= zlo;
            }

            // STEP 2: CARTESIAN → FRACTIONAL (BACK-SUBSTITUTION)
            double fz = z / Lz;
            double fy = (y - fz * yz) / Ly;
            double fx = (x - fy * xy - fz * xz) / Lx;

            // STEP 3: WRAP FRACTIONAL COORDINATES INTO [0, 1)
            fx = fmod(fx, 1.0); if(fx < 0) fx += 1.0;
            fy = fmod(fy, 1.0); if(fy < 0) fy += 1.0;
            fz = fmod(fz, 1.0); if(fz < 0) fz += 1.0;

            // STEP 4: FRACTIONAL → CARTESIAN (ORIGIN AT 0)
            x = fx * Lx + fy * xy + fz * xz;
            y =           fy * Ly + fz * yz;
            z =                     fz * Lz;
        }
        else
        {
            // ORTHOGONAL: CONVERT AND WRAP PER AXIS INDEPENDENTLY
            if(scaled_coordinates)
            {
                x = xlo + x * Lx;
                y = ylo + y * Ly;
                if(dimension == 3) z = zlo + z * Lz;
            }

            x = xlo + fmod(x - xlo, Lx);
            if(x < xlo) x += Lx;

            y = ylo + fmod(y - ylo, Ly);
            if(y < ylo) y += Ly;

            if(dimension == 3)
            {
                z = zlo + fmod(z - zlo, Lz);
                if(z < zlo) z += Lz;
            }
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
