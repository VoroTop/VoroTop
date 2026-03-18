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

////   File: output.cc


#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "filters.hh"
#include "variables.hh"


void output_lammps_dump(std::string filename)
{
    std::string output_file_name = filename + ".dump";
    std::ofstream output_file(output_file_name.c_str(), std::ofstream::out);
    if (!output_file) {
        std::cerr << "Error opening file: " << output_file_name << std::endl;
        return;
    }
    if (!output_file.is_open()) {
        std::cerr << "Error opening file: " << output_file_name << std::endl;
        return;
    }

    std::string full_line;
    data_file.open(filename_data, std::ifstream::in);

    
    ////////////////////////////////////////////////////
    ////
    ////    PRINT HEADER INFO TO FILE
    ////
    ////////////////////////////////////////////////////
    
    for (int c = 0; c < header_lines; ++c)
    {
        getline(data_file, full_line);
        output_file << full_line;
        if (full_line.find("ITEM: ATOMS") == 0)
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
    
    
    ////////////////////////////////////////////////////
    ////
    ////    PRINT PARTICLE DATA TO FILE
    ////
    ////////////////////////////////////////////////////
    
    for (int c = 0; c < number_of_particles; ++c)
    {
        // OUTPUT INITIAL DATA
        getline(data_file, full_line);
        
        output_file << full_line << '\t';
        
        // OUTPUT VORONOI TOPOLOGY
        output_file << vt_structure_types[c] << '\t';
        
        // OUTPUT RESOLVED TYPE AND CLUSTER INFORMATION, AS SPECIFIED
        if (r_switch) output_file << vt_structure_types_resolved[c] << '\t';
        if (c_switch)
        {
            output_file << cluster_index[c] << '\t';
            output_file << cluster_sizes[c] << '\t';
        }
        
        output_file << '\n';
    }
    
    data_file.close();
}


void output_eps(voro::container_2d& con, std::string filename)
{
    // COLORING SCHEME FOR PARTICLES
    //  (0) PARTICLES OMITTED
    //  (1) PARTICLES COLORED BLACK
    //  (2) PARTICLES COLORED BY NUMBER OF EDGES
    //  (3) PARTICLES COLORED BY FILTER CLASSIFICATION
    //  (4) PARTICLES COLORED BY VORONOI DISTANCE FROM CENTRAL PARTICLE

    // DRAWING VORONOI CELLS OR NOT
    //  (1) DRAW VORONOI CELLS
    //  (0) DO NOT DRAW VORONOI CELLS, ONLY PARTICLES
    int draw_voronoi_cells = 1;
    if (n_switch == 1) draw_voronoi_cells = 0;
    
    // THERE ARE SEVERAL OBJECTS WHOSE DIMENSIONS WE RECORD IN DIFFERENT UNITS:
    //  1. ENTIRE SYSTEM, MEASURED IN NATURAL UNITS (ANSTROMS, ETC)
    //  2. REGION TO BE DRAWN, WE CALL THIS THE INNER WINDOW; MIGHT BE ENTIRE
    //     SYSTEM; MEASURED IN NATURAL UNITS.
    //  3. EPS FIGURE, MEASURED IN POINTS, 72 POINTS = 1 INCH
    
    // SYSTEM MEASURED IN NATURAL UNITS GIVEN BY LAMMPS (ANGSTROMS, ETC)
    double system_width  = (xhi-xlo);
    double system_height = (yhi-ylo);
    
    // IF particles_in_eps WAS NOT SPECIFIED, THEN DRAW ALL PARTICLES.  
    if(particles_in_eps==0) particles_in_eps=number_of_particles;
    
    // DEFINING INNER WINDOW BOUNDARIES, MEASURED IN NATURAL UNITS
    double scaling_factor = sqrt(double(particles_in_eps) / number_of_particles);
    double system_x_min = xlo + (1. - scaling_factor) * system_width / 2.;
    double system_x_max = xlo + (1. + scaling_factor) * system_width / 2.;
    double system_y_min = ylo + (1. - scaling_factor) * system_height / 2.;
    double system_y_max = ylo + (1. + scaling_factor) * system_height / 2.;

    double padding_x = 0.5 * (1.-scaling_factor) * system_width;
    double padding_y = 0.5 * (1.-scaling_factor) * system_height;

    
    // THIS IS SCALING INFORMATION NECESSARY FOR THE EPS FIGURE.
    // THE EPS FIGURE WILL BE SCALED TO FIT INSIDE A RECTANGLE 
    // PROPORTIONAL TO THE SIZE OF THE SYSTEM AND TO THE NUMBER 
    // OF PARTICLES IN THE EPS FIGURE.
    const double area_per_particle = 1000.;
    double linear_factor_natural_to_eps = sqrt((double)number_of_particles * area_per_particle/
                                                 (system_width * system_height));
    
    // SIZE OF PARTICLE RADIUS AND EDGE WIDTH; BOTH ARE IN EPS UNITS OF POINTS
    double particle_radius = 7.2;
    double voronoi_cell_edge_width = particle_radius / 7.2 / 1.5;
    

    // IF NOT DRAWING PARTICLES, THEN MAKE THE VORONOI CELL EDGES THICKER
    if (particle_coloring_scheme == 0) voronoi_cell_edge_width *= 1.5;  
    
    // OPEN FILE FOR EPS OUTPUT
    std::string output_file_name = filename + ".eps";
    std::ofstream output_file(output_file_name.c_str(), std::ofstream::out);
    if (!output_file.is_open()) {
        std::cerr << "Error opening file: " << output_file_name << std::endl;
        return;
    }

    double eps_min_x = system_x_min * linear_factor_natural_to_eps;
    double eps_min_y = system_y_min * linear_factor_natural_to_eps;
    double eps_max_x = system_x_max * linear_factor_natural_to_eps;
    double eps_max_y = system_y_max * linear_factor_natural_to_eps;

    double figure_width  = linear_factor_natural_to_eps * system_width;
    double figure_height = linear_factor_natural_to_eps * system_height;


    // IF THE WINDOW IS ONLY MINIMALLY LARGER THAN THE SYSTEM, THEN WE 
    // DRAW ALL PARTICLES, POSSIBLY MULTIPLE TIMES DUE TO PERIODIC
    // BOUNDARIES.  THIS IS THE CASE WHEN THE WINDOW BEGINS WITHIN
    // A PARTICLE_RADIUS OF THE ORIGIN.
    if(padding_x*linear_factor_natural_to_eps < particle_radius || padding_y*linear_factor_natural_to_eps < particle_radius)
        particles_in_eps = number_of_particles;

    // OUTPUT HEADER
    output_file << "%!PS-Adobe-3.0 EPSF-3.0\n";
    output_file << "%%Creator: VoroTop\n";
    
    // FIGURE WILL ALWAYS GO FROM (0,0) TO (figure_width, figure_height)
    output_file << "%%BoundingBox: " << eps_min_x << " " << eps_min_y << " " << eps_max_x << " " << eps_max_y << "\n";
    
    // Define a function that draws a polygon given an array of coordinates
    output_file << "% Define a function that draws a polygon given an array of coordinates\n";
    output_file << "/polygon {\n";
    output_file << "  /coords exch def  % Store the coordinate array\n";
    output_file << "  newpath\n";
    output_file << "  \n";
    output_file << "  % First point is moveto\n";
    output_file << "  coords 0 get coords 1 get moveto\n";
    output_file << "  \n";
    output_file << "  % Process remaining points as lineto\n";
    output_file << "  2 2 coords length 1 sub {\n";
    output_file << "    /i exch def\n";
    output_file << "    coords i get coords i 1 add get lineto\n";
    output_file << "  } for\n";
    output_file << "  \n";
    output_file << "  closepath\n";
    output_file << "  stroke\n";
    output_file << "} def\n";
    output_file << "\n";

    // Define constant for circle radius
    output_file << "/CIRCLE_RADIUS 7.2 def % Constant radius for all circles\n";
    output_file << "\n";
    output_file << "% Define color palette (RGB values directly)\n";

    int max_colors = 0;

    if      (particle_coloring_scheme == 0) {}
    else if (particle_coloring_scheme == 1) {}
    else if (particle_coloring_scheme == 2)
    {
        max_colors = 12;
        output_file << "/colorpalette [\n";
        output_file << " 0.700 0.700 0.700   % color 0 - 0 edges\n";
        output_file << " 0.700 0.700 0.700   % color 0 - 1 edge \n";
        output_file << " 0.700 0.700 0.700   % color 0 - 2 edges\n";

        output_file << " 0.631 0.173 0.329   % color 1\n";
        output_file << " 0.812 0.373 0.396   % color 2\n";
        output_file << " 0.914 0.537 0.369   % color 3\n";
        output_file << " 0.980 0.902 0.647   % color 5\n";
        output_file << " 0.341 0.600 0.773   % color 9\n";
        output_file << " 0.753 0.878 0.718   % color 7\n";
        output_file << " 0.439 0.400 0.678   % color 10\n";
        output_file << " 0.925 0.957 0.694   % color 6\n";
        output_file << " 0.961 0.745 0.498   % color 4\n";
        output_file << " 0.541 0.788 0.710   % color 8\n";
        output_file << "] def\n";
    }
    else if (particle_coloring_scheme == 3)
    {
        // IF FILTER CONTAINS ONLY ONE STRUCTURE TYPE, THEN COLOR IT LIGHT YELLOW 
        // AND ALL ELSE IN DARK BLUE.
        if(filter_structure_types == 1)
        {
            max_colors = 1;
            output_file << "/colorpalette [\n";
            output_file << " 0.388 0.608 0.706   % color 0 - UNASSIGNED\n";
            output_file << " 0.961 0.914 0.725   % color 1 - CRYSTAL\n";
            output_file << "] def\n";
        }
        
        else
        {
        // IF FILTER HAS MULTIPLE STRUCTURE TYPES, THEN COLOR THEM AS FOLLOWS.
            max_colors = 5;
            output_file << "/colorpalette [\n";
            output_file << " 0.950 0.950 0.950   % color 0 - UNASSIGNED\n";
            output_file << " 0.961 0.914 0.725   % color 1 - CRYSTAL\n";
            output_file << " 0.388 0.608 0.706   % color 2 - GRAIN BOUNDARY\n";
            output_file << " 0.874 0.506 0.353   % color 3 - DISLOCATION\n";
            output_file << " 0.490 0.749 0.580   % color 4 - VACANCY\n";
            output_file << " 0.757 0.596 0.918   % color 5 - INTERSTITIAL\n";
            output_file << "] def\n";
        }
    }
    else if (particle_coloring_scheme == 4)
    {
        output_file << "/colorpalette [\n";
        output_file << " 0.700 0.700 0.700   % color 0 - CENTRAL PARTICLE\n";
        output_file << " 0.631 0.173 0.329   % color 1 - FIRST NEIGHBORS\n";
        output_file << " 0.812 0.373 0.396   % color 2 - SECOND NEIGHBORS\n";
        output_file << " 0.914 0.537 0.369   % color 3\n";
        output_file << " 0.961 0.745 0.498   % color 4\n";
        output_file << " 0.980 0.902 0.647   % color 5\n";
        output_file << " 0.925 0.957 0.694   % color 6\n";
        output_file << " 0.753 0.878 0.718   % color 7\n";
        output_file << " 0.541 0.788 0.710   % color 8\n";
        output_file << " 0.341 0.600 0.773   % color 9\n";
        output_file << " 0.439 0.400 0.678   % color 10\n";
        output_file << "] def\n";
    }
    else if (particle_coloring_scheme == 5 || particle_coloring_scheme == 6 || particle_coloring_scheme == 7 || particle_coloring_scheme == 8)
    {
        output_file << "/colorpalette [\n";
        output_file << " 0.900 0.900 0.900   % color 0\n";
        output_file << " 0.631 0.173 0.329   % color 1\n";
        output_file << " 0.812 0.373 0.396   % color 2\n";
        output_file << " 0.914 0.537 0.369   % color 3\n";
        output_file << " 0.961 0.745 0.498   % color 4\n";
        output_file << " 0.980 0.902 0.647   % color 5\n";
        output_file << " 0.925 0.957 0.694   % color 6\n";
        output_file << " 0.753 0.878 0.718   % color 7\n";
        output_file << " 0.541 0.788 0.710   % color 8\n";
        output_file << " 0.341 0.600 0.773   % color 9\n";
        output_file << " 0.439 0.400 0.678   % color 10\n";
        output_file << "] def\n";
    }
    else
    {
        std::cout << "Color palette not chosen for particles\n";
        exit(0);
    }
    output_file << "\n";

    if(particle_coloring_scheme == 1)
    {
        output_file << "/circle { % x y\n";
        output_file << "  /y exch def\n";
        output_file << "  /x exch def\n";
        output_file << "  0 0 0 setrgbcolor\n";
        output_file << "  newpath x y CIRCLE_RADIUS 0 360 arc fill\n";
        output_file << "} def\n";
    }
    else
    {
        output_file << "% Define functions for drawing filled and stroked circles\n";
        output_file << "/filledcircle { % x y colorindex\n";
        output_file << "  /colorindex exch def\n";
        output_file << "  /y exch def\n";
        output_file << "  /x exch def\n";
        output_file << "  \n";
        output_file << "  % Set color from color palette\n";
        output_file << "  colorindex 3 mul colorpalette exch get\n";
        output_file << "  colorindex 3 mul 1 add colorpalette exch get\n";
        output_file << "  colorindex 3 mul 2 add colorpalette exch get\n";
        output_file << "  setrgbcolor\n";
        output_file << "  \n";
        output_file << "  % Draw filled circle\n";
        output_file << "  newpath\n";
        output_file << "  x y CIRCLE_RADIUS 0 360 arc\n";
        output_file << "  fill\n";
        output_file << "} def\n";
        output_file << "\n";
        output_file << "/strokedcircle { % x y\n";
        output_file << "  /y exch def\n";
        output_file << "  /x exch def\n";
        output_file << "  \n";
        output_file << "  % Black stroke for all circles\n";
        output_file << "  0 0 0 setrgbcolor\n";
        output_file << "  \n";
        output_file << "  % Draw stroked circle\n";
        output_file << "  newpath\n";
        output_file << "  x y CIRCLE_RADIUS 0 360 arc\n";
        output_file << "  stroke\n";
        output_file << "} def\n";
        output_file << "\n";
        output_file << "% Define a function that draws a filled and stroked circle\n";
        output_file << "/circle { % x y colorindex\n";
        output_file << "  3 copy % duplicate all parameters\n";
        output_file << "  filledcircle % draw filled circle\n";
        output_file << "  pop % remove colorindex\n";
        output_file << "  strokedcircle % draw stroked circle\n";
        output_file << "} def\n";
    }

    output_file << "\n";
    
    // DRAW VORONOI CELLS
    if (draw_voronoi_cells == 1)
    {
        // SET EDGE WIDTH AND COLOR
        output_file << voronoi_cell_edge_width << " setlinewidth\n";
        output_file << "0 0 0 setrgbcolor\n";
        
        // ITERATE OVER ALL PARTICLES IN CONTAINER
        voro::voronoicell_neighbor_2d c(con);
        for (auto cli = con.begin(); cli < con.end(); ++cli)
        {
            if (con.compute_cell(c, cli))
            {
                int ijk = cli->ijk, q = cli->q;
                int pid = con.id[ijk][q];

                // COORDINATES OF CENTRAL PARTICLE FOR THIS VORONOI CELL
                double x = particle_coordinates[2 * pid];
                double y = particle_coordinates[2 * pid + 1];
                
                int to_draw = 0;

                // DETERMINE IF COPIES OF VORONOI CELL ARE NEEDED
                int min_xb = 0;
                int max_xb = 0;
                int min_yb = 0;
                int max_yb = 0;
                
                int k = 0;
                do {
                    double corner_x = x + 0.5 * c.pts[2 * k];
                    double corner_y = y + 0.5 * c.pts[2 * k + 1];
                    k = c.ed[2 * k];
                    
                    if (corner_x > system_x_min && corner_y > system_y_min && corner_x < system_x_max && corner_y < system_y_max) to_draw = 1;
                
                    if (corner_x < system_x_min) max_xb = +1;
                    if (corner_y < system_y_min) max_yb = +1;
                    if (corner_x > system_x_max) min_xb = -1;
                    if (corner_y > system_y_max) min_yb = -1;
                } while (k != 0);

                if(to_draw == 0) continue;

                for (int s = min_xb; s <= max_xb; ++s) for (int t = min_yb; t <= max_yb; ++t)
                {                   
                    output_file << "[";

                    int k = 0;
                    do {
                        double corner_x = x + 0.5 * c.pts[2 * k];
                        double corner_y = y + 0.5 * c.pts[2 * k + 1];
                        output_file << corner_x*linear_factor_natural_to_eps + s*figure_width << " " << corner_y*linear_factor_natural_to_eps + t*figure_height << " ";
                            k = c.ed[2 * k];
                    } while (k != 0);
                    
                    output_file << "] polygon\n";
                }
            }
            else
            {
                std::cerr << "Trouble computing Voronoi cell" << std::endl;
                exit(-1);
            }
        }
    }
    
    // IN THIS COLORING SCHEME, PARTICLES ARE NOT DRAWN
    if (particle_coloring_scheme == 0)
    {
        output_file << "showpage\n\n";
        output_file.close();
        return;
    }
    

    // DRAW PARTICLES; DOES NOT REQUIRE ITERATING THROUGH CONTAINER
    for (int pid = 0; pid < number_of_particles; ++pid)
    {        
        // COORDINATES OF PARTICLES IN NATURAL UNITS
        double x = particle_coordinates[2 * pid];
        double y = particle_coordinates[2 * pid + 1];
        
        double x_in_eps_units = linear_factor_natural_to_eps * x;
        double y_in_eps_units = linear_factor_natural_to_eps * y;
        
        int min_xb = 0;
        int max_xb = 0;
        int min_yb = 0;
        int max_yb = 0;

        // DO NOT DRAW PARTICLES THAT DO NOT INTERSECT THE INNER WINDOW.
        if (particles_in_eps != number_of_particles)
        {
            if (x_in_eps_units < eps_min_x - particle_radius) continue;
            if (y_in_eps_units < eps_min_y - particle_radius) continue;
            if (x_in_eps_units > eps_max_x + particle_radius) continue;
            if (y_in_eps_units > eps_max_y + particle_radius) continue;
        }
        
        // WHEN DRAWING THE ENTIRE SYSTEM WITH PERIODIC BOUNDARY CONDITIONS, THERE MAY
        // BE PARTICLES AT AN EDGE OF THE SYSTEM THAT SHOULD ALSO BE DRAWN ON THE
        // OPPOSITE SIDE, AS THE PARTICLE MAY INTERSECT BOTH EDGES.  IN SOME CASES THEY
        // SHOULD BE DRAWN FOUR TIMES.
        else
        {
            if (x_in_eps_units < eps_min_x + particle_radius) max_xb = +1;
            if (y_in_eps_units < eps_min_y + particle_radius) max_yb = +1;
            if (x_in_eps_units > eps_max_x - particle_radius) min_xb = -1;
            if (y_in_eps_units > eps_max_y - particle_radius) min_yb = -1;
        }
        
        int voronoi_cell_sides = cell_neighbor_count[pid];
        
        // IN THIS LOOP EACH PARTICLE IS DRAWN ONCE, TWICE, OR FOUR TIMES
        for (int s = min_xb; s <= max_xb; ++s) for (int t = min_yb; t <= max_yb; ++t)
        {
            double x = x_in_eps_units + s*figure_width;
            double y = y_in_eps_units + t*figure_height;

            // DRAW PARTICLES ALL IN BLACK
            if (particle_coloring_scheme == 1)
                output_file << x << " " << y << " circle\n";   
            
            // COLOR PARTICLES ACCORDING TO NUMBER OF EDGES
            else if (particle_coloring_scheme == 2)
            {
                if(voronoi_cell_sides<=max_colors) output_file << x << " " << y << " " << voronoi_cell_sides << " circle\n";   
                else                               output_file << x << " " << y << " " << 0                  << " circle\n";   
            }
            
            // COLOR PARTICLES ACCORDING TO FILTER INDEXES
            else if (particle_coloring_scheme == 3)
            {
                if(vt_structure_types[pid]<=max_colors) output_file << x << " " << y << " " << vt_structure_types[pid] << " circle\n";   
                else                                    output_file << x << " " << y << " " << 0                       << " circle\n";
            }
            
            // COLOR PARTICLES ACCORDING TO VORONOI DISTANCE FROM CENTRAL PARTICLE
            else if (particle_coloring_scheme == 4)
                output_file << x << " " << y << " " << (ring_index[pid]-1)%10 + 1 << " circle\n";
                              

            // COLOR PARTICLES ACCORDING TO DEFECT CLUSTER ID; THIS FEATURE IS DESIGNED
            // FOR PRIMARILY CRYSTALLINE SYSTEMS
            else if (particle_coloring_scheme == 5)
            {
                int color;
                if  (cluster_index[pid] > 0) color = 0;
                else cluster_index[pid] =    color = (-cluster_index[pid]) % 10 +1;;
                output_file << x << " " << y << " " << color << " circle\n";
            }
            
            // COLOR PARTICLES ACCORDING TO CRYSTAL CLUSTER ID; THIS FEATURE IS DESIGNED
            // FOR PRIMARILY DISORDERED SYSTEMS
            else if (particle_coloring_scheme == 6)
            {
                int color;
                if  (cluster_index[pid] < 0) color = 0;
                else cluster_index[pid] =    color = cluster_index[pid] % 10 +1;;
                output_file << x << " " << y << " " << color << " circle\n";
            }

            // COLOR PARTICLES ACCORDING TO DEFECT CLUSTER SIZE; THIS FEATURE IS DESIGNED
            // FOR PRIMARILY CRYSTALLINE SYSTEMS
            else if (particle_coloring_scheme == 7)
            {
                int color;
                if  (cluster_index[pid] > 0) color = 0;
                else                         color = cluster_sizes[pid] % 10 + 1;
                output_file << x << " " << y << " " << color << " circle\n";
            }

            // COLOR PARTICLES ACCORDING TO CRYSTAL CLUSTER SIZE; THIS FEATURE IS DESIGNED
            // FOR PRIMARILY DISORDERED SYSTEMS
            else if (particle_coloring_scheme == 8)
            {
                int color;
                if  (cluster_index[pid] < 0) color = 0;
                else                         color = cluster_sizes[pid] % 10 + 1;
                output_file << x << " " << y << " " << color << " circle\n";
            }

            else
            {
                std::cout << "EPS coloring scheme not chosen for particles\n";
                exit(0);
            }
        }
    }
    
    output_file << "showpage\n\n";
    output_file.close();
}


void ring_coloring()
{
    ////////////////////////////////////////////////////
    ////
    ////   LABELS ALL PARTICLES IN SYSTEM BY TOPOLOGICAL
    ////   DISTANCE FROM CENTRAL-MOST PARTICLE IN SYSTEM,
    ////   USEFUL FOR GENERATING COLOR ILLUSTRATIONS.
    ////
    ////////////////////////////////////////////////////
    
    // MAXIMUM RADIUS FOR CLUSTER BUILDING
    int max_rings = 100;
    
    // DETERMINE PARTICLE CLOSEST TO CENTER OF SYSTEM
    double closest_distance_squared_to_center = (xhi-xlo)*(yhi-ylo);
    int closest_particle_ID = 0;
    
    double cell_centerx = (xlo + xhi) / 2.;
    double cell_centery = (ylo + yhi) / 2.;
    
    for (int c = 0; c < number_of_particles; ++c)
    {
        double distance_squared = pow(particle_coordinates[2 * c] - cell_centerx, 2.) + pow(particle_coordinates[2 * c + 1] - cell_centery, 2.);
        if (distance_squared < closest_distance_squared_to_center)
        {
            closest_particle_ID = c;
            closest_distance_squared_to_center = distance_squared;
        }
    }
    
    // COMPUTE GROWING CLUSTER AROUND CENTRAL PARTICLE
    // WE WILL GROW THIS CLUSTER IN RINGS, AND USE THOSE
    // RINGS FOR THE DATA WE WANT.
    std::vector<int> cluster;
    cluster.reserve(number_of_particles);
    
    // FOR RECORDING NUMBER OF NEIGHBORS IN EACH RING
    std::vector<int> kneighbors(max_rings);
    
    // visited TRACKS WHICH PARTICLES HAVE ALREADY BEEN
    // COUNTED/INCLUDED.  INITIALLIZED WITH -1, AND THEN
    // UPDATED TO K WITH ALGORITHM
    std::vector<int> visited(number_of_particles, -1);
    
    cluster.push_back(closest_particle_ID);   // BEGIN WITH CURRENT PARTICLE
    visited[closest_particle_ID] = 0;
    kneighbors[0] = 1;
    
    // BUILD RINGS, ADDING UNVISITED NEIGHBORS OF PRIOR RING
    for (int k = 1; k < max_rings; ++k)
    {
        kneighbors[k] = 0;
        
        // WE ONLY NEED TO CONSIDER PARTICLES IN PRIOR RING
        // THIS COMPUTES THEIR INDICES IN THE CLUSTER
        int begin = 0;
        for (int c = 0; c < k - 1; ++c)
            begin += kneighbors[c];
        int end = begin + kneighbors[k - 1];
        
        // ITERATE OVER ALL PARTICLES IN PRIOR RING
        for (int c = begin; c < end; ++c)
        {
            // ID OF CURRENT PARTICLE
            int tempc = cluster[c];
            
            // ITERATE OVER ALL ITS NEIGHBORS
            for (int d = 0; d < cell_neighbor_count[tempc]; ++d)
            {
                if (visited[list_of_neighbors[tempc][d]] == -1)
                {
                    ring_index[list_of_neighbors[tempc][d]] = k;
                    visited[list_of_neighbors[tempc][d]] = k;
                    cluster.push_back(list_of_neighbors[tempc][d]);
                    kneighbors[k]++;
                }
            }
        }
    }
}
