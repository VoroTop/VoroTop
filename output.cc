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
    
    // CONSTRUCT COLOR PALETTE BASED ON CHOICE OF
    // PARTICLE COLORING SCHEME
    const int max_colors = 256;
    std::string color_strings[max_colors];
    
    // DO NOT DRAW PARTICLES; NO NEED TO SET PALETTE
    if (particle_coloring_scheme == 0) {}
    
    // ALL PARTICLES DRAWN BLACK; NO NEED TO SET PALETTE
    else if (particle_coloring_scheme == 1) {}
    
    // PALETTE FOR COLORING BY EDGE COUNT
    else if (particle_coloring_scheme == 2)
    {
        // DEFAULT GREY COLORING
        for (int j = 0; j < max_colors; ++j)
            color_strings[j] = "0.95 0.95 0.95";
        
        // VORONOI CELLS WITH FEWER THAN 4 OR MORE THAN 8 EDGES ARE COLORED GREY
        color_strings[4] = "0.490 0.749 0.580"; //#7DC094
        color_strings[5] = "0.874 0.506 0.353"; //#E0825A
        color_strings[6] = "0.961 0.914 0.725"; //#F6EABA
        color_strings[7] = "0.388 0.608 0.706"; //#639CB5
        color_strings[8] = "0.757 0.596 0.918"; //#C299EB
    }
    
    // PALETTE FOR COLORING BY FILTER ANALYSIS
    else if (particle_coloring_scheme == 3)
    {
        // DEFAULT GREY COLORING
        for (int j = 0; j < max_colors; ++j)
            color_strings[j] = "0.95 0.95 0.95";
        
        // PARTICLES WITH CLASSIFICATION INDICES GREATER THAN 4 ARE COLORED GREY
        color_strings[1] = "0.961 0.914 0.725"; //#F6EABA YELLOW    CRYSTAL
        color_strings[2] = "0.388 0.608 0.706"; //#639CB5 BLUE      GRAIN BOUNDARY
        color_strings[3] = "0.874 0.506 0.353"; //#E0825A RED       DISLOCATION
        color_strings[4] = "0.490 0.749 0.580"; //#7DC094 GREEN     VACANCY
        color_strings[5] = "0.757 0.596 0.918"; //#C299EB PURPLE    INTERSTITIAL
    }
    
    // PALETTE FOR VORONOI DISTANCES
    else if (particle_coloring_scheme == 4)
    {
        // PALETTE COLORS FOR VORONOI DISTANCES
        color_strings[0] = "0.700 0.700 0.700"; // CENTRAL PARTICLE
        color_strings[1] = "0.631 0.173 0.329"; // FIRST NEIGHBORS
        color_strings[2] = "0.812 0.373 0.396"; // SECOND NEIGHBORS
        color_strings[3] = "0.914 0.537 0.369"; // ...
        color_strings[4] = "0.961 0.745 0.498";
        color_strings[5] = "0.980 0.902 0.647";
        color_strings[6] = "0.925 0.957 0.694";
        color_strings[7] = "0.753 0.878 0.718";
        color_strings[8] = "0.541 0.788 0.710";
        color_strings[9] = "0.341 0.600 0.773";
        color_strings[10] = "0.439 0.400 0.678";
        
        for (int c = 11; c < max_colors; ++c)                // REPEATS EVERY 10
            color_strings[c] = color_strings[(c - 1) % 10 + 1];
    }
    
    else
    {
        std::cout << "Color palette not chosen for particles\n";
        exit(0);
    }
    
    
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
    if(eps_min_x < particle_radius || eps_min_y < particle_radius)
        particles_in_eps = number_of_particles;

    // OUTPUT HEADER
    output_file << "%!PS-Adobe-3.0 EPSF-3.0\n";
    output_file << "%%Creator: VoroTop\n";
    
    // FIGURE WILL ALWAYS GO FROM (0,0) TO (figure_width, figure_height)
    output_file << "%%BoundingBox: " << eps_min_x << " " << eps_min_y << " " << eps_max_x << " " << eps_max_y << "\n";
    
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

                double corner_x = x + 0.5 * c.pts[2 * k];
                double corner_y = y + 0.5 * c.pts[2 * k + 1];
                k = c.ed[2 * k];

                if (corner_x > padding_x)              to_draw = 1;
                if (corner_y > padding_y)              to_draw = 1;
                if (corner_x < system_x_max-padding_x) to_draw = 1;
                if (corner_y < system_y_max-padding_y) to_draw = 1;
                
                if (corner_x < -2.*padding_x)               max_xb = +1;
                if (corner_y < -2.*padding_y)               max_yb = +1;
                if (corner_x > system_x_max + 2.*padding_x) min_xb = -1;
                if (corner_y > system_y_max + 2.*padding_y) min_yb = -1;
                
                do {
                    double corner_x = x + 0.5 * c.pts[2 * k];
                    double corner_y = y + 0.5 * c.pts[2 * k + 1];
                    k = c.ed[2 * k];
                    
                    if (corner_x > padding_x)              to_draw = 1;
                    if (corner_y > padding_y)              to_draw = 1;
                    if (corner_x < system_x_max-padding_x) to_draw = 1;
                    if (corner_y < system_y_max-padding_y) to_draw = 1;
                    
                    if (corner_x < -2.*padding_x)               max_xb = +1;
                    if (corner_y < -2.*padding_y)               max_yb = +1;
                    if (corner_x > system_x_max + 2.*padding_x) min_xb = -1;
                    if (corner_y > system_y_max + 2.*padding_y) min_yb = -1;
                } while (k != 0);
            
                if(to_draw == 0) continue;

                // WHEN THE INNER_WINDOW IS SMALLER THAN THE SYSTEM, WE JUST
                for (int s = min_xb; s <= max_xb; ++s) for (int t = min_yb; t <= max_yb; ++t)
                {                   
                    output_file << "newpath\n";

                    int k = 0;
                    double corner_x = x + 0.5 * c.pts[2 * k];
                    double corner_y = y + 0.5 * c.pts[2 * k + 1];
                    
                    output_file << corner_x*linear_factor_natural_to_eps + s*figure_width << '\t' << corner_y*linear_factor_natural_to_eps + t*figure_height << " moveto\n";
                    k = c.ed[2 * k];
                    
                    do {
                        double corner_x = x + 0.5 * c.pts[2 * k];
                        double corner_y = y + 0.5 * c.pts[2 * k + 1];

                        output_file << corner_x*linear_factor_natural_to_eps + s*figure_width << '\t' << corner_y*linear_factor_natural_to_eps + t*figure_height << " lineto\n";
                            k = c.ed[2 * k];
                    } while (k != 0);
                    
                    output_file << "closepath\n";
                    output_file << "stroke\n";
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
            if (x_in_eps_units < 2*particle_radius)                 max_xb = +1;
            if (y_in_eps_units < 2*particle_radius)                 max_yb = +1;
            if (x_in_eps_units > figure_width  - 2*particle_radius) min_xb = -1;
            if (y_in_eps_units > figure_height - 2*particle_radius) min_yb = -1;
        }
        
        int voronoi_cell_sides = cell_neighbor_count[pid];
        
        // IN THIS LOOP EACH PARTICLE IS DRAWN ONCE, TWICE, OR FOUR TIMES
        for (int s = min_xb; s <= max_xb; ++s) for (int t = min_yb; t <= max_yb; ++t)
        {
            double x = x_in_eps_units + s*figure_width;
            double y = y_in_eps_units + t*figure_height;

            // DRAW PARTICLES ALL IN BLACK
            if (particle_coloring_scheme == 1)
                output_file << x << " " << y << " " << particle_radius << " 0 360 arc fill\n";
            
            // COLOR PARTICLES ACCORDING TO NUMBER OF EDGES
            else if (particle_coloring_scheme == 2)
            {
                output_file << color_strings[voronoi_cell_sides] << " setrgbcolor\n";
                output_file << x << " " << y << " " << particle_radius << " 0 360 arc fill\n";
                output_file << "0 0 0 setrgbcolor\n";
                output_file << x << " " << y << " " << particle_radius << " 0 360 arc stroke\n";
            }
            
            // COLOR PARTICLES ACCORDING TO FILTER INDEXES
            else if (particle_coloring_scheme == 3)
            {
                output_file << color_strings[vt_structure_types[pid]] << " setrgbcolor\n";
                output_file << x << " " << y << " " << particle_radius << " 0 360 arc fill\n";
                output_file << "0 0 0 setrgbcolor\n";
                output_file << x << " " << y << " " << particle_radius << " 0 360 arc stroke\n";
            }
            
            // COLOR PARTICLES ACCORDING TO VORONOI DISTANCE FROM CENTRAL PARTICLE
            else if (particle_coloring_scheme == 4)
            {
                output_file << color_strings[ring_index[pid]] << " setrgbcolor\n";
                output_file << x << " " << y << " " << particle_radius << " 0 360 arc fill\n";
                output_file << "0 0 0 setrgbcolor\n";
                output_file << x << " " << y << " " << particle_radius << " 0 360 arc stroke\n";
            }
            
            // COLOR PARTICLES ACCORDING TO CLUSTER ID; CURRENTLY THIS FEATURE IS DESIGNED
            // FOR PRIMARILY CRYSTALLINE SYSTEMS, AND SO INDEXES OF CLUSTERS ARE NEGATIVE.
            // CODE CAN BE ADJUSTED FOR SYSTEMS THAT ARE PRIMARILY DISORDERED, AND IN WHICH
            // CLUSTERS ARE PRIMARILY CRYSTALLINE, AND HENCE HAVE POSITIVE INDEXES.
            else if (particle_coloring_scheme == 5)
            {
                if (cluster_index[pid] > 0) cluster_index[pid] = 0;
                else cluster_index[pid] = -cluster_index[pid];
                output_file << color_strings[cluster_index[pid] % 10 + 1] << " setrgbcolor\n";
                output_file << x << " " << y << " " << particle_radius << " 0 360 arc fill\n";
                output_file << "0 0 0 setrgbcolor\n";
                output_file << x << " " << y << " " << particle_radius << " 0 360 arc stroke\n";
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
