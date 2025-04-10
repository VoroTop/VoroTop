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


void output_system(std::string filename)
{
    std::string output_file_name = filename + ".dump";
    std::ofstream output_file(output_file_name.c_str(), std::ofstream::out);
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
    std::string color_strings[50];
    
    // DO NOT DRAW PARTICLES; NO NEED TO SET PALETTE
    if (particle_coloring_scheme == 0) {}
    
    // ALL PARTICLES DRAWN BLACK; NO NEED TO SET PALETTE
    else if (particle_coloring_scheme == 1) {}
    
    // PALETTE FOR COLORING BY EDGE COUNT
    else if (particle_coloring_scheme == 2)
    {
        // DEFAULT GREY COLORING
        for (int j = 0; j < 50; ++j)
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
        for (int j = 0; j < 50; ++j)
            color_strings[j] = "0.95 0.95 0.95";
        color_strings[0] = "0.961 0.914 0.725"; //#F6EABA LIGHT YELLOW  CRYSTAL
        
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
        
        for (int c = 11; c < 50; ++c)                // REPEATS EVERY 10
            color_strings[c] = color_strings[(c - 1) % 10 + 1];
    }
    
    else
    {
        std::cout << "Color palette not chosen for particles\n";
        exit(0);
    }
    
    
    // THERE ARE SEVERAL OBJECTS WHOSE DIMENSIONS WE RECORD IN DIFFERENT UNITS:
    //  1. ENTIRE SYSTEM, MEASURED IN NATURAL UNITS
    //  2. REGION TO BE DRAWN, WE CALL THIS THE INNER WINDOW; MIGHT BE ENTIRE
    //     SYSTEM; MEASURED IN NATURAL UNITS.
    //  3. EPS FIGURE, MEASURED IN POINTS, 72 POINTS = 1 INCH
    
    // SYSTEM MEASURED IN NATURAL UNITS GIVEN BY LAMMPS (ANGSTROMS, ETC)
    double system_width = hi_bound[0] - origin[0];
    double system_height = hi_bound[1] - origin[1];
    
    // BY DEFAULT WE WILL DRAW THE ENTIRE SYSTEM; BUT IF THE USER CHOOSES TO DRAW
    // A SMALLER PART OF THE SYSTEM, THEN WE WILL ADJUST THESE WINDOWS ACCORDINGLY.
    // THIS FEATURE CAN BE USEFUL WHEN CONSIDERING VERY LARGE SYSTEMS.
    double inner_window_width = system_width;
    double inner_window_height = system_height;
    
    // IF particles_in_eps WAS NOT SPECIFIED, THEN DRAW ALL PARTICLES.  
    if(particles_in_eps==0) particles_in_eps=number_of_particles;
    
    // IF DRAWING SUBSYSTEM, THEN THE WINDOW WILL BE A SQUARE WHOSE AREA IS PROPORTIONAL
    // TO THE NUMBER OF PARTICLES WE WISH TO DRAW.
    if (particles_in_eps != number_of_particles)
    {
        double inner_window_area = double(particles_in_eps) / number_of_particles * system_width * system_height;
        
        inner_window_width = sqrt(inner_window_area);
        inner_window_height = inner_window_width;
    }
    
    // ALL MEASUREMENTS HERE ARE IN NATURAL UNITS
    double window_center_x = (origin[0] + hi_bound[0]) / 2.;
    double window_center_y = (origin[1] + hi_bound[1]) / 2.;
    
    double inner_window_origin_x = window_center_x - inner_window_width / 2.;
    double inner_window_origin_y = window_center_y - inner_window_height / 2.;
    
    // FIGURE DIMENSIONS ARE IN EPS UNITS OF POINTS
    double figure_width = 36. * sqrt(particles_in_eps * inner_window_width / inner_window_height);
    double figure_height = figure_width * inner_window_height / inner_window_width;
    
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
    
    // OUTPUT HEADER
    output_file << "%!PS-Adobe-3.0 EPSF-3.0\n";
    output_file << "%%Creator: VoroTop\n";
    
    // FIGURE WILL ALWAYS GO FROM (0,0) TO (figure_width, figure_height)
    output_file << "%%BoundingBox: " << 0. << " " << 0. << " " << figure_width << " " << figure_height << "\n";
    
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
                
                // IF DRAWING ONLY A LIMITED PART OF SYSTEM, THEN IGNORE VORONOI CELLS 
                // OF PARTICLES THAT ARE NOT NEAR THE WINDOW; THIS WILL SAVE DISK SPACE
                // FOR VERY LARGE SYSTEMS.
                if (particles_in_eps != number_of_particles)
                {
                    if (abs(x - window_center_x) > inner_window_width)  continue;
                    if (abs(y - window_center_y) > inner_window_height) continue;
                }
                
                int to_draw = 0;
                int min_xb = 0;
                int max_xb = 0;
                int min_yb = 0;
                int max_yb = 0;
                
                if (particles_in_eps == number_of_particles)
                {
                    // IN THIS CASE WE KNOW THAT WE NEED TO DRAW THE VORONOI.  THE ONLY
                    // QUESTION IS HOW MANY TIMES, 1, 2, OR 4, DEPENDING ON WHERE THESE
                    // CORNERS ARE LOCATED.
                    to_draw = 1;
                    
                    int k = 0;
                    double corner_x = figure_width * ((x + 0.5 * c.pts[2 * k]) - inner_window_origin_x) / inner_window_width;
                    double corner_y = figure_height * ((y + 0.5 * c.pts[2 * k + 1]) - inner_window_origin_y) / inner_window_height;
                    k = c.ed[2 * k];
                    
                    if (corner_x < 0) min_xb = -1;
                    if (corner_y < 0) min_yb = -1;
                    if (corner_x > figure_width)  max_xb = +1;
                    if (corner_y > figure_height) max_yb = +1;

                    do {
                        corner_x = figure_width * ((x + 0.5 * c.pts[2 * k]) - inner_window_origin_x) / inner_window_width;
                        corner_y = figure_height * ((y + 0.5 * c.pts[2 * k + 1]) - inner_window_origin_y) / inner_window_height;
                        k = c.ed[2 * k];
                        
                        if (corner_x < 0) min_xb = -1;
                        if (corner_y < 0) min_yb = -1;
                        if (corner_x > figure_width)  max_xb = +1;
                        if (corner_y > figure_height) max_yb = +1;
                    } while (k != 0);
                }
                
                // IF ANY PART OF THE VORONOI CELL IS IN OUR WINDOW THEN WE SHOULD DRAW IT
                else
                {
                    int k = 0;
                    double corner_x = figure_width * ((x + 0.5 * c.pts[2 * k]) - inner_window_origin_x) / inner_window_width;
                    double corner_y = figure_height * ((y + 0.5 * c.pts[2 * k + 1]) - inner_window_origin_y) / inner_window_height;
                    k = c.ed[2 * k];
                    if (0 <= corner_x && corner_x <= figure_width && 0 <= corner_y && corner_y <= figure_height) to_draw = 1;
                    
                    do {
                        corner_x = figure_width * ((x + 0.5 * c.pts[2 * k]) - inner_window_origin_x) / inner_window_width;
                        corner_y = figure_height * ((y + 0.5 * c.pts[2 * k + 1]) - inner_window_origin_y) / inner_window_height;
                        k = c.ed[2 * k];
                        if (0 <= corner_x && corner_x <= figure_width && 0 <= corner_y && corner_y <= figure_height) to_draw = 1;
                    } while (k != 0);
                }
                
                // DO NOT DRAW VORONOI CELLS THAT DO NOT INTERSECT OUR REGION.
                if (to_draw == 0) continue;
                

                // WHEN THE INNER_WINDOW IS SMALLER THAN THE SYSTEM, WE JUST
                for (int s = min_xb; s <= max_xb; ++s) for (int t = min_yb; t <= max_yb; ++t)
                {
                    output_file << "newpath\n";
                    
                    int k = 0;
                    double corner_x = figure_width * ((x + 0.5 * c.pts[2 * k]) - inner_window_origin_x) / inner_window_width;
                    double corner_y = figure_height * ((y + 0.5 * c.pts[2 * k + 1]) - inner_window_origin_y) / inner_window_height;
                    output_file << corner_x - s * figure_width << '\t' << corner_y - t * figure_height << " moveto\n";
                    k = c.ed[2 * k];
                    
                    do {
                        corner_x = figure_width * ((x + 0.5 * c.pts[2 * k]) - inner_window_origin_x) / inner_window_width;
                        corner_y = figure_height * ((y + 0.5 * c.pts[2 * k + 1]) - inner_window_origin_y) / inner_window_height;
                        output_file << corner_x - s * figure_width << '\t' << corner_y - t * figure_height << " lineto\n";
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
        
        double x_in_eps_units = figure_width * (x - inner_window_origin_x) / inner_window_width;
        double y_in_eps_units = figure_height * (y - inner_window_origin_y) / inner_window_height;
        
        int min_xb = 0;
        int max_xb = 0;
        int min_yb = 0;
        int max_yb = 0;
        
        // DO NOT DRAW PARTICLES THAT DO NOT INTERSECT THE INNER WINDOW.
        if (particles_in_eps != number_of_particles)
        {
            if (x_in_eps_units < -particle_radius)                continue;
            if (y_in_eps_units < -particle_radius)                continue;
            if (x_in_eps_units > figure_width + particle_radius) continue;
            if (y_in_eps_units > figure_height + particle_radius) continue;
        }
        
        // WHEN DRAWING THE ENTIRE SYSTEM WITH PERIODIC BOUNDARY CONDITIONS, THERE MAY
        // BE PARTICLES AT AN EDGE OF THE SYSTEM THAT SHOULD ALSO BE DRAWN ON THE
        // OPPOSITE SIDE, AS THE PARTICLE MAY INTERSECT BOTH EDGES.  IN SOME CASES THEY
        // SHOULD BE DRAWN FOUR TIMES.
        else
        {
            if (x_in_eps_units < particle_radius)                 min_xb = -1;
            if (y_in_eps_units < particle_radius)                 min_yb = -1;
            if (x_in_eps_units > figure_width - particle_radius) max_xb = +1;
            if (y_in_eps_units > figure_height - particle_radius) max_yb = +1;
        }
        
        int voronoi_cell_sides = cell_neighbor_count[pid];
        
        // IN THIS LOOP EACH PARTICLE IS DRAWN ONCE, TWICE, OR FOUR TIMES
        for (int s = min_xb; s <= max_xb; ++s) for (int t = min_yb; t <= max_yb; ++t)
        {
            // DRAW PARTICLES ALL IN BLACK
            if (particle_coloring_scheme == 1)
                output_file << x_in_eps_units - s * figure_width << " " << y_in_eps_units - t * figure_height << " " << particle_radius << " 0 360 arc fill\n";
            
            // COLOR PARTICLES ACCORDING TO NUMBER OF EDGES
            else if (particle_coloring_scheme == 2)
            {
                output_file << color_strings[voronoi_cell_sides] << " setrgbcolor\n";
                output_file << x_in_eps_units - s * figure_width << " " << y_in_eps_units - t * figure_height << " " << particle_radius << " 0 360 arc fill\n";
                output_file << "0 0 0 setrgbcolor\n";
                output_file << x_in_eps_units - s * figure_width << " " << y_in_eps_units - t * figure_height << " " << particle_radius << " 0 360 arc stroke\n";
            }
            
            // COLOR PARTICLES ACCORDING TO FILTER INDEXES
            else if (particle_coloring_scheme == 3)
            {
                output_file << color_strings[vt_structure_types[pid]] << " setrgbcolor\n";
                output_file << x_in_eps_units - s * figure_width << " " << y_in_eps_units - t * figure_height << " " << particle_radius << " 0 360 arc fill\n";
                output_file << "0 0 0 setrgbcolor\n";
                output_file << x_in_eps_units - s * figure_width << " " << y_in_eps_units - t * figure_height << " " << particle_radius << " 0 360 arc stroke\n";
            }
            
            // COLOR PARTICLES ACCORDING TO VORONOI DISTANCE FROM CENTRAL PARTICLE
            else if (particle_coloring_scheme == 4)
            {
                output_file << color_strings[ring_index[pid]] << " setrgbcolor\n";
                output_file << x_in_eps_units - s * figure_width << " " << y_in_eps_units - t * figure_height << " " << particle_radius << " 0 360 arc fill\n";
                output_file << "0 0 0 setrgbcolor\n";
                output_file << x_in_eps_units - s * figure_width << " " << y_in_eps_units - t * figure_height << " " << particle_radius << " 0 360 arc stroke\n";
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
                output_file << x_in_eps_units - s * figure_width << " " << y_in_eps_units - t * figure_height << " " << particle_radius << " 0 360 arc fill\n";
                output_file << "0 0 0 setrgbcolor\n";
                output_file << x_in_eps_units - s * figure_width << " " << y_in_eps_units - t * figure_height << " " << particle_radius << " 0 360 arc stroke\n";
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
    double closest_distance_squared_to_center = supercell_edges[0][0] + supercell_edges[1][1];
    int closest_particle_ID = 0;
    
    double system_width = hi_bound[0] - origin[0];
    double system_height = hi_bound[1] - origin[1];

    double cell_centerx = origin[0] + system_width / 2.;
    double cell_centery = origin[1] + system_height / 2.;
    
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
                if (visited[neighbors_list_char[tempc][d]] == -1)
                {
                    ring_index[neighbors_list_char[tempc][d]] = k;
                    visited[neighbors_list_char[tempc][d]] = k;
                    cluster.push_back(neighbors_list_char[tempc][d]);
                    kneighbors[k]++;
                }
            }
        }
    }
}
