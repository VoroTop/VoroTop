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

////   File: vectors.cc


#include <array>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>

#include "filters.hh"
#include "variables.hh"
#include "functions.hh"

using namespace voro;



////////////////////////////////////////////////////
////
////   COMPUTES W-VECTORS FOR ALL PARTICLES;
////   REPORTS NUMBER THAT ARE "BROKEN".
////
////////////////////////////////////////////////////

void count_and_store_neighbors_2d(container_2d& con)
{
#pragma omp parallel for num_threads(threads)
    for(auto cli=con.begin();cli<con.end();cli++)
    {
        voronoicell_neighbor_2d c(con);
        if (con.compute_cell(c,cli))
        {
            const int ijk=cli->ijk,q=cli->q;
            const int pid = con.id[ijk][q];
            cell_neighbor_count[pid]=c.p;
            c.neighbors(list_of_neighbors[pid]);
        }
        else
        {
            std::cout << "A Voronoi cell failed to compute\n";
            exit(0);
        }
    }
}


template<class Container>
void count_and_store_neighbors_3d(Container& con)
{
#pragma omp parallel for num_threads(threads)
    for(auto cli=con.begin();cli<con.end();cli++)
    {
        voronoicell_neighbor_3d c(con);
        if (con.compute_cell(c,cli))
        {
            const int ijk=cli->ijk,q=cli->q;
            const int pid = con.id[ijk][q];
            c.neighbors(list_of_neighbors[pid]);
            cell_neighbor_count[pid]=list_of_neighbors[pid].size();
        }
        else
        {
            std::cout << "A Voronoi cell failed to compute\n";
            exit(0);
        }
    }
}
template void count_and_store_neighbors_3d(container_3d&);
template void count_and_store_neighbors_3d(container_triclinic&);


void calc_distribution_2d(Filter &filter)
{
    std::vector<Filter> local_filter(threads);
    
    double system_width  = xhi-xlo;
    double system_height = yhi-ylo;

#pragma omp parallel for num_threads(threads)
    for (int pid=0; pid<number_of_particles; pid++)
    {
        unsigned int number_of_neighbors = cell_neighbor_count[pid];
        std::vector< std::pair <double,int> > unordered_neighbors;
        unordered_neighbors.reserve(number_of_neighbors);
        
        for(unsigned int q=0; q<number_of_neighbors; q++)
        {
            double dx = particle_coordinates[2*list_of_neighbors[pid][q]]   - particle_coordinates[2*pid];
            double dy = particle_coordinates[2*list_of_neighbors[pid][q]+1] - particle_coordinates[2*pid+1];
            if(dx >  system_width/2.)  dx -= system_width;
            if(dx < -system_width/2.)  dx += system_width;
            if(dy >  system_height/2.) dy -= system_height;
            if(dy < -system_height/2.) dy += system_height;
            double theta = std::atan2(dy, dx);
            
            int p_id = cell_neighbor_count[list_of_neighbors[pid][q]];
            unordered_neighbors.emplace_back(theta, p_id);
        }
        std::sort(unordered_neighbors.begin(), unordered_neighbors.end());
        
        std::vector<std::vector<int>> numbers1(number_of_neighbors, std::vector<int>(number_of_neighbors));
        std::vector<std::vector<int>> numbers2(number_of_neighbors, std::vector<int>(number_of_neighbors));

        for(unsigned int c=0; c<number_of_neighbors; c++)
        {
            for(unsigned int d=0; d<number_of_neighbors; d++)
            {
                numbers1[c][d] = unordered_neighbors[(c+d)%number_of_neighbors].second;
                numbers2[c][d] = unordered_neighbors[(c-d+number_of_neighbors)%number_of_neighbors].second;
            }
        }
        
        std::sort(numbers1.begin(), numbers1.end());
        std::sort(numbers2.begin(), numbers2.end());
        
        int                                chirality =  0;
        if     (numbers1[0] < numbers2[0]) chirality =  1;
        else if(numbers1[0] > numbers2[0]) chirality = -1;
        
        std::vector<int> canonical_code;  // CANONICAL CODE WILL BE STORED HERE
        if(chirality==-1)
        {
            canonical_code=numbers2[0];
            canonical_code.insert(canonical_code.begin(), number_of_neighbors);
        }
        else
        {
            canonical_code=numbers1[0];
            canonical_code.insert(canonical_code.begin(), number_of_neighbors);
        }
        
        const int tid=omp_get_thread_num();
        local_filter[tid].increment_or_add(canonical_code,chirality,1);
    }
    
    for(int tid=0; tid<threads; tid++)
        filter.copy_filter(local_filter[tid]);    
}


template<class Container>
void calc_distribution_3d(Container& con, Filter &filter)
{
    std::vector<Filter> local_filter(threads);

#pragma omp parallel for num_threads(threads)
    for(auto cli=con.begin();cli<con.end();cli++)
    {
        voronoicell_3d vcell;
        if (con.compute_cell(vcell,cli))
        {
            VoronoiTopology result = compute_canonical_code_3d(vcell);

            const int tid=omp_get_thread_num();
            local_filter[tid].increment_or_add(result.canonical_code,result.chirality,1);
        }
    }

    for(int tid=0; tid<threads; tid++)
        filter.copy_filter(local_filter[tid]);
}
template void calc_distribution_3d(container_3d&, Filter&);
template void calc_distribution_3d(container_triclinic&, Filter&);



void  print_topology_vectors_2d(std::string filename)
{
    constexpr int MAX_NEIGHBORS = 50;
    std::string vectors_name(filename);
    vectors_name.append(".vectors");
    std::ofstream vector_file(vectors_name.c_str(), std::ofstream::out);
    
    // BUFFER OUTPUT FROM EACH THREAD
    std::vector<stringstream> buf(threads);

    double system_width  = xhi-xlo;
    double system_height = yhi-ylo;
    
#pragma omp parallel for num_threads(threads)
    for (int pid=0; pid<number_of_particles; pid++) 
    {
        std::array<int, MAX_NEIGHBORS> p2pvector = {0};
        int maxp = 0;                               // STORE THE MAX NUMBER OF EDGES AMONG NEIGHBORS
        
        unsigned int number_of_neighbors = cell_neighbor_count[pid];
        std::vector< std::pair <double,int> > unordered_neighbors;
        unordered_neighbors.reserve(number_of_neighbors);
        
        for(unsigned int q=0; q<number_of_neighbors; q++)
        {
            p2pvector[cell_neighbor_count[list_of_neighbors[pid][q]]]++;
            if(cell_neighbor_count[list_of_neighbors[pid][q]]>maxp) maxp = cell_neighbor_count[list_of_neighbors[pid][q]];
            
            double dx = particle_coordinates[2*list_of_neighbors[pid][q]]   - particle_coordinates[2*pid];
            double dy = particle_coordinates[2*list_of_neighbors[pid][q]+1] - particle_coordinates[2*pid+1];
            if(dx >  system_width/2.)  dx -= system_width;
            if(dx < -system_width/2.)  dx += system_width;
            if(dy >  system_height/2.) dy -= system_height;
            if(dy < -system_height/2.) dy += system_height;
            double theta = std::atan2(dy, dx);

            int p_id = cell_neighbor_count[list_of_neighbors[pid][q]];
            unordered_neighbors.emplace_back(theta, p_id);
        }
        std::sort(unordered_neighbors.begin(), unordered_neighbors.end());
        
        std::vector<std::vector<int>> numbers1(number_of_neighbors, std::vector<int>(number_of_neighbors));
        std::vector<std::vector<int>> numbers2(number_of_neighbors, std::vector<int>(number_of_neighbors));

        for(unsigned int c=0; c<number_of_neighbors; c++)
        {
            for(unsigned int d=0; d<number_of_neighbors; d++)
            {
                numbers1[c][d] = unordered_neighbors[(c+d)%number_of_neighbors].second;
                numbers2[c][d] = unordered_neighbors[(c-d+number_of_neighbors)%number_of_neighbors].second;
            }
        }
        
        std::sort(numbers1.begin(), numbers1.end());
        std::sort(numbers2.begin(), numbers2.end());

        int                                chirality =  0;
        if     (numbers1[0] < numbers2[0]) chirality =  1;
        else if(numbers1[0] > numbers2[0]) chirality = -1;
        
        int symmetry_order = 1;
        for(unsigned int c=1; c<number_of_neighbors; c++)
            if(numbers1[c] == numbers1[0])
                symmetry_order++;
        if(chirality==0) symmetry_order *= 2;

        const int tid=omp_get_thread_num();
        buf[tid] << particle_ids[pid] << '\t';   // PARTICLE ID
        buf[tid] << number_of_neighbors        << '\t';   // NUMBER OF EDGES
        buf[tid] << '(';                         // NUMBER OF NEIGHBORS WITH DIFFERENT NUMBRERS OF EDGES
        for(int d=3; d<maxp; d++)
            buf[tid] << p2pvector[d] << ",";
        buf[tid] << p2pvector[maxp] << ")\t";
        
        if(chirality==-1)
        {
            buf[tid] << "(";
            buf[tid] << number_of_neighbors << ",";
            for(long unsigned int d=0; d<numbers2[0].size()-1; d++)
                buf[tid] << numbers2[0][d] << ",";
            buf[tid] << numbers2[0][numbers2[0].size()-1] << ")\t";
        }
        else
        {
            buf[tid] << "(";
            buf[tid] << number_of_neighbors << ",";
            for(long unsigned int d=0; d<numbers1[0].size()-1; d++)
                buf[tid] << numbers1[0][d] << ",";
            buf[tid] << numbers1[0][numbers1[0].size()-1] << ")\t";
        }
        buf[tid] << symmetry_order << '\t';
        buf[tid] << chirality << '\n';
    }
    
    // OUTPUT ALL DATA TO FILE
    for(int t=0; t<threads; t++)
        vector_file << buf[t].rdbuf();
    
    vector_file.close();
}




template<class Container>
void print_topology_vectors_3d(Container& con, std::string filename)
{
    std::string vectors_name(filename);
    vectors_name.append(".vectors");
    std::ofstream vector_file(vectors_name.c_str(), std::ofstream::out);

    std::vector<stringstream> buf(threads);

#pragma omp parallel for num_threads(threads)
    for(auto cli=con.begin();cli<con.end();cli++)
    {
        voronoicell_3d vcell;
        if (con.compute_cell(vcell,cli))
        {
            int ijk=cli->ijk,q=cli->q;
            int pid = con.id[ijk][q];
            
            // COMPUTE THE CANONICAL CODE, CHIRALITY, SYMMETRY, AND FACE DATA
            VoronoiTopology result = compute_canonical_code_3d(vcell);

            int edge_count       = vcell.number_of_edges();
            int face_count       = result.face_count;
            int max_face_edges   = result.max_face_edges;
            int chirality        = result.chirality;
            int symmetry_counter = result.symmetry_counter;
            
            
            const int tid=omp_get_thread_num();

            buf[tid] << particle_ids[pid] << '\t';     // PARTICLE ID
            buf[tid] << face_count        << '\t';     // NUMBER OF FACES
            buf[tid] << '(';                           // P VECTOR
            for(int d=3; d<max_face_edges; d++)
                buf[tid] << result.face_edge_counts[d] << ',';
            buf[tid] << result.face_edge_counts[max_face_edges] << ')' << '\t';
            
            buf[tid] << '(';                           // WEINBERG VECTOR
            for(int d=0; d<2*edge_count; d++)
                buf[tid] << result.canonical_code[d] << ',';
            buf[tid] << 1 << ')' << '\t';
            
            buf[tid] << symmetry_counter << '\t';      // SYMMETRIES
            buf[tid] << chirality << '\t';             // CHIRALITY
            
            buf[tid] << '\n';
        }
    }
    
    // OUTPUT ALL DATA TO FILE
    for(int t=0; t<threads; t++)
        vector_file << buf[t].rdbuf();
    
    vector_file.close();
}
template void print_topology_vectors_3d(container_3d&, std::string);
template void print_topology_vectors_3d(container_triclinic&, std::string);


VoronoiTopology compute_canonical_code_3d(voro::voronoicell_3d& vcell)
{
    const int max_epf = 256;    // MAXIMUM EDGES PER FACE
    const int max_epc = 512;    // MAXIMUM EDGES PER CELL
    const int max_vpc = 512;    // MAXIMUM VERTICES PER CELL

    int   edge_count     = vcell.number_of_edges();
    int   vertex_count   = vcell.p;    // TOTAL NUMBER OF VERTICES
    int*  vertex_degrees = vcell.nu;   // VERTEX DEGREE ARRAY
    int** ed             = vcell.ed;   // EDGE CONNECTIONS ARRAY

    std::vector<int> canonical_code(2*edge_count, 0);  // CANONICAL CODE WILL BE STORED HERE

    int face_count         = 0;
    int max_face_edges     = 3;     // EVERY CONVEX POLYHEDRON MUST HAVE AT LEAST ONE FACE WITH 3 OR MORE EDGES
    int min_face_edges     = 5;     // EVERY CONVEX POLYHEDRON MUST HAVE AT LEAST ONE FACE WITH 5 OR FEWER EDGES
    int face_edge_counts[max_epf] = {};  // NUMBER OF FACES WITH EACH NUMBER OF EDGES
    int origins[2*max_epc] = {};    // NO VORONOI CELL IN FILTER HAS MORE THAN max_epc EDGES
    int origin_c           = 0;

    // USE VORO++'S FACE_VERTICES() TO GET FACE DATA DIRECTLY,
    // AVOIDING MANUAL FACE TRACING AND THE SIGN-FLIP/RESET CYCLE.
    std::vector<int> fv;
    vcell.face_vertices(fv);
    {
        int pos = 0;
        while(pos < (int)fv.size())
        {
            int face_c = fv[pos];   // NUMBER OF EDGES/VERTICES IN THIS FACE
            face_count++;
            face_edge_counts[face_c]++;
            if(face_c > max_face_edges)
                max_face_edges = face_c;
            if(face_c < min_face_edges)
            {
                min_face_edges = face_c;
                origin_c = face_c;
                for(int c = 0; c < face_c; c++)
                    origins[c] = fv[pos + 1 + c];
            }
            else if(face_c == min_face_edges)
            {
                for(int c = 0; c < face_c; c++)
                    origins[origin_c + c] = fv[pos + 1 + c];
                origin_c += face_c;
            }
            pos += face_c + 1;
        }
    }

    // KEEPING TRACK OF THIS WILL ALLOW US TO SPEED UP SOME COMPUTATION, OF BCC
    int likely_bcc=0;
    if(face_count==14 && face_edge_counts[4]==6 && face_edge_counts[6]==8) likely_bcc=1;   // FACE_EDGE_COUNTS (0,6,0,8,0,...) OF A SIMPLE POLYHEDRON APPEARS IN 3 DIFFERENT TYPES, WITH SYMMETRIES 4, 8, AND 48
    
    
    ////////////////////////////////////////////////////////////////
    // BUILD THE CANONICAL CODE
    ////////////////////////////////////////////////////////////////

    int vertices_temp_labels[max_vpc] = {};         // TEMPORARY LABELS FOR ALL VERTICES; MAX max_vpc VERTICES

    // GENERATION COUNTER FOR O(1) RESET OF VISITED-DART STATE.
    // INSTEAD OF FLIPPING SIGNS IN ed[][] AND RESETTING BEFORE EACH
    // STARTING DART, WE INCREMENT A GENERATION COUNTER.  A DART
    // (vertex, branch) IS "VISITED" IFF visited_gen[offset] == generation.
    int dart_offset[max_vpc];
    dart_offset[0] = 0;
    for(int i=1; i<vertex_count; i++)
        dart_offset[i] = dart_offset[i-1] + vertex_degrees[i-1];
    int visited_gen[2*max_epc] = {};
    int generation = 0;

    int finished   =  0;
    int chirality  = -1;
    int symmetry_counter = 0;     // TRACKS NUMBER OF REPEATS OF A CODE, I.E. SYMMETRY ORDER

    for(int orientation=0; orientation<2 && finished==0; orientation++)
    {
        for(int q=0; q<origin_c && finished==0; q++)
        {
            // CLEAR ALL LABELS; MARK ALL BRANCHES AS NEW (O(1) RESET)
            std::fill(vertices_temp_labels, vertices_temp_labels+vertex_count, 0);
            generation++;

            int initial = origins[q];
            int next;
            int branch=0;

            if(orientation==0)
            {
                if((q+1)%min_face_edges==0) next = origins[q - min_face_edges + 1];
                else next = origins[q + 1];
            }
            else
            {
                if(q    %min_face_edges==0) next = origins[q + min_face_edges - 1];
                else next = origins[q - 1];
            }
            for(int j=0; j<vertex_degrees[origins[q]]; j++)
                if(ed[origins[q]][j]==next) branch=j;
            visited_gen[dart_offset[initial] + branch] = generation;

            int current_code_length   = 0;
            int current_highest_label = 1;
            int continue_code         = 0;    // 0: UNDECIDED; 1: GO AHEAD, DO NOT EVEN CHECK.
            if(q==0 && orientation==0)        // FIRST CODE, GO AHEAD
                continue_code=1;

            vertices_temp_labels[initial] = current_highest_label++;
            canonical_code[current_code_length]  = vertices_temp_labels[initial];
            current_code_length++;

            // BUILD EACH CODE FOLLOWING WEINBERG'S RULES FOR TRAVERSING A GRAPH TO BUILD
            // A HAMILTONIAN PATH, LABELING VERTICES ALONG THE WAY, AND RECORDING VERTICES
            // AS VISITED.
            int end_flag=0;
            while(end_flag==0)
            {
                // NEXT VERTEX HAS NOT BEEN VISITED; TAKE RIGHT-MOST BRANCH TO CONTINUE.
                if(vertices_temp_labels[next]==0)
                {
                    // LABEL THE NEW VERTEX
                    vertices_temp_labels[next] = current_highest_label++;

                    if(continue_code==0)
                    {
                        if(vertices_temp_labels[next]>canonical_code[current_code_length]) break;
                        if(vertices_temp_labels[next]<canonical_code[current_code_length])
                        {
                            symmetry_counter = 0;
                            continue_code    = 1;
                            if(orientation==1) chirality=1;
                        }
                    }

                    // BUILD THE CODE
                    canonical_code[current_code_length] = vertices_temp_labels[next];
                    current_code_length++;

                    // FIND NEXT DIRECTION TO MOVE ALONG, UPDATE, AND RELOOP
                    if(orientation==0) branch  = vcell.cycle_up  (ed[initial][vertex_degrees[initial]+branch],next);
                    else               branch  = vcell.cycle_down(ed[initial][vertex_degrees[initial]+branch],next);
                    initial = next;
                    next    = ed[initial][branch];
                    visited_gen[dart_offset[initial] + branch] = generation;
                }

                else    // NEXT VERTEX *HAS* BEEN VISITED BEFORE
                {
                    int next_branch = ed[initial][vertex_degrees[initial]+branch];
                    int branches_tested = 0;

                    while(visited_gen[dart_offset[next] + next_branch] == generation && branches_tested<vertex_degrees[next])
                    {
                        if(orientation==0) next_branch = vcell.cycle_up  (next_branch,next);
                        else               next_branch = vcell.cycle_down(next_branch,next);

                        branches_tested++;
                    }

                    if(branches_tested < vertex_degrees[next])
                    {
                        if(continue_code==0)
                        {
                            if(vertices_temp_labels[next]>canonical_code[current_code_length]) break;
                            if(vertices_temp_labels[next]<canonical_code[current_code_length])
                            {
                                symmetry_counter = 0;
                                continue_code    = 1;
                                if(orientation==1) chirality=1;
                            }
                        }

                        // BUILD THE CODE
                        canonical_code[current_code_length] = vertices_temp_labels[next];
                        current_code_length++;

                        // FIND NEXT BRANCH
                        branch  = next_branch;
                        initial = next;
                        next    = ed[initial][branch];
                        visited_gen[dart_offset[initial] + branch] = generation;
                    }

                    else
                    {
                        end_flag=1;

                        if(likely_bcc && symmetry_counter>4 && orientation==0) { chirality=0; symmetry_counter = 48; finished=1; }
                        else if(chirality==-1 && orientation==1)               { chirality=0; symmetry_counter *= 2; finished=1; }
                        else symmetry_counter++;
                    }
                }
            }
        }
    }
    canonical_code.push_back(1);
    return {canonical_code, chirality, symmetry_counter,
            face_count, max_face_edges,
            std::vector<int>(face_edge_counts, face_edge_counts + max_epf)};
}


template<class Container>
int classify_particles_by_voronoi_topology_3d(Container& con, Filter &filter)
{
#pragma omp parallel for num_threads(threads)
    for(auto cli=con.begin();cli<con.end();cli++)
    {
        voronoicell_3d vcell;
        if (con.compute_cell(vcell,cli))
        {
            const int ijk=cli->ijk,q=cli->q;
            const int pid = con.id[ijk][q];

            VoronoiTopology result = compute_canonical_code_3d(vcell);

            vt_structure_types[pid] = filter.vt_structure_type(result.canonical_code);
        }
    }

    return 0;
}
template int classify_particles_by_voronoi_topology_3d(container_3d&, Filter&);
template int classify_particles_by_voronoi_topology_3d(container_triclinic&, Filter&);


void classify_particles_by_voronoi_topology_2d(Filter &filter)
{
    double system_width  = xhi-xlo;
    double system_height = yhi-ylo;
        
#pragma omp parallel for num_threads(threads)
    for (int pid=0; pid<number_of_particles; pid++)
    {
        unsigned int number_of_neighbors = cell_neighbor_count[pid];
        std::vector< std::pair <double,int> > unordered_neighbors;
        unordered_neighbors.reserve(number_of_neighbors);
        
        for(unsigned int q=0; q<number_of_neighbors; q++)
        {
            double dx = particle_coordinates[2*list_of_neighbors[pid][q]]   - particle_coordinates[2*pid];
            double dy = particle_coordinates[2*list_of_neighbors[pid][q]+1] - particle_coordinates[2*pid+1];
            if(dx >  system_width/2.)  dx -= system_width;
            if(dx < -system_width/2.)  dx += system_width;
            if(dy >  system_height/2.) dy -= system_height;
            if(dy < -system_height/2.) dy += system_height;
            double theta = std::atan2(dy, dx);

            int p_id = cell_neighbor_count[list_of_neighbors[pid][q]];
            unordered_neighbors.emplace_back(theta, p_id);
        }
        std::sort(unordered_neighbors.begin(), unordered_neighbors.end());
        
        // CONSTRUCT THE NON-CANONICAL CODES
        std::vector<std::vector<int>> non_canonical_codes(2*number_of_neighbors, std::vector<int>(number_of_neighbors));
        for(unsigned int c=0; c<number_of_neighbors; c++)
        {
            for(unsigned int d=0; d<number_of_neighbors; d++)
            {
                non_canonical_codes[c][d] = unordered_neighbors[(c+d)%number_of_neighbors].second;
                non_canonical_codes[c+number_of_neighbors][d] = unordered_neighbors[(c-d+number_of_neighbors)%number_of_neighbors].second;
            }
        }        

        std::sort(non_canonical_codes.begin(), non_canonical_codes.end());
        non_canonical_codes[0].insert(non_canonical_codes[0].begin(), number_of_neighbors);
        vt_structure_types[pid] = filter.vt_structure_type(non_canonical_codes[0]);
    }
}

