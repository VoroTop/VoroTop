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

////   File: vectors.cc


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

void count_and_store_neighbors(container_2d& con)
{
#pragma omp parallel for num_threads(threads)
    for(container_2d::iterator cli=con.begin();cli<con.end();cli++)
    {
        voronoicell_neighbor_2d c(con);
        if (con.compute_cell(c,cli))
        {
            int ijk=cli->ijk,q=cli->q;
            int pid = con.id[ijk][q];
            cell_neighbor_count[pid]=c.p;
            c.neighbors(neighbors_list_char[pid]);
        }
        else
        {
            std::cout << "A Voronoi cell failed to compute\n";
            exit(0);
        }
    }
}


void count_and_store_neighbors(container_3d& con)
{
#pragma omp parallel for num_threads(threads)
    for(container_3d::iterator cli=con.begin();cli<con.end();cli++)
    {
        voronoicell_neighbor_3d c(con);
        if (con.compute_cell(c,cli))
        {
            int ijk=cli->ijk,q=cli->q;
            int pid = con.id[ijk][q];
            c.neighbors(neighbors_list_char[pid]);
            cell_neighbor_count[pid]=neighbors_list_char[pid].size();
        }
        else
        {
            std::cout << "A Voronoi cell failed to compute\n";
            exit(0);
        }
    }
}


void calc_distribution(container_2d& con, Filter &filter)
{
    Filter local_filter[threads];
    
    double xdim = hi_bound[0]-origin[0];
    double ydim = hi_bound[1]-origin[1];
    
#pragma omp parallel for num_threads(threads)
    for (int pid=0; pid<number_of_particles; pid++)
    {
        // THIS IS OUR LIST OF NEIGHBORS, IT IS A VECTOR.
        int p2pvector[50];
        for(int d=0; d<50; d++) p2pvector[d]=0;
        int maxp = 0;                               // STORE THE MAX NUMBER OF EDGES AMONG NEIGHBORS
        
        unsigned int nneighbors = neighbors_list_char[pid].size();
        std::vector< std::pair <double,int> > prepvect;
        prepvect.reserve(nneighbors);
        
        for(unsigned int q=0; q<nneighbors; q++)
        {
            p2pvector[cell_neighbor_count[neighbors_list_char[pid][q]]]++;
            if(cell_neighbor_count[neighbors_list_char[pid][q]]>maxp) maxp = cell_neighbor_count[neighbors_list_char[pid][q]];
            
            double dx = particle_coordinates[2*neighbors_list_char[pid][q]]   - particle_coordinates[2*pid];
            double dy = particle_coordinates[2*neighbors_list_char[pid][q]+1] - particle_coordinates[2*pid+1];
            if(dx >  xdim/2) dx -= xdim;
            if(dx < -xdim/2) dx += xdim;
            if(dy >  ydim/2) dy -= ydim;
            if(dy < -ydim/2) dy += ydim;
            double theta1 = atan(dy/dx) + 3.14159265359/2.;
            if(dx<0)            theta1 += 3.14159265359;
            
            int p_id = cell_neighbor_count[neighbors_list_char[pid][q]];
            prepvect.push_back(std::make_pair(theta1, p_id));
        }
        std::sort(prepvect.begin(), prepvect.end());
        
        std::vector< std::vector <int> > numbers1;
        std::vector< std::vector <int> > numbers2;
        
        for(unsigned int c=0; c<nneighbors; c++)
        {
            std::vector<int> v;
            for(unsigned int d=0; d<nneighbors; d++)
                v.push_back(prepvect[(c+d)%nneighbors].second);
            numbers1.push_back(v);
        }
        for(unsigned int c=0; c<nneighbors; c++)
        {
            std::vector<int> v;
            for(unsigned int d=0; d<nneighbors; d++)
                v.push_back(prepvect[(c-d+nneighbors)%nneighbors].second);
            numbers2.push_back(v);
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
            canonical_code.insert(canonical_code.begin(), nneighbors);
        }
        else
        {
            canonical_code=numbers1[0];
            canonical_code.insert(canonical_code.begin(), nneighbors);
        }
        
        int tid=omp_get_thread_num();
        local_filter[tid].increment_or_add(canonical_code,chirality,1);
    }
    
    for(int tid=0; tid<threads; tid++)
        filter.copy_filter(local_filter[tid]);
    
    return;
}


void calc_distribution(container_3d& con, Filter &filter)
{
    Filter local_filter[threads];
    
#pragma omp parallel for num_threads(threads)
    for(container_3d::iterator cli=con.begin();cli<con.end();cli++)
    {
        voronoicell_3d vcell;
        if (con.compute_cell(vcell,cli))
        {
            std::vector<int> canonical_code;                // CANONICAL CODE WILL BE STORED HERE
            int chirality = compute_canonical_code(canonical_code, vcell);
            
            int tid=omp_get_thread_num();
            local_filter[tid].increment_or_add(canonical_code,chirality,1);
        }
    }
    
    for(int tid=0; tid<threads; tid++)
        filter.copy_filter(local_filter[tid]);
}






////////////////////////////////////////////////////
////
////    CALCULATES THE NUMBER OF SIDES AND THE LIST
////    OF NEIGHBORS OF EACH PARTICLE. NECESSARY
////    FOR CALCULATING P-VECTORS AND FOR SOME DRAWING
////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
////
////    CALCULATES THE P-VECTOR FOR EACH PARTICLE.
////    REQUIRES PREVIOUSLY CALLING calc_neighbors()
////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
////
////    SHORTENS EACH P-VECTOR TO NOT INCLUDE THE
////    SYMMETRY ORDER OR CHIRALITY
////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
////
////    CALCULATES THE P-VECTOR FOR A SINGLE PARTICLE.
////    USES INFORMATION COMPUTED IN calc_neighbors(),
////    WHICH COMPUTES THE NUMBER OF NEIGHBORS OF EACH
////    PARTICLE, AND A LIST OF THOSE NEIGHBORS.
////
////////////////////////////////////////////////////




void  print_topology_vectors2d(std::string filename)
{
    std::string vectors_name(filename);
    vectors_name.append(".vectors");
    std::ofstream vector_file(vectors_name.c_str(), std::ofstream::out);
    
    // BUFFER OUTPUT FROM EACH THREAD
    std::stringstream buf[threads];
    
    double xdim = hi_bound[0]-origin[0];
    double ydim = hi_bound[1]-origin[1];
    
    
#pragma omp parallel for num_threads(threads)
    for (int pid=0; pid<number_of_particles; pid++) {
        {
            // THIS IS OUR LIST OF NEIGHBORS, IT IS A VECTOR.
            int p2pvector[50];
            for(int d=0; d<50; d++) p2pvector[d]=0;
            int maxp = 0;                               // STORE THE MAX NUMBER OF EDGES AMONG NEIGHBORS
            
            
            unsigned int nneighbors = neighbors_list_char[pid].size();
            std::vector< std::pair <double,int> > prepvect;
            prepvect.reserve(nneighbors);
            
            for(unsigned int q=0; q<nneighbors; q++)
            {
                p2pvector[cell_neighbor_count[neighbors_list_char[pid][q]]]++;
                if(cell_neighbor_count[neighbors_list_char[pid][q]]>maxp) maxp = cell_neighbor_count[neighbors_list_char[pid][q]];
                
                double dx = particle_coordinates[2*neighbors_list_char[pid][q]]   - particle_coordinates[2*pid];
                double dy = particle_coordinates[2*neighbors_list_char[pid][q]+1] - particle_coordinates[2*pid+1];
                if(dx >  xdim/2) dx -= xdim;
                if(dx < -xdim/2) dx += xdim;
                if(dy >  ydim/2) dy -= ydim;
                if(dy < -ydim/2) dy += ydim;
                double theta1 = atan(dy/dx) + 3.14159265359/2.;
                if(dx<0)            theta1 += 3.14159265359;
                
                int p_id = cell_neighbor_count[neighbors_list_char[pid][q]];
                prepvect.push_back(std::make_pair(theta1, p_id));
            }
            std::sort(prepvect.begin(), prepvect.end());
            
            std::vector< std::vector <int> > numbers1;
            std::vector< std::vector <int> > numbers2;
            
            for(unsigned int c=0; c<nneighbors; c++)
            {
                std::vector<int> v;
                for(unsigned int d=0; d<nneighbors; d++)
                    v.push_back(prepvect[(c+d)%nneighbors].second);
                numbers1.push_back(v);
            }
            for(unsigned int c=0; c<nneighbors; c++)
            {
                std::vector<int> v;
                for(unsigned int d=0; d<nneighbors; d++)
                    v.push_back(prepvect[(c-d+nneighbors)%nneighbors].second);
                numbers2.push_back(v);
            }
            
            std::sort(numbers1.begin(), numbers1.end());
            std::sort(numbers2.begin(), numbers2.end());
            
            int                                chirality =  0;
            if     (numbers1[0] < numbers2[0]) chirality =  1;
            else if(numbers1[0] > numbers2[0]) chirality = -1;
            
            int symmetry_order = 1;
            for(unsigned int c=0; c<nneighbors; c++)    // FIX. THIS SEEMS WRONG. SHOULD INDEX START AT c=1?
                if(numbers1[c] == numbers1[0])
                    symmetry_order = c+1;
            if(chirality==0) symmetry_order *= 2;
            
            
            int tid=omp_get_thread_num();
            buf[tid] << particle_ids[pid] << '\t';   // PARTICLE ID
            buf[tid] << nneighbors        << '\t';   // NUMBER OF EDGES
            buf[tid] << '(';                         // NUMBER OF NEIGHBORS WITH DIFFERENT NUMBRERS OF EDGES
            for(int d=3; d<maxp; d++)
                buf[tid] << p2pvector[d] << ",";
            buf[tid] << p2pvector[maxp] << ")\t";
            
            if(chirality==-1)
            {
                buf[tid] << "(";
                buf[tid] << nneighbors << ",";
                for(long unsigned int d=0; d<numbers2[0].size()-1; d++)
                    buf[tid] << numbers2[0][d] << ",";
                buf[tid] << numbers2[0][numbers2[0].size()-1] << ")\t";
            }
            else
            {
                buf[tid] << "(";
                buf[tid] << nneighbors << ",";
                for(long unsigned int d=0; d<numbers1[0].size()-1; d++)
                    buf[tid] << numbers1[0][d] << ",";
                buf[tid] << numbers1[0][numbers1[0].size()-1] << ")\t";
            }
            buf[tid] << symmetry_order << '\t';
            buf[tid] << chirality << '\n';
            
            //std::vector<std::pair>::iterator it = vector.begin();
            //prepvect
        }
    }
    
    // OUTPUT ALL DATA TO FILE
    for(int t=0; t<threads; t++)
        vector_file << buf[t].rdbuf();
    
    //    std::cout << "Here\n";
    //    while(1==1);
    
    
    vector_file.close();
    return;      // 0 IF ALL CELLS COMPUTED, POSITIVE OTHERWISE
}




int  print_topology_vectors(container_3d& con, std::string filename)
{
    std::string vectors_name(filename);
    vectors_name.append(".vectors");
    std::ofstream vector_file(vectors_name.c_str(), std::ofstream::out);
    
    int counter[threads]; for(int c=0; c<threads; c++) counter[threads]=0;
    std::stringstream buf[threads];
    
#pragma omp parallel for num_threads(threads)
    for(container_3d::iterator cli=con.begin();cli<con.end();cli++)
    {
        voronoicell_3d vcell;
        if (con.compute_cell(vcell,cli))
        {
            int ijk=cli->ijk,q=cli->q;
            int pid = con.id[ijk][q];
            
            const int max_epf = 256;    // MAXIMUM EDGES PER FACE
            const int max_epc = 512;    // MAXIMUM EDGES PER CELL
            const int max_vpc = 512;    // MAXIMUM VERTICES PER CELL
            
            int   edge_count     = vcell.number_of_edges();
            int   vertex_count   = vcell.p;    // TOTAL NUMBER OF VERTICES
            int*  vertex_degrees = vcell.nu;   // VERTEX DEGREE ARRAY
            int** ed             = vcell.ed;   // EDGE CONNECTIONS ARRAY
            
            int face_count         = 0;
            int max_face_edges     = 3;     // EVERY CONVEX POLYHEDRON MUST HAVE AT LEAST ONE FACE WITH 3 OR MORE EDGES
            int min_face_edges     = 5;     // EVERY CONVEX POLYHEDRON MUST HAVE AT LEAST ONE FACE WITH 5 OR FEWER EDGES
            int pvector[max_epf]   = {};    // RECORDS NUMBER OF FACES WITH EACH NUMBER OF EDGES, NO FACE IN FILTER HAS MORE THAN max_epf-1 EDGES
            int origins[2*max_epc] = {};    // NO VORONOI CELL IN FILTER HAS MORE THAN max_epc EDGES
            int origin_c           = 0;
            
            // DETERMINE VERTICES ON FACES WITH MINIMAL EDGES, AND FACES WITH DIFFERENT NUMBERS OF EDGES
            for(int i=0;i<vertex_count;i++)
            {
                for(int j=0;j<vertex_degrees[i];j++)
                {
                    int k = ed[i][j];
                    if(k >= 0)
                    {
                        int face[max_epf]={};  // NO SINGLE FACE WILL HAVE MORE THAN max_epf EDGES
                        int face_c=0;
                        
                        ed[i][j]=-1-k;      // INDICATE THAT WE HAVE CHECKED THIS VERTEX
                        int l=vcell.cycle_up(ed[i][vertex_degrees[i]+j],k);
                        face[face_c++]=k;
                        do {
                            int m=ed[k][l];
                            ed[k][l]=-1-m;
                            l=vcell.cycle_up(ed[k][vertex_degrees[k]+l],m);
                            k=m;
                            
                            face[face_c++]=m;
                        } while (k!=i);
                        
                        // KEEP TRACK OF MINIMAL AND MAXIMAL FACE EDGES
                        if(face_c>max_face_edges)
                            max_face_edges = face_c;
                        if(face_c<min_face_edges)
                        {
                            min_face_edges = origin_c = face_c;
                            for(int c=0; c<face_c; c++)
                                origins[c] = face[c];
                        }
                        else if(face_c==min_face_edges)
                        {
                            for(int c=0; c<face_c; c++)
                                origins[origin_c+c] = face[c];
                            origin_c += face_c;
                        }
                        pvector[face_c]++;
                        face_count++;
                    }
                }
            }
            
            // RESET EDGES
            for(int i=0;i<vertex_count;i++)
                for(int j=0;j<vertex_degrees[i];j++)
                    ed[i][j]=-1-ed[i][j];
            
            // KEEPING TRACK OF THIS WILL ALLOW US TO SPEED UP SOME COMPUTATION, OF BCC
            int likely_bcc=0;
            if(face_count==14 && pvector[4]==6 && pvector[6]==8) likely_bcc=1;   // THIS PVECTOR (0,6,0,8,0,...) OF A SIMPLE POLYHEDRON APPEARS IN 3 DIFFERENT TYPES, WITH SYMMETRIES 4, 8, AND 48
            
            
            ////////////////////////////////////////////////////////////////
            // BUILD THE CANONICAL CODE
            ////////////////////////////////////////////////////////////////
            
            using WeinbergVector = std::vector<int>;
            WeinbergVector canonical_code(2*edge_count,0);  // CANONICAL CODE WILL BE STORED HERE
            int vertices_temp_labels[max_vpc] = {};         // TEMPORARY LABELS FOR ALL VERTICES; MAX max_vpc VERTICES
            
            int finished   =  0;
            int chirality  = -1;
            int symmetry_counter = 0;     // TRACKS NUMBER OF REPEATS OF A CODE, I.E. SYMMETRY ORDER
            
            for(int orientation=0; orientation<2 && finished==0; orientation++)
            {
                for(int q=0; q<origin_c && finished==0; q++)
                {
                    // CLEAR ALL LABELS; MARK ALL BRANCHES OF ALL VERTICES AS NEW
                    std::fill(vertices_temp_labels, vertices_temp_labels+vertex_count, 0);
                    
                    for(int i=0;i<vertex_count;i++)
                        for(int j=0;j<vertex_degrees[i];j++)
                            if(ed[i][j]<0) ed[i][j]=-1-ed[i][j];
                    
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
                    ed[initial][branch] = -1-next;
                    
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
                            ed[initial][branch] = -1-next;
                        }
                        
                        else    // NEXT VERTEX *HAS* BEEN VISITED BEFORE
                        {
                            int next_branch = ed[initial][vertex_degrees[initial]+branch];
                            int branches_tested = 0;
                            
                            while(ed[next][next_branch] < 0 && branches_tested<vertex_degrees[next])
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
                                ed[initial][branch] = -1-next;
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
            
            
            int tid=omp_get_thread_num();
            counter[tid]++;
            
            buf[tid] << particle_ids[pid] << '\t';     // PARTICLE ID
            buf[tid] << face_count        << '\t';     // NUMBER OF FACES
            buf[tid] << '(';                           // P VECTOR
            for(int d=3; d<max_face_edges; d++)
                buf[tid] << pvector[d] << ',';
            buf[tid] << pvector[max_face_edges] << ')' << '\t';
            
            buf[tid] << '(';                           // WEINBERG VECTOR
            for(int d=0; d<2*edge_count; d++)
                buf[tid] << canonical_code[d] << ',';
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
    return 0;
}


int compute_canonical_code(vector<int>& canonical_code, voro::voronoicell_3d& vcell)
{
    const int max_epf = 256;    // MAXIMUM EDGES PER FACE
    const int max_epc = 512;    // MAXIMUM EDGES PER CELL
    const int max_vpc = 512;    // MAXIMUM VERTICES PER CELL
    
    int   edge_count     = vcell.number_of_edges();
    int   vertex_count   = vcell.p;    // TOTAL NUMBER OF VERTICES
    int*  vertex_degrees = vcell.nu;   // VERTEX DEGREE ARRAY
    int** ed             = vcell.ed;   // EDGE CONNECTIONS ARRAY
    
    canonical_code.resize(2*edge_count);  // CANONICAL CODE WILL BE STORED HERE
    
    int face_count         = 0;
    int max_face_edges     = 3;     // EVERY CONVEX POLYHEDRON MUST HAVE AT LEAST ONE FACE WITH 3 OR MORE EDGES
    int min_face_edges     = 5;     // EVERY CONVEX POLYHEDRON MUST HAVE AT LEAST ONE FACE WITH 5 OR FEWER EDGES
    int pvector[max_epf]   = {};    // RECORDS NUMBER OF FACES WITH EACH NUMBER OF EDGES, NO FACE IN FILTER HAS MORE THAN max_epf-1 EDGES
    int origins[2*max_epc] = {};    // NO VORONOI CELL IN FILTER HAS MORE THAN max_epc EDGES
    int origin_c           = 0;
    
    // DETERMINE VERTICES ON FACES WITH MINIMAL EDGES, AND FACES WITH DIFFERENT NUMBERS OF EDGES
    for(int i=0;i<vertex_count;i++)
    {
        for(int j=0;j<vertex_degrees[i];j++)
        {
            int k = ed[i][j];
            if(k >= 0)
            {
                int face[max_epf]={};  // NO SINGLE FACE WILL HAVE MORE THAN max_epf EDGES
                int face_c=0;
                
                ed[i][j]=-1-k;      // INDICATE THAT WE HAVE CHECKED THIS VERTEX
                int l=vcell.cycle_up(ed[i][vertex_degrees[i]+j],k);
                face[face_c++]=k;
                do {
                    int m=ed[k][l];
                    ed[k][l]=-1-m;
                    l=vcell.cycle_up(ed[k][vertex_degrees[k]+l],m);
                    k=m;
                    
                    face[face_c++]=m;
                } while (k!=i);
                
                // KEEP TRACK OF MINIMAL AND MAXIMAL FACE EDGES
                if(face_c>max_face_edges)
                    max_face_edges = face_c;
                if(face_c<min_face_edges)
                {
                    min_face_edges = origin_c = face_c;
                    for(int c=0; c<face_c; c++)
                        origins[c] = face[c];
                }
                else if(face_c==min_face_edges)
                {
                    for(int c=0; c<face_c; c++)
                        origins[origin_c+c] = face[c];
                    origin_c += face_c;
                }
                pvector[face_c]++;
                face_count++;
            }
        }
    }
    
    // RESET EDGES
    for(int i=0;i<vertex_count;i++)
        for(int j=0;j<vertex_degrees[i];j++)
            ed[i][j]=-1-ed[i][j];
    
    // KEEPING TRACK OF THIS WILL ALLOW US TO SPEED UP SOME COMPUTATION, OF BCC
    int likely_bcc=0;
    if(face_count==14 && pvector[4]==6 && pvector[6]==8) likely_bcc=1;   // THIS PVECTOR (0,6,0,8,0,...) OF A SIMPLE POLYHEDRON APPEARS IN 3 DIFFERENT TYPES, WITH SYMMETRIES 4, 8, AND 48
    
    
    ////////////////////////////////////////////////////////////////
    // BUILD THE CANONICAL CODE
    ////////////////////////////////////////////////////////////////
    
    int vertices_temp_labels[max_vpc] = {};         // TEMPORARY LABELS FOR ALL VERTICES; MAX max_vpc VERTICES
    
    int finished   =  0;
    int chirality  = -1;
    int symmetry_counter = 0;     // TRACKS NUMBER OF REPEATS OF A CODE, I.E. SYMMETRY ORDER
    
    for(int orientation=0; orientation<2 && finished==0; orientation++)
    {
        for(int q=0; q<origin_c && finished==0; q++)
        {
            // CLEAR ALL LABELS; MARK ALL BRANCHES OF ALL VERTICES AS NEW
            std::fill(vertices_temp_labels, vertices_temp_labels+vertex_count, 0);
            
            for(int i=0;i<vertex_count;i++)
                for(int j=0;j<vertex_degrees[i];j++)
                    if(ed[i][j]<0) ed[i][j]=-1-ed[i][j];
            
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
            ed[initial][branch] = -1-next;
            
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
                    ed[initial][branch] = -1-next;
                }
                
                else    // NEXT VERTEX *HAS* BEEN VISITED BEFORE
                {
                    int next_branch = ed[initial][vertex_degrees[initial]+branch];
                    int branches_tested = 0;
                    
                    while(ed[next][next_branch] < 0 && branches_tested<vertex_degrees[next])
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
                        ed[initial][branch] = -1-next;
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
    return chirality;
}


int classify_particles_by_voronoi_topology(container_3d& con, Filter &filter)
{
#pragma omp parallel for num_threads(threads)
    for(container_3d::iterator cli=con.begin();cli<con.end();cli++)
    {
        voronoicell_3d vcell;
        if (con.compute_cell(vcell,cli))
        {
            int ijk=cli->ijk,q=cli->q;
            int pid = con.id[ijk][q];
            
            std::vector<int> canonical_code;                // CANONICAL CODE WILL BE STORED HERE
            compute_canonical_code(canonical_code, vcell);
            
            vt_structure_types[pid] = filter.vt_structure_type(canonical_code);
        }
    }
    
    return 0;
}


// FIX: WE DON'T CARE ABOUT CHIRALITY, SYMMETRY, ETC, SO NEED TO
// MAKE AND SORT ONLY ONE LIST OF NUMBERS INSTEAD OF TWO
void classify_particles_by_voronoi_topology(Filter &filter)
{
    double xdim = hi_bound[0]-origin[0];
    double ydim = hi_bound[1]-origin[1];
    
#pragma omp parallel for num_threads(threads)
    for (int pid=0; pid<number_of_particles; pid++)
    {
        // THIS IS OUR LIST OF NEIGHBORS, IT IS A VECTOR.
        int p2pvector[50];
        for(int d=0; d<50; d++) p2pvector[d]=0;
        int maxp = 0;                               // STORE THE MAX NUMBER OF EDGES AMONG NEIGHBORS
        
        unsigned int nneighbors = neighbors_list_char[pid].size();
        std::vector< std::pair <double,int> > prepvect;
        prepvect.reserve(nneighbors);
        
        for(unsigned int q=0; q<nneighbors; q++)
        {
            p2pvector[cell_neighbor_count[neighbors_list_char[pid][q]]]++;
            if(cell_neighbor_count[neighbors_list_char[pid][q]]>maxp) maxp = cell_neighbor_count[neighbors_list_char[pid][q]];
            
            double dx = particle_coordinates[2*neighbors_list_char[pid][q]]   - particle_coordinates[2*pid];
            double dy = particle_coordinates[2*neighbors_list_char[pid][q]+1] - particle_coordinates[2*pid+1];
            if(dx >  xdim/2) dx -= xdim;
            if(dx < -xdim/2) dx += xdim;
            if(dy >  ydim/2) dy -= ydim;
            if(dy < -ydim/2) dy += ydim;
            double theta1 = atan(dy/dx) + 3.14159265359/2.;
            if(dx<0)            theta1 += 3.14159265359;
            
            int p_id = cell_neighbor_count[neighbors_list_char[pid][q]];
            prepvect.push_back(std::make_pair(theta1, p_id));
        }
        std::sort(prepvect.begin(), prepvect.end());
        
        std::vector< std::vector <int> > numbers1;
        std::vector< std::vector <int> > numbers2;
        
        for(unsigned int c=0; c<nneighbors; c++)
        {
            std::vector<int> v;
            for(unsigned int d=0; d<nneighbors; d++)
                v.push_back(prepvect[(c+d)%nneighbors].second);
            numbers1.push_back(v);
        }
        for(unsigned int c=0; c<nneighbors; c++)
        {
            std::vector<int> v;
            for(unsigned int d=0; d<nneighbors; d++)
                v.push_back(prepvect[(c-d+nneighbors)%nneighbors].second);
            numbers2.push_back(v);
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
            canonical_code.insert(canonical_code.begin(), nneighbors);
        }
        else
        {
            canonical_code=numbers1[0];
            canonical_code.insert(canonical_code.begin(), nneighbors);
        }
        vt_structure_types[pid] = filter.vt_structure_type(canonical_code);
    }
    
    return;
}

