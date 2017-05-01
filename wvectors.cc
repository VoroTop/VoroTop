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

////   File: wvectors.cc


#include <vector>
#include <cstring>
#include <fstream>
#include <iostream>

#include "voro++.hh"
#include "filters.hh"
#include "vorotop.hh"

using namespace voro;



std::vector<int> calc_wvector(voronoicell_base &vcell)
{
    return calc_wvector(vcell, 0);
}



////////////////////////////////////////////////////
////
////   FAST CALCULATION OF CANONICAL W-VECTOR.
////   NEXT GENERATION, TAKES AND MODIFIES A REFERENCE
////
////////////////////////////////////////////////////

std::vector<int> calc_wvector(voronoicell_base &vcell, bool extended)
{
    int   edge_count     = vcell.number_of_edges();                    
    int   vertex_count   = vcell.p;    // TOTAL NUMBER OF VERTICES
    int*  vertex_degrees = vcell.nu;   // VERTEX DEGREE ARRAY
    int** ed             = vcell.ed;   // EDGE CONNECTIONS ARRAY
    
    bool stable = 1;
    int face_count     = 0;
    int max_face_edges = 3;         // EVERY CONVEX POLYHEDRON MUST HAVE AT LEAST ONE FACE WITH 3 OR MORE EDGES
    int min_face_edges = 5;         // EVERY CONVEX POLYHEDRON MUST HAVE AT LEAST ONE FACE WITH 5 OR FEWER EDGES
    int pvector[1024]  = {};        // RECORDS NUMBER OF FACES WITH EACH NUMBER OF EDGES, NO FACE HAS MORE THAN 1024 EDGES
    int origins[4096]  = {};        // NO VORONOI CELL WILL HAVE MORE THAN 2048 EDGES
    int origin_c=0;
    
    // DETERMINE VERTICES ON FACES WITH MINIMAL EDGES, AND FACES WITH DIFFERENT NUMBERS OF EDGES
    for(int i=0;i<vertex_count;i++)
    {
        if(vertex_degrees[i]>3) stable = 0;
        
        for(int j=0;j<vertex_degrees[i];j++)
        {
            int k = ed[i][j];
            if(k >= 0)
            {
                int face[1024]={};  // NO SINGLE FACE WILL HAVE MORE THAN 1024 EDGES
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
    
    std::vector<int> canonical_code(2*edge_count);             // CANONICAL CODE WILL BE STORED HERE
    std::vector<int> vertices_temp_labels(vertex_count,0);     // TEMPORARY LABELS FOR ALL VERTICES

    int finished   =  0;
    int chirality  = -1;
    int symmetry_counter = 0;     // TRACKS NUMBER OF REPEATS OF A CODE, I.E. SYMMETRY ORDER
    
    for(int orientation=0; orientation<2 && finished==0; orientation++)
    {
        for(int q=0; q<origin_c && finished==0; q++)
        {
            // CLEAR ALL LABELS; MARK ALL BRANCHES OF ALL VERTICES AS NEW
            std::fill(vertices_temp_labels.begin(), vertices_temp_labels.end(), 0);
            for(int i=0;i<vertex_count;i++)
                for(int j=0;j<vertex_degrees[i];j++)
                    if(ed[i][j]<0) ed[i][j]=-1-ed[i][j];
            
            int initial = origins[q];
            int next;
            int branch;
            
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
    
    if(extended==0)
        canonical_code.push_back(1);
    
    else // extended==1
    {
        for(unsigned int i=3; i<max_face_edges+1; i++)
            canonical_code.push_back(pvector[i]);
        canonical_code.push_back(face_count);
        canonical_code.push_back(symmetry_counter);
        canonical_code.push_back(chirality);
        canonical_code.push_back(stable);
        canonical_code.push_back(max_face_edges+1);
        canonical_code.push_back(2*edge_count);
    }
    
    return canonical_code;
}





