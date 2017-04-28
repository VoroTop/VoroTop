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
    unsigned int vertex_count   = vcell.p;  // TOTAL NUMBER OF VERTICES
    int*         vertex_degrees = vcell.nu; // VERTEX DEGREE ARRAY
    int** ed                    = vcell.ed; // EDGE CONNECTIONS ARRAY
    
    unsigned int face_count = 0;
    std::vector<int> pvector(5,0);  // RECORDS NUMBER OF FACES WITH EACH NUMBER OF EDGES
    int min_face_edges = 5;         // EVERY CONVEX POLYHEDRON MUST HAVE AT LEAST ONE FACE WITH 5 OR FEWER EDGES
    bool stable = 1;
    std::vector<int> origins;
    
    
    // DETERMINE NUMBER OF FACES WITH DIFFERENT NUMBERS OF EDGES
    for(int i=0;i<vertex_count;i++)
    {
        if(vertex_degrees[i]>3) stable = 0;
        
        for(int j=0;j<vertex_degrees[i];j++)
        {
            int k = ed[i][j];
            if(k >= 0)  // ENSURES WE HAVE NOT YET CHECKED OUT THIS EDGE
            {
                std::vector<int> face;
                
                ed[i][j]=-1-k;      // RECORD THAT WE HAVE CHECKED OUT THIS VERTEX
                int l=vcell.cycle_up(ed[i][vertex_degrees[i]+j],k);
                face.push_back(k);
                do {
                    int m=ed[k][l];
                    ed[k][l]=-1-m;
                    l=vcell.cycle_up(ed[k][vertex_degrees[k]+l],m);
                    k=m;
                    
                    face.push_back(m);
                } while (k!=i);
                
                // KEEPING TRACK OF ALL MINIMAL FACE EDGES
                if(face.size()<min_face_edges)
                {
                    origins = face;
                    min_face_edges = face.size();
                }
                else if(face.size()==min_face_edges)
                    origins.insert(origins.end(), face.begin(), face.end());
                
                if(face.size()<pvector.size()) pvector[face.size()]++;
                else
                {
                    pvector.resize(face.size()+1,0);
                    pvector[face.size()]++;
                }
                
                face_count++;
            }
        }
    }
    
    // RESET ALL EDGES THAT WERE CHANGED WHILE DETERMINING P-VECTOR
    for(int i=0;i<vertex_count;i++)
        for(int j=0;j<vertex_degrees[i];j++)
            ed[i][j]=-1-ed[i][j];
    
    // KEEPING TRACK OF THIS WILL ALLOW US TO SPEED UP SOME COMPUTATION, OF BCC
    int likely_bcc=0;
    if(face_count==14 && pvector[4]==6 && pvector[6]==8) likely_bcc=1;   // THIS PVECTOR (0,6,0,8,0,...) OF A SIMPLE POLYHEDRON APPEARS IN 3 DIFFERENT TYPES, WITH SYMMETRIES 4, 8, AND 48
    
    
    ////////////////////////////////////////////////////////////////
    // BUILD THE CANONICAL CODE
    ////////////////////////////////////////////////////////////////
    
    int edge_count = face_count + vertex_count - 2;     // USE EULER'S FORMULA TO COMPUTE NUMBER OF EDGES
    std::vector<int> wvector(2*edge_count);
    
    int vertices_temp_labels[1024];      // TEMPORARY LABELS FOR ALL VERTICES STORED HERE
    int vertex_visited      [1024][16];  // TRACKS WHICH VERTEX/EDGE PAIRS WE HAVE VISITED
    
    
    // IF E=edge_count IS THE NUMBER OF EDGES IN THE POLYHEDRON, THEN WE WILL
    // CONSTRUCT 4E CODES, EACH OF LENGTH 2E.
    int finished   =  0;
    int chirality  = -1;
    int symmetry_counter=0;     // TRACKS NUMBER OF REPEATS OF A CODE, I.E. SYMMETRY ORDER
    
    for(int orientation=0; orientation<2 && finished==0; orientation++)
    {
//         std::cout << orientation << '\n';
        for(int q=0; q<vertex_count; q++)
        {
//            std::cout << "Vertex " << q << '\n';
            for(int r=0; r<vertex_degrees[q]; r++)
            {
                // CLEAR ALL LABELS; MARK ALL BRANCHES OF ALL VERTICES AS NEW (0)
                //std::cout << "About to clear all " << vertex_count << " temp_lables\n";
                for(int k=0; k<vertex_count; k++)
                    vertices_temp_labels[k] = 0;
                for(int i=0;i<vertex_count;i++)
                    for(int j=0;j<vertex_degrees[i];j++)
                        if(ed[i][j]<0) ed[i][j]=-1-ed[i][j];

//                std::cout << "Vertex " << q << ", branch " << r << '\n';
                int initial = q;
                int branch  = r;
                int next    = ed[initial][branch];
                ed[initial][branch] = -1-next;

                int current_code_length   = 0;
                int current_highest_label = 1;
                int continue_code         = 0;     // 0: UNDECIDED; 1: GO AHEAD, DO NOT EVEN CHECK.
                if(q==0 && orientation==0)         // FIRST CODE, GO AHEAD
                    continue_code=1;
                

                
                vertices_temp_labels[initial] = current_highest_label++;
                wvector[current_code_length]  = vertices_temp_labels[initial];
                current_code_length++;
                
                int end_flag=0;
                while(end_flag==0)
                {
//                    std::cout << "inside, and next is " << next << "\n";
///////
                    // NEXT VERTEX HAS NOT BEEN VISITED; TAKE RIGHT-MOST BRANCH TO CONTINUE.
                    if(vertices_temp_labels[next]==0)
                    {
//                        std::cout << "Not visited\n";
                        //   LABEL THE NEW VERTEX
                        vertices_temp_labels[next] = current_highest_label++;
                        
                        if(continue_code==0)
                        {
                            if(vertices_temp_labels[next]>wvector[current_code_length]) break;
                            if(vertices_temp_labels[next]<wvector[current_code_length])
                            {
                                symmetry_counter = 0;
                                continue_code    = 1;
                                if(orientation==1) chirality=1;
                            }
                        }
                        
                        //   BUILD THE CODE
                        wvector[current_code_length] = vertices_temp_labels[next];
                        current_code_length++;
                        
                        //   FIND THE NEXT DIRECTION TO MOVE ALONG
                        //   UPDATE INITIAL AND BRANCH TO RELOOP
                        branch  = vcell.cycle_up(ed[initial][vertex_degrees[initial]+branch],next);
                        initial = next;
                        next    = ed[initial][branch];
                        ed[initial][branch] = -1-next;
                    }
                    
                    else    // NEXT VERTEX *HAS* BEEN VISITED BEFORE
                    {
//                        std::cout << "Visited " << vertices_temp_labels[next] << "\n";

                        int next_branch = ed[initial][vertex_degrees[initial]+branch];  // BEGIN ON RETURN BRANCH
                        int branches_tested = 0;
                        
                        while(ed[next][next_branch] < 0 && branches_tested<vertex_degrees[next])
                        {
                            next_branch = vcell.cycle_up(next_branch,next);
                            branches_tested++;
                        }

                        if(branches_tested < vertex_degrees[next])
                        {
                            if(continue_code==0)
                            {
                                if(vertices_temp_labels[next]>wvector[current_code_length]) break;
                                if(vertices_temp_labels[next]<wvector[current_code_length])
                                {
                                    symmetry_counter=0;
                                    continue_code=1;
                                    if(orientation==1) chirality=1;
                                }
                            }
                            
                            // BUILD THE CODE
                            wvector[current_code_length] = vertices_temp_labels[next];
                            current_code_length++;
                            
                            // FIND NEXT BRANCH, EASY IN THIS CASE
                            initial = next;
                            branch  = next_branch;
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
///////
                }
            }
        }
        
        // AFTER MAKING ALL CODES FOR ONE ORIENTATION, FLIP ORIENTATION OF EDGES AT EACH
        // VERTEX, AND REPEAT THE ABOVE FOR THE OPPOSITE ORIENTATION.
        
        if(orientation==3)
        {
            for(int i=0;i<vertex_count;i++)
            {
                for(int j=0; j<vertex_degrees[i]; j++)
                    ed[i][j+vertex_degrees[i]] = vertex_degrees[ed[i][j]]-ed[i][j+vertex_degrees[i]]-1;
                
                for(int j=vertex_degrees[i]/2; j--;)
                {
                    int temp = ed[i][j];
                    ed[i][j] = ed[i][vertex_degrees[i]-j-1];
                    ed[i][vertex_degrees[i]-j-1] = temp;
                    
                    temp = ed[i][j+vertex_degrees[i]];
                    ed[i][j+vertex_degrees[i]] = ed[i][vertex_degrees[i]-j-1+vertex_degrees[i]];
                    ed[i][vertex_degrees[i]-j-1+vertex_degrees[i]] = temp;
                }
            }
            
            for(int i=0;i<vertex_count;i++)
                for(int j=0; j<vertex_degrees[i]; j++)
                    ed[ed[i][j]][ed[i][vertex_degrees[i]+j]]=i;
        }
    }
    
//    for(int c=0; c<wvector.size(); c++)
//        std::cout << wvector[c] << ",";
//    std::cout << '\n';
    
    if(extended==0)
        wvector.push_back(1);
    
    else // extended==1
    {
        for(unsigned int i=3; i<pvector.size(); i++)
            wvector.push_back(pvector[i]);
        wvector.push_back(face_count);
        wvector.push_back(symmetry_counter);
        wvector.push_back(chirality);
        wvector.push_back(stable);
        wvector.push_back(pvector.size());
        wvector.push_back(2*edge_count);
    }
    
    return wvector;
}





////////////////////////////////////////////////////
////
////   FAST CALCULATION OF CANONICAL W-VECTOR.
////   NEXT GENERATION, TAKES AND MODIFIES A REFERENCE
////
////////////////////////////////////////////////////

std::vector<int> calc_wvectorOLD(voronoicell_base &vcell, bool extended)
{
    unsigned int vertex_count   = vcell.p;  // TOTAL NUMBER OF VERTICES
    int*         vertex_degrees = vcell.nu; // VERTEX DEGREE ARRAY
    int** ed                    = vcell.ed; // EDGE CONNECTIONS ARRAY
    
    unsigned int face_count = 0;
    std::vector<int> pvector(5,0);  // RECORDS NUMBER OF FACES WITH EACH NUMBER OF EDGES
    int min_face_edges = 5;         // EVERY CONVEX POLYHEDRON MUST HAVE AT LEAST ONE FACE WITH 5 OR FEWER EDGES
    bool stable = 1;
    
    // FOR EACH VERTEX, RECORD SMALLEST NUMBER OF EDGES OF ALL ADJACENT FACES.
    // THIS WILL REDUCE NUMBER OF COMPUTATIONS NECESSARY TO DO LATER.  WE BEGIN
    // BY SETTING THIS TO 6, SINCE EVERY VORONOI CELL MUST HAVE AT LEAST ONE
    // FACE WITH 5 OR FEWER EDGES, WE WON'T CHECK ANYTHING WITH MORE THAN 5.
    std::vector<int> lowest_adjacent_face_degree(vertex_count,6);
    
    // DETERMINE NUMBER OF FACES WITH DIFFERENT NUMBERS OF EDGES
    for(int i=0;i<vertex_count;i++)
    {
        if(vertex_degrees[i]>3) stable = 0;
        
        for(int j=0;j<vertex_degrees[i];j++)
        {
            int k = ed[i][j];
            if(k >= 0)  // ENSURES WE HAVE NOT YET CHECKED OUT THIS EDGE
            {
                unsigned int side_count = 1;
                std::vector<int> temp;
                
                ed[i][j]=-1-k;      // RECORD THAT WE HAVE CHECKED OUT THIS VERTEX
                int l=vcell.cycle_up(ed[i][vertex_degrees[i]+j],k);
                temp.push_back(k);
                do {
                    int m=ed[k][l];
                    ed[k][l]=-1-m;
                    l=vcell.cycle_up(ed[k][vertex_degrees[k]+l],m);
                    k=m;
                    
                    side_count++;
                    temp.push_back(m);
                } while (k!=i);
                
                for(int c=0; c<temp.size(); c++)
                    if(side_count<lowest_adjacent_face_degree[temp[c]]) lowest_adjacent_face_degree[temp[c]] = side_count;
                
                if(side_count<pvector.size()) pvector[side_count]++;
                else
                {
                    pvector.resize(side_count+1,0);
                    pvector[side_count]++;
                }
                
                if(side_count<min_face_edges)
                    min_face_edges = side_count;
                
                face_count++;
            }
        }
    }
    
    // RESET ALL EDGES THAT WERE CHANGED WHILE DETERMINING P-VECTOR
    for(int i=0;i<vertex_count;i++)
        for(int j=0;j<vertex_degrees[i];j++)
            ed[i][j]=-1-ed[i][j];
    
    // KEEPING TRACK OF THIS WILL ALLOW US TO SPEED UP SOME COMPUTATION, OF BCC
    int likely_bcc=0;
    if(face_count==14 && pvector[4]==6 && pvector[6]==8) likely_bcc=1;   // THIS PVECTOR (0,6,0,8,0,...) OF A SIMPLE POLYHEDRON APPEARS IN 3 DIFFERENT TYPES, WITH SYMMETRIES 4, 8, AND 48
    
    
    ////////////////////////////////////////////////////////////////
    // BUILD THE CANONICAL CODE
    ////////////////////////////////////////////////////////////////
    
    int current_code_length;                        // FOR KEEPING TRACK OF WHERE WE ARE IN BUILDING A CODE
    int edge_count = face_count + vertex_count - 2;            // USE EULER'S FORMULA TO COMPUTE NUMBER OF EDGES
    std::vector<int> wvector(2*edge_count);
    
    int vertices_temp_labels[1024];      // TEMPORARY LABELS FOR ALL VERTICES STORED HERE
    int vertex_visited      [1024][64];  // TRACKS WHICH VERTEX/EDGE PAIRS WE HAVE VISITED
    
    
    // IF E=edge_count IS THE NUMBER OF EDGES IN THE POLYHEDRON, THEN WE WILL
    // CONSTRUCT 4E CODES, EACH OF LENGTH 2E.
    int first_code =  1;        // KEEPS TRACK OF WHETHER THIS IS FIRST CODE RECORDED
    int finished   =  0;
    int chirality  = -1;
    int symmetry_counter=0;     // TRACKS NUMBER OF REPEATS OF A CODE, I.E. SYMMETRY ORDER
    
    for(int orientation=0; orientation<2 && finished==0; orientation++)
    {
        // FOR EACH ORIENTATION, WE CONSTRUCT THREE CODES FOR EACH
        // STARTING VERTEX, ONE FOR EACH EDGE LEAVING IT.
        for(int i=0; i<vertex_count; i++)
        {
            // THIS VERTEX IS ADJACENT TO MANY FACES WITH MANY EDGES,
            // SO IT CAN'T GENERATE A LEXICOGRAPHICALLY SMALLEST CODE
            if(lowest_adjacent_face_degree[i] > min_face_edges) continue;
            
            // BEGIN CODE BY CHOOSING ONE DIRECTION ALONG WHICH TO TRAVEL.
            for(int j=0; j<vertex_degrees[i] && finished==0; j++)
            {
                // REQUIRE THAT NEXT VERTEX BE ADJACENT TO A SMALL FACE
                if(lowest_adjacent_face_degree[ed[i][j]] > min_face_edges) continue;
                
                current_code_length=0;
                int current_highest_label=1;
                int continue_code = 0;    // 0: UNDECIDED; 1: GO AHEAD, DO NOT EVEN CHECK.
                if(first_code==1) continue_code=1;
                
                // CLEAR ALL LABELS; MARK ALL BRANCHES OF ALL VERTICES AS NEW
                for(int k=0; k<vertex_count; k++)
                {
                    vertices_temp_labels[k] = 0;
                    for(int l=0; l<vertex_degrees[k]; l++)
                        vertex_visited[k][l] = 0;
                }
                
                vertices_temp_labels[i] = current_highest_label++;
                wvector[current_code_length]=vertices_temp_labels[i];
                current_code_length++;
                
                int initial = i;   // VERTEX WE ARE LEAVING, i IS IN {0,1,2,...,p-1}
                int branch  = j;   // EDGE WE ARE LEAVING ALONG; INDEX RELATIVE TO I IN ed[][]; j IS IN {0,1,...,nu[i]}
                vertex_visited[initial][branch]=1;
                
                
                // THIS SECTION BUILDS EACH CODE, FOLLOWING THE WEINBERG RULES FOR TRAVERSING A
                // GRAPH MAKING A HAMILTONIAN PATH, LABELING THE VERTICES ALONG THE WAY, AND
                // RECORDING THE VERTICES VISITED.
                int end_flag=0;
                while(end_flag==0)
                {
                    // NEXT VERTEX HAS NOT BEEN VISITED BEFORE; TAKE RIGHT-MOST BRANCH TO CONTINUE.
                    if(vertices_temp_labels[ed[initial][branch]]==0)
                    {
                        //   LABEL THE NEW VERTEX
                        vertices_temp_labels[ed[initial][branch]] = current_highest_label++;
                        
                        if(continue_code==0)
                        {
                            if(vertices_temp_labels[ed[initial][branch]]>wvector[current_code_length]) break;
                            if(vertices_temp_labels[ed[initial][branch]]<wvector[current_code_length])
                            { symmetry_counter=0; continue_code=1; if(orientation==1) chirality=1; }
                        }
                        
                        //   BUILD THE CODE
                        wvector[current_code_length]=vertices_temp_labels[ed[initial][branch]];
                        current_code_length++;
                        
                        //   FIND THE NEXT DIRECTION TO MOVE ALONG
                        //   UPDATE INITIAL AND BRANCH TO RELOOP
                        int old_branch = branch;
                        branch = ed[initial][branch+vertex_degrees[initial]]+1 == vertex_degrees[ed[initial][branch]]?0:ed[initial][branch+vertex_degrees[initial]]+1;
                        initial = ed[initial][old_branch];
                        vertex_visited[initial][branch]=1;
                    }
                    
                    else    // NEXT VERTEX *HAS* BEEN VISITED BEFORE
                    {
                        // IF RETURN BRANCH HAS NOT BEEN TAKEN
                        if(vertex_visited[ed[initial][branch]][ed[initial][branch+vertex_degrees[initial]]]==0)
                        {
                            if(continue_code==0)
                            {
                                if(vertices_temp_labels[ed[initial][branch]]>wvector[current_code_length]) break;
                                if(vertices_temp_labels[ed[initial][branch]]<wvector[current_code_length])
                                { symmetry_counter=0; continue_code=1; if(orientation==1) chirality=1; }
                            }
                            
                            // BUILD THE CODE
                            wvector[current_code_length]=vertices_temp_labels[ed[initial][branch]];
                            current_code_length++;
                            
                            // FIND NEXT BRANCH, EASY IN THIS CASE
                            int old_branch = branch;
                            branch  = ed[initial][branch+vertex_degrees[initial]];
                            initial = ed[initial][old_branch];
                            vertex_visited[initial][branch]=1;
                        }
                        
                        // IF RETURN BRANCH *HAS* BEEN TAKEN ALREADY
                        else
                        {
                            // FINDS NEXT OPEN BRANCH IF EXISTS
                            int index = ed[initial][branch+vertex_degrees[initial]]+1 == vertex_degrees[ed[initial][branch]]?0:ed[initial][branch+vertex_degrees[initial]]+1;
                            while(index != ed[initial][branch+vertex_degrees[initial]] && vertex_visited[ed[initial][branch]][index]!=0)
                                index = index+1 == vertex_degrees[ed[initial][branch]]?0:index+1;
                            
                            // IF SUCH A BRANCH EXISTS
                            if(vertex_visited[ed[initial][branch]][index]==0)
                            {
                                // BUILD THE CODE
                                if(continue_code==0)
                                {
                                    if(vertices_temp_labels[ed[initial][branch]]>wvector[current_code_length]) break;
                                    if(vertices_temp_labels[ed[initial][branch]]<wvector[current_code_length])
                                    { symmetry_counter=0; continue_code=1; if(orientation==1) chirality=1; }
                                }
                                
                                if(continue_code==1)    // I THINK WE CAN REMOVE THIS CONDITIONAL
                                    wvector[current_code_length]=vertices_temp_labels[ed[initial][branch]];
                                current_code_length++;
                                
                                int old_branch = branch;
                                branch  = index;
                                initial = ed[initial][old_branch];
                                vertex_visited[initial][branch]=1;
                            }
                            
                            // NO OPEN BRANCH EXISTS; CODE COMPLETED.
                            else
                            {
                                // I THINK WE'RE DONE
                                end_flag=1;
                                first_code=0;   // WE HAVE NOW GOTTEN AT LEAST ONE CODE
                                
                                if(likely_bcc && symmetry_counter>4 && orientation==0) { chirality=0; symmetry_counter = 48; finished=1; }
                                else if(chirality==-1 && orientation==1)               { chirality=0; symmetry_counter *= 2; finished=1; }
                                else    symmetry_counter++;
                            }
                        }
                    }
                }
            }
        }
        
        // AFTER MAKING ALL CODES FOR ONE ORIENTATION, FLIP ORIENTATION OF EDGES AT EACH
        // VERTEX, AND REPEAT THE ABOVE FOR THE OPPOSITE ORIENTATION.
        if(orientation==0)
        {
            for(int i=0;i<vertex_count;i++)
            {
                for(int j=0; j<vertex_degrees[i]; j++)
                    ed[i][j+vertex_degrees[i]] = vertex_degrees[ed[i][j]]-ed[i][j+vertex_degrees[i]]-1;
                
                for(int j=vertex_degrees[i]/2; j--;)
                {
                    int temp = ed[i][j];
                    ed[i][j] = ed[i][vertex_degrees[i]-j-1];
                    ed[i][vertex_degrees[i]-j-1] = temp;
                    
                    temp = ed[i][j+vertex_degrees[i]];
                    ed[i][j+vertex_degrees[i]] = ed[i][vertex_degrees[i]-j-1+vertex_degrees[i]];
                    ed[i][vertex_degrees[i]-j-1+vertex_degrees[i]] = temp;
                }
            }
        }
    }
    
    if(extended==0)
        wvector.push_back(1);
    
    else // extended==1
    {
        for(unsigned int i=3; i<pvector.size(); i++)
            wvector.push_back(pvector[i]);
        wvector.push_back(face_count);
        wvector.push_back(symmetry_counter);
        wvector.push_back(chirality);
        wvector.push_back(stable);
        wvector.push_back(pvector.size());
        wvector.push_back(2*edge_count);
    }
    
    return wvector;
}





