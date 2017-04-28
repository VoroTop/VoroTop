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

////   File: output.cc


#include <queue>
#include <cstring>
#include <fstream>
#include <iostream>

#include "mtrand.h"
#include "import.hh"
#include "voro++.hh"
#include "filters.hh"
#include "vorotop.hh"
#include "functions.hh"

using namespace voro;



bool vector_size(const std::vector<int>& a,const std::vector<int>& b)
{
    return (a.size() > b.size());
}



////////////////////////////////////////////////////
////
////   MAKE A CFG FILE USING A FILTER.  IN A FILTER
////   STRUCTURE TYPES ARE NUMBERED 1,2,3,...  IF A
////   TOPOLOGY IS NOT IN THE FILTER, THE STRUCTURE
////   TYPE IS 0; IF THE VORONOI CELL CANNOT BE
////   COMPUTED, THE STRUCTURE TYPE IS -1.
////
////////////////////////////////////////////////////

void create_cfg_file(std::string filename, Filter &filter)
{
    ////////////////////////////////////////////////////
    ////
    ////   OPEN OUTPUT FILE
    ////
    ////////////////////////////////////////////////////
    
    std::string cfg_file_name(filename);
    cfg_file_name.append(".cfg");
    std::ofstream cfg_file;
    cfg_file.open(cfg_file_name.c_str());
    
    
    ////////////////////////////////////////////////////
    ////
    ////    PRINT HEADER INFO TO FILE
    ////
    ////////////////////////////////////////////////////
    
    cfg_file << "Number of particles = " << number_of_particles << "\n";
    cfg_file << "A = " << cfg_lscale << " Angstrom (basic length-scale)\n";
    for(int c=0; c<3; c++)
        for(int d=0; d<3; d++)
            if(supercell_edges[c][d]!=0)
                cfg_file << "H0(" << c+1 << "," << d+1 << ") = " << supercell_edges[c][d] << '\n';
    
    if(file_type==1 || no_velocity==1) cfg_file << ".NO_VELOCITY.\n";
    
    int auxiliary_counter=0;
    
    int entry_count = attribute_labels.size()+3;
    
    if(f_switch) entry_count+=1;  // FOR VORONOI TOPOLOGY
    if(df_switch)entry_count+=2;  // FOR DISTRIBUTION TYPE AND COUNT
    if(c_switch) entry_count+=2;  // FOR CLUSTER ID AND COUNT
    
    cfg_file << "entry_count = " << entry_count << "\n";
    
    
    // ATTRIBUTES FROM INPUT DATA
    for(int d=0; d<attribute_labels.size(); d++)
        cfg_file << "auxiliary[" << auxiliary_counter++ << "] = " << attribute_labels[d] << "\n";
    
    // NEW VORONOI ANALYSIS ATTRIBUTES
    if     ( f_switch && !df_switch)
        cfg_file << "auxiliary[" << auxiliary_counter++ << "] = vt\n";
    
    else if(!f_switch &&  df_switch)
    {
        cfg_file << "auxiliary[" << auxiliary_counter++ << "] = distribution type\n";
        cfg_file << "auxiliary[" << auxiliary_counter++ << "] = distribution count\n";
    }
    
    else if( f_switch &&  df_switch)
    {
        cfg_file << "auxiliary[" << auxiliary_counter++   << "] = vt\n";
        cfg_file << "auxiliary[" << auxiliary_counter++ << "] = distribution type\n";
        cfg_file << "auxiliary[" << auxiliary_counter++ << "] = distribution count\n";
    }
    
    if(c_switch)
    {
        cfg_file << "auxiliary[" << auxiliary_counter++   << "] = cluster index\n";
        cfg_file << "auxiliary[" << auxiliary_counter++   << "] = cluster size\n";
    }
    
    cfg_file << cfg_atomic_mass << "\n";
    cfg_file << cfg_chem_symbol << "\n";
    
    
    ////////////////////////////////////////////////////
    ////
    ////   ITERATE, DETERMINE TYPES, AND OUTPUT DATA
    ////
    ////////////////////////////////////////////////////
    
    filter.sort_by_wvector();
    for(int c=0; c<number_of_particles; c++)
    {
        double x = xcoord[c];
        double y = ycoord[c];
        double z = zcoord[c];
        
        if(file_type==1 && scaled_coordinates==0)
        {
            double a = supercell_edges[0][0];
            double b = supercell_edges[1][1];
            double c = supercell_edges[2][2];
            double d = supercell_edges[1][0];
            double e = supercell_edges[2][1];
            double f = supercell_edges[2][0];
            
            double newx = x/a - y*d/a/b + z*(d*e-b*f)/a/b/c;
            double newy = y/b - z*e/b/c;
            double newz = z/c;
            
            if(newx<0) newx+=1.; if(newx>=1) newx-=1.;
            if(newy<0) newy+=1.; if(newy>=1) newy-=1.;
            if(newz<0) newz+=1.; if(newz>=1) newz-=1.;
            
            x=newx;
            y=newy;
            z=newz;
        }
        
        int index = filter.wvector_index(all_wvectors[c]);
        int dtype = filter.get_entry_type (index);          // DISTRIBUTION TYPE
        int count = filter.get_entry_count(index);          // DISTRIBUTION COUNT
        int vtype = dtype;                                  // VORONOI TOPOLOGY TYPE
        if(vtype > filter.get_max_ff_type())
            vtype = 0;
        
        ////////////////////////////////////////////////////
        ////
        ////   OUTPUT ALL DATA TO FILE
        ////
        ////////////////////////////////////////////////////
        
        // OUTPUT COORDINATES
        cfg_file << x << '\t' << y << '\t' << z << '\t';
        
        // OUTPUT DATA IN INITIAL FILE
        for(int d=0; d<attribute_labels.size(); d++) cfg_file << particle_data[c][d] << '\t';
        
        // OUTPUT TOPOLOGICAL TYPE, DISTRIBUTION DATA, CLUSTER DATA AS SPECIFIED
        if(f_switch) cfg_file << vtype << '\t';
        if(df_switch)cfg_file << dtype << '\t' << count << '\t';
        if(c_switch) cfg_file << cluster_index[c] << '\t';
        if(c_switch) cfg_file << cluster_sizes[c] << '\t';  // REQUIRES FIXING, PLACEHOLDER
        
        cfg_file << '\n';
    }
}



////////////////////////////////////////////////////
////
////   MAKE A CFG FILE USING FILTER, AND ALSO CLUSTERING
////   SETS OF NEIGHBORING PARTICLES; USEFUL FOR ANALYSIS OF
////   DEFECTS, SUCH AS THOSE IN PHASE-TRANSITIONS.
////
////////////////////////////////////////////////////

void   cluster_analysis(Filter &filter)
{
    filter.sort_for_clustering();
    for(int counter=0; counter<number_of_particles; counter++)
        cluster_index  [counter] = filter.wvector_type(all_wvectors[counter]);
    
    int total_defect_particles  = 0;
    int total_crystal_particles = 0;
    
    std::vector <std::vector <int> > clusters_good;
    std::vector <std::vector <int> > clusters_bad;
    
    // BUILD DEFECT CLUSTERS IN O(N) TIME;
    std::vector <int> visited(number_of_particles,0);
    
    for(int d=0; d<number_of_particles; d++)
    {
        if(cluster_index[d]==0) total_defect_particles++;
        else                    visited[d]=1;
    }
    
    for(int d=0; d<number_of_particles; d++)
    {
        if(visited[d]==0)
        {
            std::vector<int> cluster;
            std::queue<int> q;
            q.push(d);
            visited[d]=1;
            
            while(!q.empty())
            {
                int w = q.front();
                cluster.push_back(w);
                q.pop();
                
                for(int e=0; e<neighbors_list_c[w]; e++)
                {
                    if(visited[neighbors_list[w][e]]==0)
                    {
                        visited[neighbors_list[w][e]]=1;
                        q.push(neighbors_list[w][e]);
                    }
                }
            }
            clusters_bad.push_back(cluster);
        }
    }
    sort(clusters_bad.begin(),clusters_bad.end(),vector_size);
    
    // BUILD CRYSTALLINE CLUSTERS IN O(N) TIME;
    std::fill(visited.begin(),visited.end(), 0);
    for(int d=0; d<number_of_particles; d++)
    {
        if(cluster_index[d]>0)  total_crystal_particles++;
        else                    visited[d]=1;
    }
    
    for(int d=0; d<number_of_particles; d++)
    {
        if(visited[d]==0)
        {
            std::vector<int> cluster;
            std::queue<int> q;
            q.push(d);
            visited[d]=1;
            
            while(!q.empty())
            {
                int w = q.front();
                cluster.push_back(w);
                q.pop();
                
                for(int e=0; e<neighbors_list_c[w]; e++)
                {
                    if(visited[neighbors_list[w][e]]==0)
                    {
                        visited[neighbors_list[w][e]]=1;
                        q.push(neighbors_list[w][e]);
                    }
                }
            }
            clusters_good.push_back(cluster);
        }
    }
    sort(clusters_good.begin(),clusters_good.end(),vector_size);
    

    // ASSIGN EACH PARTICLE TWO NUMBERS
    // 1. CLUSTER INDEX.  -1,-2,-3,... INDICATE DEFECT CLUSTERS
    //    IN DESCENDING ORDER OF SIZE; 1,2,3,... INDICATE CRYSTAL
    //    CLUSTERS IN DESCENDING ORDER OF SIZE.
    // 2. NUMBER OF PARTICLES IN CLUSTER CONTAINING THIS PARTICLE
    for(int c=0; c<clusters_bad.size(); c++)
    {
        int size = clusters_bad[c].size();
        for(int d=0; d<size; d++)
        {
            cluster_index[clusters_bad[c][d]]=-c-1;
            cluster_sizes[clusters_bad[c][d]]=size;
        }
    }
    
    for(int c=0; c<clusters_good.size(); c++)
    {
        int size = clusters_good[c].size();
        for(int d=0; d<size; d++)
        {
            cluster_index[clusters_good[c][d]]=c+1;
            cluster_sizes[clusters_good[c][d]]=size;
        }
    }
    
    
    ////////////////////////////////////////////////////
    ////
    ////   OUTPUT BASIC ANALYSIS OF CLUSTERING DATA TO SCREEN
    ////
    ////////////////////////////////////////////////////
    
    int max_cluster_size = 100;
    std::vector<int> c_sizes(max_cluster_size,0);
    double sum_of_squares = 0;

    for(int c=0; c<clusters_bad.size(); c++)
    {
        if(clusters_bad[c].size()<max_cluster_size)
            c_sizes[clusters_bad[c].size()]++;
        sum_of_squares += clusters_bad[c].size()*clusters_bad[c].size();
    }
    
    std::cout << timestep               << '\t';
    std::cout << number_of_particles    << '\t';
    std::cout << total_defect_particles << '\t';
    std::cout << clusters_bad.size()    << '\t';
    std::cout << clusters_good.size()   << '\t';
    
    if(!clusters_bad.empty())
    {
        double average = double(total_defect_particles)/double(clusters_bad.size());
        double stdev   = sqrt(double(sum_of_squares)/double(clusters_bad.size()) - average*average);
        
        std::cout << (clusters_bad. back()).size() << '\t';  // SIZE OF SMALLEST DEFECT CLUSTER
        std::cout << (clusters_bad.front()).size() << '\t';  // SIZE OF LARGEST DEFECT CLUSTER
        std::cout << average                       << '\t';  // AVERAGE SIZE OF CLUSTER
        std::cout << stdev                         << '\t';  // STANDARD DEVIATION OF CLUSTER SIZE
    }
    else
        std::cout << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t';
    
    for(int c=1; c<max_cluster_size; c++)                    // NUMBER OF CLUSTERS WITH 1, 2, 3, ..., max_cluster_size-1 PARTICLES
        std::cout << c_sizes[c] << '\t';

    std::cout << '\n';
}



////////////////////////////////////////////////////
////
////   COMPUTES AND PRINTS DISTRIBUTION OF TYPES IN
////   SYSTEM AFTER RANDOM GAUSSIAN PERTURBATIONS.
////   COMPUTES perturbation_samples RANDOM
////   PERTURBATIONS OF SYSTEM, EACH OF STANDARD
////   DEVIATION perturbation_size, AND SUMS THE
////   DISTRIBUTION OVER ALL SYSTEMS.  THIS CAN
////   BE USEFUL IN COMPUTING TYPES ASSOCIATED TO
////   REALISTIC VERSIONS OF IDEAL SYSTEMS.
////
////////////////////////////////////////////////////

void calc_gaussian_distribution(container_periodic& con, particle_order& vo, Filter &filter)
{
    MTRand_open mtrand; // random float in (0, 1), using Mersenne Twister
    
    for(int loop = 0; loop < perturbation_samples; loop++)
    {
        // NEW, PERTURBED VERSION OF THE PRIMARY SYSTEM
        particle_order voP;
        container_periodic conP (supercell_edges[0][0],supercell_edges[1][0],supercell_edges[1][1],
                                 supercell_edges[2][0],supercell_edges[2][1],supercell_edges[2][2],
                                 n_x,n_y,n_z,8);
        
        int pid;
        double x,y,z;
        c_loop_order_periodic vlo(con,vo);
        
        // ADD PERTURBATION OF PARTICLES IN PRIMARY CONTAINER TO PERTURBED CONTAINER
        if(vlo.start()) do
        {
            pid = vlo.id[vlo.ijk][vlo.q];
            vlo.pos(x,y,z);
            
            double rand1 = mtrand();
            double rand2 = mtrand();
            double rand3 = mtrand();
            double rand4 = mtrand();
            
            // GENERATE GAUSSIAN NOISE USING Box–Muller TRANSFORM
            double randx = sqrt(-2.*log(rand1))*cos(2.*M_PI*rand2)*perturbation_size;
            double randy = sqrt(-2.*log(rand1))*sin(2.*M_PI*rand2)*perturbation_size;
            double randz = sqrt(-2.*log(rand3))*cos(2.*M_PI*rand4)*perturbation_size;
            
            x += randx;
            y += randy;
            z += randz;
            
            // ADD ``PERTURBED'' COPY
            conP.put(voP,pid,x,y,z);
            
            // WE SHOULD FIX THIS FOR PERIODIC BOUNDARY CONDITIONS; DIFFICULT FOR
            // NON-ORTHORHOMBIC CELLS
            
        } while(vlo.inc());
        
        
        // COMPUTE DISTRIBUTION OF TYPES IN THE PERTURBED VERSION
        voronoicell c;
        c_loop_order_periodic vloP(conP,voP);
        if(vloP.start()) do if(conP.compute_cell(c,vloP))
            filter.increment_or_add(calc_wvector(c),1);
        while(vloP.inc());
    }
    filter.sort_by_count();
}



////////////////////////////////////////////////////
////
////   PRINT WEINBERG VECTOR FOR EACH PARTICLE
////
////////////////////////////////////////////////////

int print_wvectors(std::string filename)
{
    std::string wvectors_name(filename);
    wvectors_name.append(".wvectors");
    std::ofstream wvector_file;
    wvector_file.open(wvectors_name.c_str());
    
    for(int c=0; c<number_of_particles; c++)
    {
        std::vector<int> extended_wvector = all_wvectors[c];
        int  length = extended_wvector.back();                          // LENGTH OF WVECTOR
        int plength = extended_wvector.end()[-2];//back()-1;            // LENGTH OF PVECTOR
        
        wvector_file << extended_wvector[length+plength-3] << '\t';     // NUMBER OF FACES
        
        wvector_file << '(';                                            // P VECTOR
        for(int c=length; c<length+plength-4; c++)
            wvector_file << extended_wvector[c] << ',';
        wvector_file << extended_wvector[length+plength-4] << ')' << '\t';
        
        wvector_file << '(';                                            // WEINBERG VECTOR
        for(int c=0; c<length; c++)
            wvector_file << extended_wvector[c] << ',';
        wvector_file << 1 << ')' << '\t';
        
        wvector_file << extended_wvector[length+plength-2] << '\t';     // SYMMETRIES
        wvector_file << extended_wvector[length+plength-1] << '\t';     // CHIRALITY
        wvector_file << extended_wvector[length+plength]   << '\t';     // STABLE
        
        wvector_file << '\n';
    }
    
    return 0;
}





