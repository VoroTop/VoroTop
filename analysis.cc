////////////////////////////////////////////////////////
////                                                ////
////   ******************************************   ////
////   *                                        *   ////
////   *     VoroTop: Voronoi Cell Topology     *   ////
////   *   Visualization and Analysis Toolkit   *   ////
////   *             (Version 0.4)              *   ////
////   *                                        *   ////
////   *           Emanuel A. Lazar             *   ////
////   *      University of Pennsylvania        *   ////
////   *           December 5, 2014             *   ////
////   *                                        *   ////
////   ******************************************   ////
////                                                ////
////////////////////////////////////////////////////////

////   File: analysis.cc


#include <queue>
#include <random>
#include <iostream>

#include "filters.hh"
#include "variables.hh"



bool vector_size(const std::vector<int>& a,const std::vector<int>& b)
{
    return (a.size() > b.size());
}


////////////////////////////////////////////////////
////
////   DETERMINE "CLUSTERS" OF NON-CRYSTALLINE PARTICLES;
////   BY DEFAULT, ONLY TYPES NOT APPEARING IN A FILTER
////   ARE CLASSIFIED AS DEFECTS AND ARE CLUSTERED.  IF
////   AN OPTIONAL STRUCTURE TYPE IS PROVIDED AT THE COMMAND
////   LINE THEN ALL PARTICLES WITH OTHER STRUCTURE TYPES ARE
////   CONSIDERED DEFECTS FOR THIS CLUSTERING PURPOSE.
////
////////////////////////////////////////////////////

void   cluster_analysis(Filter &filter)
{
    int total_defect_particles  = 0;
    int total_crystal_particles = 0;
    
    std::vector <std::vector <int> > clusters_good;
    std::vector <std::vector <int> > clusters_bad;

    // IF OPTION TO RESOLVE INDETERMINATE TYPES IS USED,
    // THEN THE RESULTS OF THAT ANALYSIS ARE USED HERE.
    if(clustering_default_switch == 0)
    {
        if(r_switch) for(int c=0; c<number_of_particles; c++) cluster_index[c] = resolved_types    [c];
        else         for(int c=0; c<number_of_particles; c++) cluster_index[c] = vt_structure_types[c];
    }
    else
    {
        if(r_switch)
        {
            for(int c=0; c<number_of_particles; c++)
            {
                if(resolved_types[c]==clustering_default)     cluster_index[c] = 1;
                else                                          cluster_index[c] = 0;
            }
        }
        else
        {
            for(int c=0; c<number_of_particles; c++)
            {
                if(vt_structure_types[c]==clustering_default) cluster_index[c] = 1;
                else                                          cluster_index[c] = 0;
            }
        }
    }
    
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
    
    // ASSIGN AND RECORD FOR EACH PARTICLE A CLUSTER ID. NEGATIVE
    // INDICES INDICATE DEFECT CLUSTERS; POSITIVE INDICES INDICATE
    // CRYSTALLINE CLUSTERS. ALSO RECORD FOR EACH PARTICLE THE
    // NUMBER OF PARTICLES IN ITS CLUSTER.
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
    /*
    int max_cluster_size = 100;
    std::vector<int> c_sizes(max_cluster_size,0);
    double sum_of_squares = 0;
    
    for(int c=0; c<clusters_bad.size(); c++)
    {
        if(clusters_bad[c].size()<max_cluster_size)
            c_sizes[clusters_bad[c].size()]++;
        sum_of_squares += clusters_bad[c].size()*clusters_bad[c].size();
    }
    
    std::cout << name_of_data_file      << '\t';
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
    */
}



////////////////////////////////////////////////////
////
////   COMPUTES DISTRIBUTION OF TYPES IN SYSTEM
////   AFTER RANDOM GAUSSIAN PERTURBATIONS.
////   COMPUTES perturbation_samples RANDOM
////   PERTURBATIONS OF SYSTEM, EACH OF STANDARD
////   DEVIATION perturbation_size, AND SUMS THE
////   DISTRIBUTION OVER ALL SYSTEMS.  THIS CAN
////   BE USEFUL IN COMPUTING TYPES ASSOCIATED TO
////   REALISTIC VERSIONS OF IDEAL SYSTEMS.
////
////////////////////////////////////////////////////

void calc_gaussian_distribution(voro::container_periodic& con, voro::particle_order& vo, Filter &filter)
{
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0., perturbation_size);
    
    for(int loop = 0; loop < perturbation_samples; loop++)
    {
        // NEW, PERTURBED VERSION OF THE PRIMARY SYSTEM
        voro::particle_order voP;
        voro::container_periodic conP (supercell_edges[0][0],supercell_edges[1][0],supercell_edges[1][1],
                                       supercell_edges[2][0],supercell_edges[2][1],supercell_edges[2][2],
                                       n_x,n_y,n_z,8);
        
        voro::c_loop_order_periodic vlo(con,vo);
        
        // ADD PERTURBATION OF PARTICLES IN PRIMARY CONTAINER TO PERTURBED CONTAINER
        if(vlo.start()) do
        {
            int pid = vlo.id[vlo.ijk][vlo.q];
            double x,y,z;
            vlo.pos(x,y,z);
            
            x += distribution(generator);
            y += distribution(generator);
            z += distribution(generator);
            
            // ADD ``PERTURBED'' COPY
            conP.put(voP,pid,x,y,z);
            
            // WE SHOULD FIX THIS FOR PERIODIC BOUNDARY CONDITIONS; DIFFICULT FOR
            // NON-ORTHORHOMBIC CELLS
            
        } while(vlo.inc());
        
        
        // COMPUTE DISTRIBUTION OF TYPES IN THE PERTURBED VERSION
        voro::voronoicell vcell;
        voro::c_loop_order_periodic vloP(conP,voP);
        if(vloP.start()) do if(conP.compute_cell(vcell,vloP))
            filter.increment_or_add(calc_wvector(vcell),1);
        while(vloP.inc());
    }

    filter.relabel_data_types();
}



////////////////////////////////////////////////////
////
////   RESOLVES INDETERMINATE TYPES.  FOR NOW, WE
////   TAKE TYPE 2, AND SEE IF THEY'RE CLOSER TO 1
////   OR 3.  NOT SURE HOW WE SHOULD EXECUTE THIS.
////
////////////////////////////////////////////////////

void resolve_indeterminate_types(voro::container_periodic& con, voro::particle_order& vo, Filter &filter)
{
    std::random_device r;
    
    std::vector <int> closer_to_resolved_primary(number_of_particles,0);
    std::vector <int> closer_to_resolved_secondary(number_of_particles,0);

    // CHOICES OF SAMPLES AND SIZE ARE LARGELY ARBITRARY
    std::default_random_engine generator(r());
    std::normal_distribution<double> distribution(0., 0.05);
    
    std::cout<< "Resolving indeterminate types\n";
    
    // THIS PART OF THE LOOP SHOULD BE PARALLELIZED, SINCE EACH
    // CAN BE MADE INDEPENDENTLY
    for(int loop = 0; loop < resolve_trials; loop++)
    {
        std::cout << "Loop " << loop+1 << " of " << resolve_trials << '\n';
        // NEW, PERTURBED VERSION OF THE PRIMARY SYSTEM
        voro::particle_order voP;
        voro::container_periodic conP (supercell_edges[0][0],supercell_edges[1][0],supercell_edges[1][1],
                                       supercell_edges[2][0],supercell_edges[2][1],supercell_edges[2][2],
                                       n_x,n_y,n_z,8);
        
        voro::c_loop_order_periodic vlo(con,vo);
        
        // ADD PERTURBATION OF PARTICLES IN PRIMARY CONTAINER TO PERTURBED CONTAINER
        if(vlo.start()) do
        {
            int pid = vlo.id[vlo.ijk][vlo.q];
            double x,y,z;
            vlo.pos(x,y,z);
            
            double psize = pow(volumes[pid],1./3.);
            
            x += distribution(generator)*psize;
            y += distribution(generator)*psize;
            z += distribution(generator)*psize;
            
            // ADD ``PERTURBED'' COPY
            conP.put(voP,pid,x,y,z);
        } while(vlo.inc());
        
        // COMPUTE DISTRIBUTION OF TYPES IN THE PERTURBED VERSION
        voro::voronoicell vcell;
        voro::c_loop_order_periodic vloP(conP,voP);
        
        if(vloP.start()) do
        {
            int pid = vloP.id[vloP.ijk][vloP.q];
            
            if(filter.is_indeterminate(vt_structure_types[pid]))
            {
                if(conP.compute_cell(vcell,vloP))
                {
                    int pid = vloP.id[vloP.ijk][vloP.q];
                    int ntype = filter.wvector_type(calc_wvector(vcell));

                    if     (ntype == filter.resolved_types[vt_structure_types[pid]].first)  closer_to_resolved_primary  [pid]++;
                    else if(ntype == filter.resolved_types[vt_structure_types[pid]].second) closer_to_resolved_secondary[pid]++;
                }
            }
        }
        while(vloP.inc());
    }
    
    for(int c=0; c<number_of_particles; c++) if(filter.is_indeterminate(vt_structure_types[c]))
    {
        if(closer_to_resolved_secondary[c]>closer_to_resolved_primary[c]) resolved_types[c] = filter.resolved_types[vt_structure_types[c]].second;
        else resolved_types[c] = filter.resolved_types[vt_structure_types[c]].first;
    }
    else
        resolved_types[c] = vt_structure_types[c];
}

    
void calc_structure_types(Filter &filter)
{
    for(int c=0; c<number_of_particles; c++)
        vt_structure_types[c] = filter.wvector_type(all_wvectors[c]);
}


