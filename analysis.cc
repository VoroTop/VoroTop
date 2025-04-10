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

////   File: analysis.cc


#include <set>
#include <queue>
#include <random>
#include <iostream>

#include "filters.hh"
#include "variables.hh"
#include "functions.hh"


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

void cluster_analysis2()
{
    int total_defect_particles  = 0;
    int total_crystal_particles = 0;
    
    std::vector <std::vector <int>> clusters_good;
    std::vector <std::vector <int>> clusters_bad;

    // IF OPTION TO RESOLVE INDETERMINATE TYPES IS USED,
    // THEN THE RESULTS OF THAT ANALYSIS ARE USED HERE.
    if(clustering_default_switch == 0)
    {
        if(r_switch) for(int c=0; c<number_of_particles; c++) cluster_index[c] = vt_structure_types_resolved[c];
        else         for(int c=0; c<number_of_particles; c++) cluster_index[c] = vt_structure_types[c];
    }
    else
    {
        if(r_switch)
        {
            for(int c=0; c<number_of_particles; c++)
            {
                if(vt_structure_types_resolved[c]==clustering_default)     cluster_index[c] = 1;
                else                                                       cluster_index[c] = 0;
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
                
                int neighbor_count = cell_neighbor_count[w];
                for(int e=0; e<neighbor_count; e++)
                {
                    if(visited[neighbors_list_char[w][e]]==0)
                    {
                        visited[neighbors_list_char[w][e]]=1;
                        q.push(neighbors_list_char[w][e]);
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
                
                int neighbor_count = cell_neighbor_count[w];
                for(int e=0; e<neighbor_count; e++)
                {
                    if(visited[neighbors_list_char[w][e]]==0)
                    {
                        visited[neighbors_list_char[w][e]]=1;
                        q.push(neighbors_list_char[w][e]);
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
    for(unsigned int c=0; c<clusters_bad.size(); c++)
    {
        int size = clusters_bad[c].size();
        for(int d=0; d<size; d++)
        {
            cluster_index[clusters_bad[c][d]]=-c-1;
            cluster_sizes[clusters_bad[c][d]]=size;
        }
    }
    
    for(unsigned int c=0; c<clusters_good.size(); c++)
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
    
    for(unsigned int c=0; c<clusters_bad.size(); c++)
    {
        int cluster_size = clusters_bad[c].size();
        if(cluster_size < max_cluster_size)
            c_sizes[cluster_size]++;
        sum_of_squares += cluster_size * cluster_size;
    }
    
    std::cout << filename_data          << '\t';
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
////   DISTRIBUTION OF TYPES IN RANDOM PERTURBATIONS 
////   OF A GIVEN 2D SYSTEM.  USEFUL FOR ANALYZING 
////   REALISTIC VERSIONS OF IDEAL SYSTEMS.
////
////////////////////////////////////////////////////

void calc_gaussian_distribution_2d(Filter &filter)
{
    std::mt19937 generator(std::random_device{}());
    std::normal_distribution<double> distribution(0., perturbation_size);
    
    // NEW, PERTURBED VERSION OF THE ORIGINAL SYSTEM
    voro::container_2d con2d_perturbed (origin[0],hi_bound[0],origin[1],hi_bound[1],n_x,n_y,true,true,4,threads);
        
    for(int loop = 0; loop < perturbation_samples; loop++)
    {
        for(int c=0; c<number_of_particles; c++)
        {
            double x = particle_coordinates[2*c]   + distribution(generator);
            double y = particle_coordinates[2*c+1] + distribution(generator);
            con2d_perturbed.put(c,x,y);
        }
        
        count_and_store_neighbors_2d(con2d_perturbed);
        calc_distribution_2d(filter);
        con2d_perturbed.clear();
    }
}


////////////////////////////////////////////////////
////
////   DISTRIBUTION OF TYPES IN RANDOM PERTURBATIONS 
////   OF A GIVEN 3D SYSTEM.  USEFUL FOR ANALYZING 
////   REALISTIC VERSIONS OF IDEAL SYSTEMS.
////
////////////////////////////////////////////////////

void calc_gaussian_distribution_3d(Filter &filter) {
    std::mt19937 generator(std::random_device{}());  // Seed with truly random seed
    std::normal_distribution<double> distribution(0., perturbation_size);

    // Create a container for perturbed particles
    voro::container_3d con3d_perturbed(supercell_edges[0][0], supercell_edges[1][0], supercell_edges[1][1],
                                       supercell_edges[2][0], supercell_edges[2][1], supercell_edges[2][2],
                                       n_x, n_y, n_z, true, true, true, 8, threads);

    for (int loop = 0; loop < perturbation_samples; loop++) {
        // Clear the container for new perturbed positions
        con3d_perturbed.clear();

        // INSERT PERTURBED PARTICLE POSITIONS INTO CONTAINER
        for (int c = 0; c < number_of_particles; c++) {
            double x = particle_coordinates[3 * c] + distribution(generator);
            double y = particle_coordinates[3 * c + 1] + distribution(generator);
            double z = particle_coordinates[3 * c + 2] + distribution(generator);
            con3d_perturbed.put(c, x, y, z);
        }

        std::vector<Filter> local_filter(threads);

        #pragma omp parallel for num_threads(threads)
        for (auto cli = con3d_perturbed.begin(); cli < con3d_perturbed.end(); ++cli) {
            voro::voronoicell_3d vcell;
            if (con3d_perturbed.compute_cell(vcell, cli)) {
                std::vector<int> canonical_code;  // CANONICAL CODE WILL BE STORED HERE
                int chirality = compute_canonical_code_3d(canonical_code, vcell);

                int tid = omp_get_thread_num();
                local_filter[tid].increment_or_add(canonical_code, chirality, 1);
            }
        }

        for (int tid = 0; tid < threads; tid++) {
            filter.copy_filter(local_filter[tid]);
        }
    }
}




void pair_correlation_analysis(void)
{
    validate_max_radius(max_radius);
    // HOW FAR AWAY FROM CENTRAL PARTICLE; DEFAULT
    // SETTING IS 20; MAX VALUE IS 127.
    
    // FOR NORMALIZING COUNTS
    std::vector<double>                 normalization     (max_radius+1);
    if(u_switch==1)
        for(int k=1; k<=max_radius; k++)
            normalization[k]=1.;
    
    else
    {
        if(dimension==2)
        {
            double c0=18.77;
            double c1=12.5838;
            double c2=-0.488085;
            double c3=-24.8677;
            
            for(int k=1; k<=max_radius; k++)
                normalization[k] = c0+c1*k + c2*pow(double(k),0.5) + c3*pow(double(k),0.25);
        }
        if(dimension==3)
        {
            normalization[1]=15.535534;
            normalization[2]=69.95265;
            normalization[3]=183.9172;
            normalization[4]=366.30152;
            normalization[5]=621.8936;
            normalization[6]=953.61698;
            normalization[7]=1363.4318;
            normalization[8]=1852.768;
            normalization[9]=2422.6866;
            normalization[10]=3074.022;
            normalization[11]=3807.453;
            normalization[12]=4623.514;
            normalization[13]=5522.7118;
            normalization[14]=6505.4498;
            normalization[15]=7572.0404;
            normalization[16]=8722.7828;
            normalization[17]=9957.9538;
            normalization[18]=11277.806;
            normalization[19]=12682.522;
            normalization[20]=14172.318;
        }
    }

    // TWO-DIMENSIONAL ARRAYS FOR TALLYING DATA OVER THREADS
    vector<vector<unsigned long long> > neighbor_counts   (threads, vector<unsigned long long> (max_radius+1, 0));
    vector<vector<unsigned long long> > neighbor_counts_sq(threads, vector<unsigned long long> (max_radius+1, 0));

    // ARRAYS FOR TALLYING DATA OVER ALL PARTICLES
    std::vector<unsigned long long int> total_neighbor_counts   (max_radius+1);
    std::vector<unsigned long long int> total_neighbor_counts_sq(max_radius+1);

    // FOR EACH PARTICLE, COMPUTE GROWING CLUSTER AROUND IT
    // COUNTING ITS NUMBER OF NEIGHBORS IN EACH "RING"
    vector<vector<unsigned int> > cluster   (threads, vector<unsigned int> (number_of_particles,  0));
    vector<vector<signed char > > visited   (threads, vector<signed char > (number_of_particles, -1));
    vector<vector<int         > > kneighbors(threads, vector<int         > (max_radius+1, 0));

#pragma omp parallel for num_threads(threads)
    for(int p=0; p<number_of_particles; p++)
    {
        int tid=omp_get_thread_num();
        for(int q=0; q<number_of_particles; q++)
            cluster[tid][q]=visited[tid][q]=0;
        for(int q=0; q<=max_radius; q++)
            kneighbors[tid][q]=0;
        
        std::fill(cluster[tid].begin(), cluster[tid].end(),  0);
        std::fill(visited[tid].begin(), visited[tid].end(), -1);

        // visited TRACKS WHICH PARTICLES HAVE BEEN INCLUDED.
        // INITIALLIZED WITH -1, AND THEN UPDATED TO DISTANCE k.
        int counter=0;
        visited[tid][p]=0;
        cluster[tid][counter++]=p;
        kneighbors[tid][0]=1;
        
        int begin= 0;
        int end  = 1;
        
        // BUILD EACH RING, ADDING UNVISITED NEIGHBORS OF PRIOR RING
        for(int k=1; k<=max_radius; k++)
        {
            // ITERATE OVER ALL PARTICLES IN PRIOR RING
            for(int c=begin; c<end; c++)
            {
                // ID OF CURRENT PARTICLE
                int tempc      = cluster[tid][c];
                int nneighbors = neighbors_list_char[tempc].size();
                
                // ITERATE OVER ALL ITS NEIGHBORS
                for(int d=0; d<nneighbors; d++)
                {
                    if(visited[tid][neighbors_list_char[tempc][d]] == -1)
                    {
                        visited[tid][neighbors_list_char[tempc][d]]=k;
                        cluster[tid][counter++]=neighbors_list_char[tempc][d];
                        kneighbors[tid][k]++;
                    }
                }
            }
            
            begin = end;
            end   = end+kneighbors[tid][k];
            
            neighbor_counts[tid][k]    += kneighbors[tid][k];
            neighbor_counts_sq[tid][k] += kneighbors[tid][k]*kneighbors[tid][k];
        }
    }

    for(int tid=0; tid<threads; tid++)
        for(int k=0; k<max_radius+1; k++)
    {
        total_neighbor_counts   [k] += neighbor_counts   [tid][k];
        total_neighbor_counts_sq[k] += neighbor_counts_sq[tid][k];
    }
    ////////////////////////////////////////////////////
    ////
    ////   OUTPUT BASIC ANALYSIS OF VPCF DATA
    ////
    ////////////////////////////////////////////////////
    
    double total=1;         // WILL BE USED TO COMPUTE WHAT FRACTION OF SYSTEM IS COVERED
                            // BY CLUSTERS OF RADIUS k.  THIS NUMBER SHOULD NOT BE CLOSE
                            // TO 1, OR ELSE THIS INDICATES "WRAP-AROUND" EFFECTS.
    
    for(int k=1; k<=max_radius; k++)
    {
        double avg =      (double(total_neighbor_counts[k])   /double(number_of_particles))/normalization[k];
        double variance = (double(total_neighbor_counts_sq[k])/double(number_of_particles) - avg*avg*normalization[k]*normalization[k])/normalization[k]/normalization[k];
        if(abs(variance)<1e-15) variance=0;
        double std = sqrt(variance);
        
        total += double(total_neighbor_counts[k]) / double(number_of_particles);
        std::cout << k << '\t' << avg << '\t' << variance << '\t' << std << '\t' << total/double(number_of_particles) << '\n';
    }
    std::cout << '\n';
}



void pair_correlation_analysis2(void)
{
    // HOW FAR AWAY FROM CENTRAL PARTICLE; DEFAULT
    // SETTING IS 20; MAX VALUE IS 127.
    
    // FOR NORMALIZING COUNTS
    std::vector<double>                 normalization     (max_radius+1);
    if(u_switch==1)
        for(int k=1; k<=max_radius; k++)
            normalization[k]=1.;
    
    else
    {
        if(dimension==2)
        {
            double c0=18.77;
            double c1=12.5838;
            double c2=-0.488085;
            double c3=-24.8677;
            
            for(int k=1; k<=max_radius; k++)
                normalization[k] = c0+c1*k + c2*pow(double(k),0.5) + c3*pow(double(k),0.25);
        }
        if(dimension==3)
        {
            normalization[1]=15.535534;
            normalization[2]=69.95265;
            normalization[3]=183.9172;
            normalization[4]=366.30152;
            normalization[5]=621.8936;
            normalization[6]=953.61698;
            normalization[7]=1363.4318;
            normalization[8]=1852.768;
            normalization[9]=2422.6866;
            normalization[10]=3074.022;
            normalization[11]=3807.453;
            normalization[12]=4623.514;
            normalization[13]=5522.7118;
            normalization[14]=6505.4498;
            normalization[15]=7572.0404;
            normalization[16]=8722.7828;
            normalization[17]=9957.9538;
            normalization[18]=11277.806;
            normalization[19]=12682.522;
            normalization[20]=14172.318;
        }
    }
    

    // TWO-DIMENSIONAL ARRAYS FOR TALLYING DATA OVER THREADS
    vector<vector<unsigned long long> > neighbor_counts   (threads, vector<unsigned long long> (max_radius+1, 0));
    vector<vector<unsigned long long> > neighbor_counts_sq(threads, vector<unsigned long long> (max_radius+1, 0));

    // FOR EACH PARTICLE, COMPUTE GROWING CLUSTER AROUND IT
    // COUNTING ITS NUMBER OF NEIGHBORS IN EACH "RING"

#pragma omp parallel for num_threads(threads)
    for(int p=0; p<number_of_particles; p++)
    {
        int tid=omp_get_thread_num();
        
        std::set <int> ring0;
        std::set <int> ring1;
        std::set <int> ring2;

        ring0.insert(p);
        for(int s=0; s<cell_neighbor_count[p]; s++)
            ring1.insert(neighbors_list_char[p][s]);
        
        for(int k=0; k<=max_radius; k++)
        {
            // MOVE RING1 TO RING0, MOVE RING2 TO RING1, AND REPEAT
            // BUILD EACH RING, ADDING UNVISITED NEIGHBORS OF PRIOR RING
            neighbor_counts[tid][k]    += ring0.size();
            neighbor_counts_sq[tid][k] += ring0.size()*ring0.size();

            
            std::set <int> big;
            for(std::set<int>::iterator it=ring1.begin(); it != ring1.end(); it++)
                for(int s=0; s<cell_neighbor_count[*it]; s++)
                    big.insert(neighbors_list_char[*it][s]);
            ring0.insert(ring1.begin(),ring1.end());
            
            std::set_difference(big.begin(), big.end(), ring0.begin(), ring0.end(),
                std::inserter(ring2, ring2.end()));

            ring0 = std::move(ring1);
            ring1 = std::move(ring2);
            ring2.clear();
        }
    }

    // ARRAYS FOR TALLYING DATA OVER ALL PARTICLES
    std::vector<unsigned long long int> total_neighbor_counts   (max_radius+1);
    std::vector<unsigned long long int> total_neighbor_counts_sq(max_radius+1);

    for(int tid=0; tid<threads; tid++)
        for(int k=0; k<max_radius+1; k++)
    {
        total_neighbor_counts   [k] += neighbor_counts   [tid][k];
        total_neighbor_counts_sq[k] += neighbor_counts_sq[tid][k];
    }
    ////////////////////////////////////////////////////
    ////
    ////   OUTPUT BASIC ANALYSIS OF VPCF DATA
    ////
    ////////////////////////////////////////////////////
    
    
    double total=1;         // WILL BE USED TO COMPUTE WHAT FRACTION OF SYSTEM IS COVERED
    // BY CLUSTERS OF RADIUS k.  THIS NUMBER SHOULD NOT BE CLOSE
    // TO 1, OR ELSE THIS INDICATES "WRAP-AROUND" EFFECTS.
    
    for(int k=1; k<=max_radius; k++)
    {
        double avg =      (double(total_neighbor_counts[k])   /double(number_of_particles))/normalization[k];
        double variance = (double(total_neighbor_counts_sq[k])/double(number_of_particles) - avg*avg*normalization[k]*normalization[k])/normalization[k]/normalization[k];
        double std = sqrt(variance);
        
        total += double(total_neighbor_counts[k]) / double(number_of_particles);
        std::cout << k << '\t' << avg << '\t' << variance << '\t' << std << '\t' << total/double(number_of_particles) << '\n';
    }
    std::cout << '\n';
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

void cluster_analysis(void)
{
    const int DEFAULT_MAX_CLUSTER_SIZE = 100;

    int total_defect_particles  = 0;
    int total_crystal_particles = 0;
    
    std::vector <std::vector <int>> clusters_good;
    std::vector <std::vector <int>> clusters_bad;
    
    
    if(clustering_default_switch == 0)
        for(int c=0; c<number_of_particles; c++) cluster_index[c] = vt_structure_types[c];
    
    else
    {
        for(int c=0; c<number_of_particles; c++)
        {
            if(vt_structure_types[c]==clustering_default) cluster_index[c] = 1;
            else                                          cluster_index[c] = 0;
        }
    }
    
    // BUILD DEFECT CLUSTERS IN O(N) TIME; FOR USE
    // IN SYSTEMS THAT ARE PRIMARILY CRYSTALLINE.
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
                
                int neighbor_count = cell_neighbor_count[w];
                for(int e=0; e<neighbor_count; e++)
                {
                    if(visited[neighbors_list_char[w][e]]==0)
                    {
                        visited[neighbors_list_char[w][e]]=1;
                        q.push(neighbors_list_char[w][e]);
                    }
                }
            }
            clusters_bad.push_back(cluster);
        }
    }
    sort(clusters_bad.begin(),clusters_bad.end(),vector_size);
    
    
    // BUILD CRYSTALLINE CLUSTERS IN O(N) TIME; FOR USE
    // IN SYSTEMS THAT ARE PRIMARILY DISORDERED.
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

                int neighbor_count = cell_neighbor_count[w];
                for(int e=0; e<neighbor_count; e++)
                {
                    if(visited[neighbors_list_char[w][e]]==0)
                    {
                        visited[neighbors_list_char[w][e]]=1;
                        q.push(neighbors_list_char[w][e]);
                    }
                }
            }
            clusters_good.push_back(cluster);
        }
    }
    sort(clusters_good.begin(),clusters_good.end(),vector_size);
    
    
    ////////////////////////////////////////////////////
    // ASSIGN AND RECORD FOR EACH PARTICLE A CLUSTER ID, STORED IN
    // cluster_index[]. NEGATIVE INDICES INDICATE DEFECT CLUSTERS;
    // POSITIVE INDICES INDICATE CRYSTALLINE CLUSTERS. ALSO RECORD
    // FOR EACH PARTICLE THE NUMBER OF PARTICLES IN ITS CLUSTER,
    // STORED IN cluster_sizes[]
    for(unsigned int c=0; c<clusters_bad.size(); c++)
    {
        unsigned int size = clusters_bad[c].size();
        for(unsigned int d=0; d<size; d++)
        {
            cluster_index[clusters_bad[c][d]]=-c-1;
            cluster_sizes[clusters_bad[c][d]]=size;
        }
    }
    
    for(unsigned int c=0; c<clusters_good.size(); c++)
    {
        unsigned int size = clusters_good[c].size();
        for(unsigned int d=0; d<size; d++)
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
    
    unsigned int max_cluster_size = DEFAULT_MAX_CLUSTER_SIZE;
    std::vector<int> c_sizes(max_cluster_size,0);
    double sum_of_squares = 0;
    
    for(unsigned int c=0; c<clusters_bad.size(); c++)
    {
        if(clusters_bad[c].size()<max_cluster_size)
            c_sizes[clusters_bad[c].size()]++;
        sum_of_squares += clusters_bad[c].size()*clusters_bad[c].size();
    }
    
    std::cout << filename_data          << '\t';
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
    
    for(unsigned int c=1; c<max_cluster_size; c++)       // NUMBER OF CLUSTERS WITH 1, 2, 3, ..., max_cluster_size-1 PARTICLES
        std::cout << c_sizes[c] << '\t';
    
    std::cout << '\n';
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

void defect_cluster_analysis(void)  // EXPERIMENTAL
{
    int total_defect_particles  = 0;
    std::vector <std::vector <int>> clusters_bad;
    
    int typeA_vacancy_count=0;
    int typeB_vacancy_count=0;
    int typeC_vacancy_count=0;
    int dislocation_count  =0;
    int grain_boundary_count=0;
    
    // IF USING clustering_default_switch FEATURE, THEN WE WANT TO
    // GROUP PARTICULAR KINDS OF DEFECT PARTICLES TOGETHER.  FOR EXAMPLE
    // SIX VACANCY DEFECTS CAN BE COMBINED AS A VACANCY, ETC.
    if(clustering_default_switch == 0)
        for(int c=0; c<number_of_particles; c++)
        {
            cluster_index[c] = vt_structure_types[c];
            if(cluster_index[c]==1) cluster_index[c]=0;
        }
    
    else
    {
        for(int c=0; c<number_of_particles; c++)
        {
            if(vt_structure_types[c]==clustering_default) cluster_index[c] = 1;
            else                                          cluster_index[c] = 0;
        }
    }
    
    // BUILD DEFECT CLUSTERS IN O(N) TIME; FOR USE
    // IN SYSTEMS THAT ARE PRIMARILY CRYSTALLINE.
    std::vector <int> visited(number_of_particles,0);
    for(int d=0; d<number_of_particles; d++)
    {
        if(cluster_index[d]==2 || cluster_index[d]==3 || cluster_index[d]==4) total_defect_particles++;
        //if(cluster_index[d]==0) total_defect_particles++;
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

                int neighbor_count = cell_neighbor_count[w];
                for(int e=0; e<neighbor_count; e++)
                {
                    if(visited[neighbors_list_char[w][e]]==0)
                    {
                        visited[neighbors_list_char[w][e]]=1;
                        q.push(neighbors_list_char[w][e]);
                    }
                }
            }
            clusters_bad.push_back(cluster);
        }
    }
    sort(clusters_bad.begin(),clusters_bad.end(),vector_size);
    
    std::cout << "We have now " << clusters_bad.size() << " defect clusters\n";
    
    
    for(unsigned int c=0; c<clusters_bad.size(); c++)
    {
        bool all_type_3=1;
        bool all_type_4=1;
        bool no_type_2 =1;
        bool no_type_4 =1;
        
        unsigned int size = clusters_bad[c].size();
        std::cout << "Cluster with " << size << " particles " << '\n';
        for(unsigned int d=0; d<size; d++)
        {
            if(vt_structure_types[clusters_bad[c][d]]==2) no_type_2 =0;
            if(vt_structure_types[clusters_bad[c][d]]==4) no_type_4 =0;
            if(vt_structure_types[clusters_bad[c][d]]!=3) all_type_3=0;
            if(vt_structure_types[clusters_bad[c][d]]!=4) all_type_4=0;
            std::cout << clusters_bad[c][d] << '\t' << vt_structure_types[clusters_bad[c][d]] << '\n';
            
            //            cluster_index[clusters_bad[c][d]]=-c-1;
            //            cluster_sizes[clusters_bad[c][d]]=size;
        }
        if(size==6 && all_type_4==1) typeA_vacancy_count++;
        if(size==6 && no_type_2 ==1)
        {
            for(unsigned int d=0; d<size; d++)
                vt_structure_types[clusters_bad[c][d]]=4;
            typeB_vacancy_count++;
        }
        if(size==3 && all_type_4==1) typeC_vacancy_count++;
        if(size==2 && all_type_3 ==1) dislocation_count++;
        
        if(size>2  && no_type_4 ==1)
        {
            for(unsigned int d=0; d<size; d++)
                vt_structure_types[clusters_bad[c][d]]=2;
            grain_boundary_count++;
        }
        std::cout << '\n';
    }
    
    ////////////////////////////////////////////////////
    // ASSIGN AND RECORD FOR EACH PARTICLE A CLUSTER ID, STORED IN
    // cluster_index[]. NEGATIVE INDICES INDICATE DEFECT CLUSTERS;
    // POSITIVE INDICES INDICATE CRYSTALLINE CLUSTERS. ALSO RECORD
    // FOR EACH PARTICLE THE NUMBER OF PARTICLES IN ITS CLUSTER,
    // STORED IN cluster_sizes[]
    for(unsigned int c=0; c<clusters_bad.size(); c++)
    {
        unsigned int size = clusters_bad[c].size();
        for(unsigned int d=0; d<size; d++)
        {
            cluster_index[clusters_bad[c][d]]=-c-1;
            cluster_sizes[clusters_bad[c][d]]=size;
        }
    }
    
    
    ////////////////////////////////////////////////////
    ////
    ////   OUTPUT BASIC ANALYSIS OF CLUSTERING DATA TO SCREEN
    ////
    ////////////////////////////////////////////////////
    
    unsigned int max_cluster_size = 100;
    std::vector<int> c_sizes(max_cluster_size,0);
    double sum_of_squares = 0;
    
    for(unsigned int c=0; c<clusters_bad.size(); c++)
    {
        if(clusters_bad[c].size()<max_cluster_size)
            c_sizes[clusters_bad[c].size()]++;
        sum_of_squares += clusters_bad[c].size()*clusters_bad[c].size();
    }
    
    std::cout << filename_data          << '\t';
    std::cout << number_of_particles    << '\t';
    std::cout << total_defect_particles << '\t';
    std::cout << clusters_bad.size()    << '\t';
    
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
    
    for(unsigned int c=1; c<max_cluster_size; c++)       // NUMBER OF CLUSTERS WITH 1, 2, 3, ..., max_cluster_size-1 PARTICLES
        std::cout << c_sizes[c] << '\t';
    
    std::cout << '\n';
    std::cout << "We have " << typeA_vacancy_count << " type A vacncies\n";
    std::cout << "We have " << typeB_vacancy_count << " type B vacncies\n";
    std::cout << "We have " << typeC_vacancy_count << " type C vacncies\n";
    std::cout << "We have " << dislocation_count << " dislocations\n";
    std::cout << "We have " << grain_boundary_count << " grain boundaries\n";
    std::cout << '\n';
}



