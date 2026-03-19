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
#include <cmath>
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
////   DISTRIBUTION OF TYPES IN RANDOM PERTURBATIONS
////   OF A GIVEN 2D SYSTEM VIA CELL REPLICATION.
////
////   THE SYSTEM IS REPLICATED INTO A LARGE SUPERCELL
////   (AT LEAST 5x5) WITH GAUSSIAN PERTURBATIONS
////   APPLIED TO EACH PARTICLE.  THIS ENSURES NEIGHBOR
////   INDEPENDENCE FOR SMALL UNIT CELLS.
////
////////////////////////////////////////////////////

void calc_gaussian_distribution_2d(Filter &filter)
{
    int n = std::max(5, (int)ceil(sqrt((double)perturbation_samples / number_of_particles)));
    long long total_particles = (long long)number_of_particles * n * n;

    // MEMORY CHECK: EACH PARTICLE NEEDS ~16 BYTES FOR COORDINATES PLUS
    // OVERHEAD FOR CONTAINERS AND NEIGHBOR LISTS.  CAP AT ~500 MILLION.
    if(total_particles > 500000000LL)
    {
        std::cerr << "Error: supercell would contain " << total_particles
                  << " particles, which is too large for available memory.\n"
                  << "Reduce the target sample count or use a larger input system.\n";
        exit(1);
    }

    std::cout << "Perturbation sampling: " << number_of_particles << " particles in input, "
              << "target " << perturbation_samples << " samples" << std::endl;
    std::cout << "Replicating to " << n << "x" << n
              << " supercell (" << total_particles << " particles)" << std::endl;

    double Lx = xhi - xlo;
    double Ly = yhi - ylo;
    double ox = xlo, oy = ylo;

    std::mt19937 generator(std::random_device{}());
    std::normal_distribution<double> noise(0., perturbation_size);

    std::vector<double> rep_coords(total_particles * 2);
    int idx = 0;

    for(int ki = 0; ki < n; ki++)
    for(int kj = 0; kj < n; kj++)
    {
        double dx = ki * Lx;
        double dy = kj * Ly;

        for(int p = 0; p < number_of_particles; p++)
        {
            rep_coords[2*idx]     = (particle_coordinates[2*p]     - ox) + dx + noise(generator);
            rep_coords[2*idx + 1] = (particle_coordinates[2*p + 1] - oy) + dy + noise(generator);
            idx++;
        }
    }

    double new_Lx = n * Lx;
    double new_Ly = n * Ly;

    // GRID BLOCK SIZES FOR SUPERCELL
    int total_blocks = total_particles / 4 + 1;
    double lpb = pow(new_Lx * new_Ly / (double)total_blocks, 1./2.);
    int rep_nx = (int)(new_Lx / lpb) + 1;
    int rep_ny = (int)(new_Ly / lpb) + 1;

    // TEMPORARILY SWAP GLOBAL STATE FOR EXISTING 2D FUNCTIONS
    int orig_nop = number_of_particles;
    double* orig_coords = particle_coordinates;
    double orig_xlo = xlo, orig_xhi = xhi;
    double orig_ylo = ylo, orig_yhi = yhi;

    number_of_particles = total_particles;
    particle_coordinates = rep_coords.data();
    xlo = 0; xhi = new_Lx;
    ylo = 0; yhi = new_Ly;

    list_of_neighbors.resize(total_particles);
    cell_neighbor_count.resize(total_particles);

    voro::container_2d con_rep(0, new_Lx, 0, new_Ly,
                               rep_nx, rep_ny, true, true, 4, threads);
    for(int i = 0; i < total_particles; i++)
        con_rep.put(i, rep_coords[2*i], rep_coords[2*i+1]);

    count_and_store_neighbors_2d(con_rep);
    calc_distribution_2d(filter);

    // RESTORE GLOBAL STATE
    number_of_particles = orig_nop;
    particle_coordinates = orig_coords;
    xlo = orig_xlo; xhi = orig_xhi;
    ylo = orig_ylo; yhi = orig_yhi;
    list_of_neighbors.resize(orig_nop);
    cell_neighbor_count.resize(orig_nop);
}


////////////////////////////////////////////////////
////
////   DISTRIBUTION OF TYPES IN RANDOM PERTURBATIONS
////   OF A GIVEN 3D SYSTEM VIA CELL REPLICATION.
////
////   THE SYSTEM IS REPLICATED INTO A LARGE SUPERCELL
////   (AT LEAST 5x5x5) WITH GAUSSIAN PERTURBATIONS
////   APPLIED TO EACH PARTICLE.  THIS ENSURES NEIGHBOR
////   INDEPENDENCE FOR SMALL UNIT CELLS.
////
////////////////////////////////////////////////////

void calc_gaussian_distribution_3d(Filter &filter)
{
    int n = std::max(5, (int)ceil(pow((double)perturbation_samples / number_of_particles, 1.0/3.0)));
    long long total_particles = (long long)number_of_particles * n * n * n;

    // MEMORY CHECK: EACH PARTICLE NEEDS ~24 BYTES FOR COORDINATES PLUS
    // OVERHEAD FOR CONTAINERS AND TOPOLOGY COMPUTATION.  CAP AT ~500 MILLION.
    if(total_particles > 500000000LL)
    {
        std::cerr << "Error: supercell would contain " << total_particles
                  << " particles, which is too large for available memory.\n"
                  << "Reduce the target sample count or use a larger input system.\n";
        exit(1);
    }

    std::cout << "Perturbation sampling: " << number_of_particles << " particles in input, "
              << "target " << perturbation_samples << " samples" << std::endl;
    std::cout << "Replicating to " << n << "x" << n << "x" << n
              << " supercell (" << total_particles << " particles)" << std::endl;

    double Lx = xhi - xlo;
    double Ly = yhi - ylo;
    double Lz = zhi - zlo;

    // ORIGIN OFFSET FOR ORTHOGONAL SYSTEMS (TRICLINIC ALREADY AT ORIGIN)
    double ox = triclinic_crystal_system ? 0.0 : xlo;
    double oy = triclinic_crystal_system ? 0.0 : ylo;
    double oz = triclinic_crystal_system ? 0.0 : zlo;

    std::mt19937 generator(std::random_device{}());
    std::normal_distribution<double> noise(0., perturbation_size);

    std::vector<double> rep_coords(total_particles * 3);
    int idx = 0;

    for(int ki = 0; ki < n; ki++)
    for(int kj = 0; kj < n; kj++)
    for(int kk = 0; kk < n; kk++)
    {
        // SHIFT = ki*a + kj*b + kk*c (LATTICE VECTOR OFFSETS)
        double dx = ki * Lx + kj * xy + kk * xz;
        double dy =           kj * Ly + kk * yz;
        double dz =                     kk * Lz;

        for(int p = 0; p < number_of_particles; p++)
        {
            rep_coords[3*idx]     = (particle_coordinates[3*p]     - ox) + dx + noise(generator);
            rep_coords[3*idx + 1] = (particle_coordinates[3*p + 1] - oy) + dy + noise(generator);
            rep_coords[3*idx + 2] = (particle_coordinates[3*p + 2] - oz) + dz + noise(generator);
            idx++;
        }
    }

    // SUPERCELL BOX DIMENSIONS
    double new_Lx = n * Lx;
    double new_Ly = n * Ly;
    double new_Lz = n * Lz;
    double new_xy = n * xy;
    double new_xz = n * xz;
    double new_yz = n * yz;

    // GRID BLOCK SIZES FOR SUPERCELL
    int total_blocks = total_particles / 4 + 1;
    double lpb = pow(new_Lx * new_Ly * new_Lz / (double)total_blocks, 1./3.);
    int rep_nx = (int)(new_Lx / lpb) + 1;
    int rep_ny = (int)(new_Ly / lpb) + 1;
    int rep_nz = (int)(new_Lz / lpb) + 1;

    // CREATE CONTAINER, ADD PARTICLES, AND COMPUTE TOPOLOGIES
    auto compute_topologies = [&](auto& con_rep) {
        for(int i = 0; i < total_particles; i++)
            con_rep.put(i, rep_coords[3*i], rep_coords[3*i+1], rep_coords[3*i+2]);

        std::vector<Filter> local_filter(threads);

        #pragma omp parallel for num_threads(threads)
        for(auto cli = con_rep.begin(); cli < con_rep.end(); ++cli)
        {
            voro::voronoicell_3d vcell;
            if(con_rep.compute_cell(vcell, cli))
            {
                VoronoiTopology result = compute_canonical_code_3d(vcell);
                int tid = omp_get_thread_num();
                local_filter[tid].increment_or_add(result.canonical_code, result.chirality, 1);
            }
        }

        for(int tid = 0; tid < threads; tid++)
            filter.copy_filter(local_filter[tid]);
    };

    if(triclinic_crystal_system)
    {
        voro::container_triclinic con_rep(new_Lx, new_xy, new_Ly,
                                          new_xz, new_yz, new_Lz,
                                          rep_nx, rep_ny, rep_nz, 8, threads);
        compute_topologies(con_rep);
    }
    else
    {
        voro::container_3d con_rep(0, new_Lx, 0, new_Ly, 0, new_Lz,
                                   rep_nx, rep_ny, rep_nz,
                                   true, true, true, 8, threads);
        compute_topologies(con_rep);
    }
}




void pair_correlation_analysis(void)
{
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
            if(max_radius>20)
            {
                std::cout << "The normalization constants for Voronoi distances greater than 20 are not yet implemented.\n";
                std::cout << "Will compute the (normalized) Voronoi pair correlation function only up to Voronoi distance 20.\n";
                std::cout << "Please use -u for the unnormalized Voronoi pair correlation function, which can be computed for any maximal distance.\n";
                max_radius=20;
            }
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
                int nneighbors = list_of_neighbors[tempc].size();
                
                // ITERATE OVER ALL ITS NEIGHBORS
                for(int d=0; d<nneighbors; d++)
                {
                    if(visited[tid][list_of_neighbors[tempc][d]] == -1)
                    {
                        visited[tid][list_of_neighbors[tempc][d]]=k;
                        cluster[tid][counter++]=list_of_neighbors[tempc][d];
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
    // SETTING IS 10; MAX VALUE IS 127.
    
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
            ring1.insert(list_of_neighbors[p][s]);
        
        for(int k=0; k<=max_radius; k++)
        {
            // MOVE RING1 TO RING0, MOVE RING2 TO RING1, AND REPEAT
            // BUILD EACH RING, ADDING UNVISITED NEIGHBORS OF PRIOR RING
            neighbor_counts[tid][k]    += ring0.size();
            neighbor_counts_sq[tid][k] += ring0.size()*ring0.size();

            
            std::set <int> big;
            for(std::set<int>::iterator it=ring1.begin(); it != ring1.end(); it++)
                for(int s=0; s<cell_neighbor_count[*it]; s++)
                    big.insert(list_of_neighbors[*it][s]);
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
    int total_crystal_particles = 0;
    int total_defect_particles  = 0;
    
    std::vector <std::vector <int>> clusters_crystalline;
    std::vector <std::vector <int>> clusters_defect;
    
    for(int c=0; c<number_of_particles; c++) cluster_index[c] = vt_structure_types[c];
    
    
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
                    if(visited[list_of_neighbors[w][e]]==0)
                    {
                        visited[list_of_neighbors[w][e]]=1;
                        q.push(list_of_neighbors[w][e]);
                    }
                }
            }
            clusters_defect.push_back(cluster);
        }
    }
    sort(clusters_defect.begin(),clusters_defect.end(),vector_size);
    
    
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
                    if(visited[list_of_neighbors[w][e]]==0)
                    {
                        visited[list_of_neighbors[w][e]]=1;
                        q.push(list_of_neighbors[w][e]);
                    }
                }
            }
            clusters_crystalline.push_back(cluster);
        }
    }
    sort(clusters_crystalline.begin(),clusters_crystalline.end(),vector_size);
    
    
    ////////////////////////////////////////////////////
    // ASSIGN AND RECORD FOR EACH PARTICLE A CLUSTER ID, STORED IN
    // cluster_index[]. NEGATIVE INDICES INDICATE DEFECT CLUSTERS;
    // POSITIVE INDICES INDICATE CRYSTALLINE CLUSTERS. ALSO RECORD
    // FOR EACH PARTICLE THE NUMBER OF PARTICLES IN ITS CLUSTER,
    // STORED IN cluster_sizes[]
    for(unsigned int c=0; c<clusters_defect.size(); c++)
    {
        unsigned int size = clusters_defect[c].size();
        for(unsigned int d=0; d<size; d++)
        {
            cluster_index[clusters_defect[c][d]]=-c-1;
            cluster_sizes[clusters_defect[c][d]]=size;
        }
    }
    
    for(unsigned int c=0; c<clusters_crystalline.size(); c++)
    {
        unsigned int size = clusters_crystalline[c].size();
        for(unsigned int d=0; d<size; d++)
        {
            cluster_index[clusters_crystalline[c][d]]=c+1;
            cluster_sizes[clusters_crystalline[c][d]]=size;
        }
    }
    
    // THE PRIMARY CLUSTER ANALYSIS HAS BEEN PERFORMED; IF A LAMMPS 
    // DUMP FILE IS REQUESTED, OR ELSE IF AN EPS FILE IS REQUESTED, 
    // THEN WE DO NOT OUTPUT THIS DATA, AND INSTEAD USE IT FOR THE 
    // LAMMPS DUMP FILE OR EPS FILE.  IF NEITHER OF THESE IS REQUESTED,
    // THEN WE OUTPUT THE DATA TO THE SCREEN.
    if(l_switch==1 || e_switch==1) return;
    
    
    ////////////////////////////////////////////////////
    ////
    ////   OUTPUT BASIC ANALYSIS OF CLUSTERING DATA TO SCREEN
    ////
    ////////////////////////////////////////////////////
    
    double sum_of_squares_defect = 0;    
    std::map<int, int> clusters_defect_countMap;
    std::map<int, int> clusters_crystal_countMap;
    
    for(unsigned int c=0; c<clusters_defect.size(); c++)
    {
       sum_of_squares_defect += clusters_defect[c].size()*clusters_defect[c].size();
       clusters_defect_countMap[clusters_defect[c].size()]++;
    }

    double sum_of_squares_crystal = 0;    
    for(unsigned int c=0; c<clusters_crystalline.size(); c++)
    {
        sum_of_squares_crystal += clusters_crystalline[c].size()*clusters_crystalline[c].size();
        clusters_crystal_countMap[clusters_crystalline[c].size()]++;
    }
   
    std::cout << '\n';
    std::cout << "File:                " << filename_data          << '\n';
    std::cout << "                     " << '\n';
    std::cout << "Number of particles: " << number_of_particles    << '\n';
    std::cout << "Crystal particles:   " << total_crystal_particles << '\n';
    std::cout << "Crystal clusters:    " << clusters_crystalline.size()   << '\n';
    std::cout << "Defect particles:    " << total_defect_particles << '\n';
    std::cout << "Defect clusters:     " << clusters_defect.size()    << '\n';
    std::cout << "                     " << '\n';

    if(!clusters_crystalline.empty())
    {
        std::cout << "Crystalline clusters " << '\n';
        std::cout << " Average size:       " << (double(total_crystal_particles) / clusters_crystalline.size()) << '\n';
        std::cout << " Standard deviation: " << sqrt(sum_of_squares_crystal / clusters_crystalline.size() - (total_crystal_particles / clusters_crystalline.size()) * (total_crystal_particles / clusters_crystalline.size())) << '\n';
        std::cout << " Smallest cluster:   " << (clusters_crystalline.back()).size() << '\n';
        std::cout << " Largest cluster:    " << (clusters_crystalline.front()).size() << '\n';
        std::cout << "                     " << '\n';

        std::cout << " Size \t Count " << '\n';
        std::cout << " ===================" << '\n';
        for (const auto& pair : clusters_crystal_countMap) {
            std::cout << " " << pair.first << '\t' << pair.second << '\n';
        }
        std::cout << '\n';    
    }
    else
    {
        std::cout << "No crystal clusters found." << '\n';
        std::cout << '\n';
    }
        
    if(!clusters_defect.empty())
    {
        std::cout << "Defect clusters      " << '\n';
        std::cout << " Average size:       " << (double(total_defect_particles) / clusters_defect.size()) << '\n';
        std::cout << " Standard deviation: " << sqrt(sum_of_squares_defect / clusters_defect.size() - (total_defect_particles / clusters_defect.size()) * (total_defect_particles / clusters_defect.size())) << '\n';
        std::cout << " Smallest cluster:   " << (clusters_defect.back()).size() << '\n';
        std::cout << " Largest cluster:    " << (clusters_defect.front()).size() << '\n';
        std::cout << "                     " << '\n';

        std::cout << " Size \t Count " << '\n';
        std::cout << " ===================" << '\n';
        for (const auto& pair : clusters_defect_countMap) {
            std::cout << " " << pair.first << '\t' << pair.second << '\n';
        }
        std::cout << '\n';    
    }
    else
        std::cout << "No defect clusters found." << '\n';
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
    std::vector <std::vector <int>> clusters_defect;
    
    int typeA_vacancy_count=0;
    int typeB_vacancy_count=0;
    int typeC_vacancy_count=0;
    int dislocation_count  =0;
    int grain_boundary_count=0;
    
    // IF USING clustering_default_switch FEATURE, THEN WE WANT TO
    // GROUP PARTICULAR KINDS OF DEFECT PARTICLES TOGETHER.  FOR EXAMPLE
    // SIX VACANCY DEFECTS CAN BE COMBINED AS A VACANCY, ETC.
    for(int c=0; c<number_of_particles; c++)
    {
        cluster_index[c] = vt_structure_types[c];
        if(cluster_index[c]==1) cluster_index[c]=0;
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
                    if(visited[list_of_neighbors[w][e]]==0)
                    {
                        visited[list_of_neighbors[w][e]]=1;
                        q.push(list_of_neighbors[w][e]);
                    }
                }
            }
            clusters_defect.push_back(cluster);
        }
    }
    sort(clusters_defect.begin(),clusters_defect.end(),vector_size);
    
    std::cout << "We have now " << clusters_defect.size() << " defect clusters\n";
    
    
    for(unsigned int c=0; c<clusters_defect.size(); c++)
    {
        bool all_type_3=1;
        bool all_type_4=1;
        bool no_type_2 =1;
        bool no_type_4 =1;
        
        unsigned int size = clusters_defect[c].size();
        std::cout << "Cluster with " << size << " particles " << '\n';
        for(unsigned int d=0; d<size; d++)
        {
            if(vt_structure_types[clusters_defect[c][d]]==2) no_type_2 =0;
            if(vt_structure_types[clusters_defect[c][d]]==4) no_type_4 =0;
            if(vt_structure_types[clusters_defect[c][d]]!=3) all_type_3=0;
            if(vt_structure_types[clusters_defect[c][d]]!=4) all_type_4=0;
            std::cout << clusters_defect[c][d] << '\t' << vt_structure_types[clusters_defect[c][d]] << '\n';
            
            //            cluster_index[clusters_defect[c][d]]=-c-1;
            //            cluster_sizes[clusters_defect[c][d]]=size;
        }
        if(size==6 && all_type_4==1) typeA_vacancy_count++;
        if(size==6 && no_type_2 ==1)
        {
            for(unsigned int d=0; d<size; d++)
                vt_structure_types[clusters_defect[c][d]]=4;
            typeB_vacancy_count++;
        }
        if(size==3 && all_type_4==1) typeC_vacancy_count++;
        if(size==2 && all_type_3 ==1) dislocation_count++;
        
        if(size>2  && no_type_4 ==1)
        {
            for(unsigned int d=0; d<size; d++)
                vt_structure_types[clusters_defect[c][d]]=2;
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
    for(unsigned int c=0; c<clusters_defect.size(); c++)
    {
        unsigned int size = clusters_defect[c].size();
        for(unsigned int d=0; d<size; d++)
        {
            cluster_index[clusters_defect[c][d]]=-c-1;
            cluster_sizes[clusters_defect[c][d]]=size;
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
    
    for(unsigned int c=0; c<clusters_defect.size(); c++)
    {
        if(clusters_defect[c].size()<max_cluster_size)
            c_sizes[clusters_defect[c].size()]++;
        sum_of_squares += clusters_defect[c].size()*clusters_defect[c].size();
    }
    
    std::cout << filename_data          << '\t';
    std::cout << number_of_particles    << '\t';
    std::cout << total_defect_particles << '\t';
    std::cout << clusters_defect.size()    << '\t';
    
    if(!clusters_defect.empty())
    {
        double average = double(total_defect_particles)/double(clusters_defect.size());
        double stdev   = sqrt(double(sum_of_squares)/double(clusters_defect.size()) - average*average);
        
        std::cout << (clusters_defect. back()).size() << '\t';  // SIZE OF SMALLEST DEFECT CLUSTER
        std::cout << (clusters_defect.front()).size() << '\t';  // SIZE OF LARGEST DEFECT CLUSTER
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



