// voro_faces: dump the face structure of every Voronoi cell (voro++ dev branch).
// Input  (stdin or file): line 1 "N Lx Ly Lz px py pz" (box [0,L), periodic flags 0/1), or, for a fully periodic
//                         triclinic box, "N T bx bxy by bxz byz bz" (lattice vectors a = (bx,0,0), b = (bxy,by,0),
//                         c = (bxz,byz,bz), the LAMMPS restricted-triclinic convention);
//                         then N lines "x y z" (particle i = line order, 0-based).
// Output (stdout): per cell, one line
//   i nfaces  { neighbour nverts v0 v1 ... }*nfaces
// with vertex indices local to the cell, each face's cycle as returned by voro++.
// Negative neighbours are container walls.
#include "voro++.hh"
#include <cstdio>
#include <vector>
#include <omp.h>
using namespace voro;

template <class C>
static void cells(C &con, int n, std::vector<std::vector<int>> &out_nb, std::vector<std::vector<int>> &out_fv,
                  std::vector<char> &ok, int threads) {
    #pragma omp parallel num_threads(threads)
    {
        voronoicell_neighbor_3d c(con);
        #pragma omp for schedule(dynamic, 256)
        for (auto cli = con.begin(); cli < con.end(); ++cli)
            if (con.compute_cell(c, cli)) {
                int pid = con.id[cli->ijk][cli->q];
                c.neighbors(out_nb[pid]);
                c.face_vertices(out_fv[pid]);
                ok[pid] = 1;
            }
    }
}

int main(int argc, char **argv) {
    FILE *fp = argc > 1 ? fopen(argv[1], "r") : stdin;
    char head[512];
    if (!fgets(head, sizeof head, fp)) return 1;
    int n, px = 1, py = 1, pz = 1; double Lx, Ly, Lz, bxy = 0, bxz = 0, byz = 0;
    char tag[8] = "";
    bool tri = sscanf(head, "%d %7s", &n, tag) == 2 && tag[0] == 'T';
    if (tri) {
        if (sscanf(head, "%d %*s %lf %lf %lf %lf %lf %lf", &n, &Lx, &bxy, &Ly, &bxz, &byz, &Lz) != 7) return 1;
    } else if (sscanf(head, "%d %lf %lf %lf %d %d %d", &n, &Lx, &Ly, &Lz, &px, &py, &pz) != 7) return 1;
    std::vector<double> x(3 * n);
    for (int i = 0; i < n; i++)
        if (fscanf(fp, "%lf %lf %lf", &x[3*i], &x[3*i+1], &x[3*i+2]) != 3) return 1;
    int threads = omp_get_max_threads();
    double dens = n / (Lx * Ly * Lz);
    int nx = std::max(1, int(Lx * pow(dens / 5, 1. / 3))), ny = std::max(1, int(Ly * pow(dens / 5, 1. / 3))),
        nz = std::max(1, int(Lz * pow(dens / 5, 1. / 3)));
    std::vector<std::vector<int>> out_nb(n), out_fv(n);
    std::vector<char> ok(n, 0);
    if (tri) {
        container_triclinic con(Lx, bxy, Ly, bxz, byz, Lz, nx, ny, nz, 8, threads);
        for (int i = 0; i < n; i++) con.put(i, x[3*i], x[3*i+1], x[3*i+2]);
        cells(con, n, out_nb, out_fv, ok, threads);
    } else {
        container_3d con(0, Lx, 0, Ly, 0, Lz, nx, ny, nz, px, py, pz, 8, threads);
        for (int i = 0; i < n; i++) con.put(i, x[3*i], x[3*i+1], x[3*i+2]);
        cells(con, n, out_nb, out_fv, ok, threads);
    }
    for (int i = 0; i < n; i++) {
        if (!ok[i]) { printf("%d 0\n", i); continue; }
        printf("%d %zu", i, out_nb[i].size());
        const std::vector<int> &fv = out_fv[i];
        size_t k = 0;
        for (size_t f = 0; f < out_nb[i].size(); f++) {
            int m = fv[k];
            printf(" %d %d", out_nb[i][f], m);
            for (int t = 1; t <= m; t++) printf(" %d", fv[k + t]);
            k += m + 1;
        }
        printf("\n");
    }
    return 0;
}
