#include "DomainDecomposition.h"

// this function returns the beginning and ending node numbers for this
// MPI rank along the specified direction

node_range DomainDecomposition1D( int n, int proc_dim, int proc_coord)
{
    int nlocal;
    int deficit;
    node_range limits;

    nlocal = n / proc_dim;
    limits.beg = proc_coord* nlocal + 1;
    deficit = n%proc_dim;
    limits.beg = limits.beg + std::min(proc_coord,deficit);

    if (proc_coord < deficit) nlocal++;

    limits.end = limits.beg + nlocal - 1;

    if ((limits.end > n) || (proc_coord == proc_dim-1)) limits.end = n;

    // C++ begins counting from 0
    limits.beg--;
    limits.end--;

    return limits;
}

/** 
This function calculates the size of the (local) buffer in this process after 3D domain partitioning
and also calculates other parameters related to the 3D partitioning that are used for parallel meshing
within the sub-domain.
Example domain decomposition in 3D:
  dims[0] = 4    (partitions along X)
  dims[1] = 3    (partitions along Y)
  dims[2] = 1    (partitions along Z)
Suppose the global node size of the RAW data is 11 x 6 x 1. Then this routine splits these nodes
among the MPI ranks in each direction (approximately equal divisions)
\verbatim
So for say, rank = 7:    x_range.beg = 6, x_range.end = 8
                         y_range.beg = 2, y_range.end = 3
                         z_range.beg = 0, z_range.end = 0
                         MX = 3, MY = 2, MZ = 1
                         local_origin_x = 6
                         local_origin_y = 2
                         local_origin_z = 0
 ---+--------------+--------------+--------------+-----------+
  5 |              |              |              |           | 
    |  rank = 2    |  rank = 5    |  rank = 8    | rank = 11 | 
  4 |              |              |              |           | 
 ---+--------------+--------------+--------------+-----------+
  3 |              |              |              |           | 
    |  rank = 1    |  rank = 4    |  rank = 7    | rank = 10 |
  2 |              |              |              |           |
 ---+--------------+--------------+--------------+-----------+
  1 |              |              |              |           | 
    |  rank = 0    |  rank = 3    |  rank = 6    | rank = 9  | 
  0 |              |              |              |           | 
 ---+--------------+--------------+--------------+-----------+
    |  0   1   2   |  3   4   5   |  6   7   8   |  9   10   | 
\endverbatim
*/
void DomainDecomposition3D(// inputs
                    int      & myid,           // MPI rank
                    MPI_Comm & CART_COMM,      // MPI communicator name
                    int      * dims,           // number of partitions of the domain along X, Y and Z, dims[0], dims[1] and dims[2]
                    int      * coords,         // (X,Y,Z) "coordinates" of this partition
                    int      & nodes_x, 
                    int      & nodes_y, 
                    int      & nodes_z,
                    double   & delta,
                    double   & x_min,
                    double   & y_min,
                    double   & z_min,
                    // outputs
                    node_range & x_range, 
                    node_range & y_range, 
                    node_range & z_range,
                    double & local_origin_x, 
                    double & local_origin_y, 
                    double & local_origin_z,
                    int & LX,
                    int & LY,
                    int & LZ)
{
    // user decides how the domain is to be partitioned based on some knowledge about the geometry
    MPI_Barrier(CART_COMM);
    if(myid==0) std::cout << std::endl;
    MPI_Barrier(CART_COMM);
    if(myid==0) std::cout << "Domain partitioning (along X, Y and Z): " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
    MPI_Barrier(CART_COMM);
    if(myid==0) std::cout << std::endl;
    MPI_Barrier(CART_COMM);

    // Decompose processors to get corresponding node ranges along X Y and Z
    x_range = DomainDecomposition1D( nodes_x, dims[0], coords[0]);
    y_range = DomainDecomposition1D( nodes_y, dims[1], coords[1]);
    z_range = DomainDecomposition1D( nodes_z, dims[2], coords[2]);

    // calculate the global (X,Y,Z) coordinate of the "local" origin for this particular task (myid)
    local_origin_x = x_range.beg;
    local_origin_y = y_range.beg;
    local_origin_z = z_range.beg;

    // dimensions for the (local) buffers inside each process
    LX = x_range.end - x_range.beg + 1;
    LY = y_range.end - y_range.beg + 1;
    LZ = z_range.end - z_range.beg + 1;

    std::cout << "Rank " << myid << " X: " << x_range.beg << " " << x_range.end << " LX = " << LX
                                 << " Y: " << y_range.beg << " " << y_range.end << " LY = " << LY
                                 << " Z: " << z_range.beg << " " << z_range.end << " LZ = " << LZ << std::endl;
}
