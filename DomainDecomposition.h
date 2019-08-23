#ifndef DOMAINDECOMPOSITION_H
#define DOMAINDECOMPOSITION_H

// include files

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <mpi.h>

//    define a struct to store the beginning and ending node numbers inside a process
      struct node_range
      {
        int beg;
        int end;
      };

node_range DomainDecomposition1D( int n, int proc_dim, int proc_coord);

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
                    int & LZ);

#endif
