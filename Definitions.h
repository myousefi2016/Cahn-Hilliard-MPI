#ifndef DEFINITIONS_H
#define DEFINITIONS_H

//    C++ headers

      #include <iostream>     // cout()
      #include <cmath>        // pow()
      #include <ctime>        // clock_t, clock(), CLOCKS_PER_SEC
      #include <mpi.h>        // MPI 

      #include "DomainDecomposition.h"

//    data structures

//    MPI 

      int numprocs;          // total number of processors
      int myid;              // processor id
      const int ndims = 3;   // number of dimensions of the Cartesian space
      int dims[ndims];       // number of domain partitions along X, Y and Z
      MPI_Comm CART_COMM;    // new Cartesian communicator after automatic 3D domain decomposition
      int coords[ndims];     // 3D coordinates of process "myid" after domain decomposition
      int nbr_WEST;          // id of neighbor in location (i-1)
      int nbr_EAST;          // id of neighbor in location (i+1)
      int nbr_SOUTH;         // id of neighbor in location (j-1)
      int nbr_NORTH;         // id of neighbor in location (j+1)
      int nbr_BOTTOM;        // id of neighbor in location (k-1)
      int nbr_TOP;           // id of neighbor in location (k+1)

//    parameters

// set a 3D volume
// To compile it with nvcc execute: nvcc -O2 -o set3d set3d.cu
//define the data set size (cubic volume)
      int NX = 128;
      int NY = 128;
      int NZ = 128;

      double delta = 1.0;
      double dt = 0.01;
      double e_AA = -(2.0/9.0);
      double e_BB = -(2.0/9.0);
      double e_AB = (2.0/9.0);
      int t_f = 25000;
      int t_freq = 10;
      double kappa = 0.5;
      double D = 1.0;

      double x_min = 0;    // global minimum X coordinate
      double x_max = NX-1;
      double y_min = 0;    // global minimum Y coordinate
      double y_max = NY-1;
      double z_min = 0;    // global minimum Z coordinate
      double z_max = NZ-1;

//    local parameters (obtained after MPI domain decomposition)

      int LX;   // number of lattice points along X
      int LY;   // number of lattice points along Y
      int LZ;   // number of lattice points along Z

      double local_origin_x;  // X coordinate of the local domain origin
      double local_origin_y;  // Y coordinate of the local domain origin
      double local_origin_z;  // Z coordinate of the local domain origin

#endif
