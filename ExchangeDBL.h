#ifndef EXCHANGEDBL_H
#define EXCHANGEDBL_H

#include <iostream>
#include <mpi.h>      // MPI header files

void ExchangeDBL  (const int      & nn,                // number of ghost cell layers
                   const int      & MX,                // number of voxels along X in this process
                   const int      & MY,                // number of voxels along Y in this process
                   const int      & MZ,                // number of voxels along Z in this process
                   const int      & myid,              // my process id
                   const MPI_Comm & CART_COMM,         // Cartesian topology communicator
                   const int      & nbr_WEST,          // process id of my western neighbor
                   const int      & nbr_EAST,          // process id of my eastern neighbor
                   const int      & nbr_SOUTH,         // process id of my southern neighbor
                   const int      & nbr_NORTH,         // process id of my northern neighbor
                   const int      & nbr_BOTTOM,        // process id of my bottom neighbor
                   const int      & nbr_TOP,           // process id of my top neighbor
                      double      * color);             // pointer to the 3D array being exchanged (of type double)

#endif
