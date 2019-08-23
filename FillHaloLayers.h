#ifndef FILLHALOLAYERS_H
#define FILLHALOLAYERS_H

#include <iostream>
#include <mpi.h>      // MPI header files
#include "ExchangeDBL.h"

void FillHaloLayers(const int       nn,              // ghost layer thickness
                           const int       LX,              // number of nodes along X (local for this MPI process)
                           const int       LY,              // number of nodes along Y (local for this MPI process)
                           const int       LZ,              // number of nodes along Z (local for this MPI process)
                           const int       myid,            // MPI process id or rank
                           const MPI_Comm  CART_COMM,       // Cartesian communicator
                           const int       nbr_WEST,        // neighboring MPI process to my west
                           const int       nbr_EAST,        // neighboring MPI process to my east
                           const int       nbr_SOUTH,       // neighboring MPI process to my south
                           const int       nbr_NORTH,       // neighboring MPI process to my north
                           const int       nbr_BOTTOM,      // neighboring MPI process to my bottom
                           const int       nbr_TOP,         // neighboring MPI process to my top
                                 double    *color);              // color

#endif
