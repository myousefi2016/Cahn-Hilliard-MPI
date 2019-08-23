#ifndef MPISETUP_H
#define MPISETUP_H

      #include <iostream>
      #include <cstdlib>
      #include <mpi.h>
      #include <sstream>
      #include <iomanip>

void MpiSetup (int argc,               // argc is short for argument count (set depending on how many arguments the user enters at the command line) 
               char *argv[],           // pointer to a char array  ... argv is short for argument values ... array size set based on the value of argc
               int* numprocs,          // pointer to an integer - number of distinct MPI processes on which this code will be executed 
               int* myid,              // pointer to an integer - process ID
               const int ndims,        // number of spatial dimensions for domain partitioning
               int* dims,              // pointer to --> dims[0] - number of partitions of the domain along X, Y and Z...dims[0], dims[1], dims[2]
               int* coords,            // pointer to --> coords[0] - coordinates of this process within the Cartesian topology ...coords[0], coords[1], coords[2]
               MPI_Comm & CART_COMM,   // name of the Cartesian communicator
               int* nbr_WEST,          // pointer to --> ID of neighboring process to my west   (i-1,j,k)
               int* nbr_EAST,          // pointer to --> ID of neighboring process to my east   (i+1,j,k)
               int* nbr_SOUTH,         // pointer to --> ID of neighboring process to my south  (i,j-1,k)
               int* nbr_NORTH,         // pointer to --> ID of neighboring process to my north  (i,j+1,k)
               int* nbr_BOTTOM,        // pointer to --> ID of neighboring process to my bottom (i,j,k-1)
               int* nbr_TOP);           // pointer to --> ID of neighboring process to my top    (i,j,k+1)


#endif
