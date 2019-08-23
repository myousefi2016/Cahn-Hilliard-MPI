#include "MpiSetup.h"

/**
Function to set up and initialize MPI for this run
This function initializes MPI on all processes and assigns a unique ID for
each MPI process in the simulation. In addition, the domain is partitioned
along X, Y and Z and a new communicator (CART_COMM) is created for the
purpose of exchanging information between processes.
For each MPI process, this also defines the IDs of neighboring processes
along X, Y and Z.
At the moment, the number of partitions along X, Y and Z are provided by the
user as command line arguments when he executes this code.
\verbatim
Example domain decomposition in 2D:      dims[0] = 4    (partitions along X)
                                         dims[1] = 3    (partitions along Y)
                                         dims[2] = 1    (partitions along Z)
 The grid below shows how MPI ranks are placed and the coordinates of each MPI rank
       X           Y           Z
       |           |           |
       |           |           |
  ( coords[0],  coords[1],  coords[2] )
 +--------------+--------------+--------------+--------------+
 |              |              |              |              | 
 |              |              |              |              | 
 |  (0,2,0)     |  (1,2,0)     |  (2,2,0)     |  (3,2,0)     | 
 |              |              |              |              |
 |  rank = 2    |  rank = 5    |  rank = 8    |  rank = 11   | 
 +--------------+--------------+--------------+--------------+
 |              |              |              |              | 
 |              |              |              |              |
 |  (0,1,0)     |  (1,1,0)     |  (2,1,0)     |  (3,1,0)     |       
 |              |              |              |              |
 |  rank = 1    |  rank = 4    |  rank = 7    |  rank = 10   | 
 +--------------+--------------+--------------+--------------+
 |              |              |              |              | 
 |              |              |              |              | 
 |  (0,0,0)     |  (1,0,0)     |  (2,0,0)     |  (3,0,0)     | 
 |              |              |              |              | 
 |  rank = 0    |  rank = 3    |  rank = 6    |  rank = 9    | 
 +--------------+--------------+--------------+--------------+
\endverbatim
*/
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
               int* nbr_TOP)           // pointer to --> ID of neighboring process to my top    (i,j,k+1)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);                     // number and value of command-line arguments
    MPI_Comm_size(MPI_COMM_WORLD, numprocs);    // get the total number of MPI processes
    MPI_Comm_rank(MPI_COMM_WORLD, myid);        // get my process ID

    // get and print processor name
    int len;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(name, &len);
    std::cout << "Processor name = " << name << std::endl;

    // boolean variables representing whether the domain is periodic along X, Y and Z
    int periods[3];
    periods[0] = 1;   // 0 = not periodic along X ... to make it periodic, use a value of 1
    periods[1] = 1;   // 0 = not periodic along Y ... to make it periodic, use a value of 1
    periods[2] = 1;   // 0 = not periodic along Z ... to make it periodic, use a value of 1

//  // automatically partition the 3D domain along X, Y and Z
//  // this is not being currently used because dims[0], dims[1] and dims[2] are explicitly specified by the user
//  dims[0] = 0;
//  dims[1] = 0;
//  dims[2] = 0;
//  MPI_Dims_create(*numprocs,3,dims);

    // user specified domain partitioning from the command line
    // if you wish to partition the domain automatically, comment out the next three lines and uncomment the block above
    dims[0] = atoi(argv[1]);  // convert character to integer - domain partitions along X
    dims[1] = atoi(argv[2]);  // convert character to integer - domain partitions along Y
    dims[2] = atoi(argv[3]);  // convert character to integer - domain partitions along Z

    // create a new communicator (CART_COMM) with Cartesian topology
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &CART_COMM);

    // get my position in the new communicator
    MPI_Comm_rank(CART_COMM,myid);

    // get process ids of my 6 neighbors along X, Y and Z
    MPI_Cart_shift(CART_COMM, 0, 1, nbr_WEST,   nbr_EAST);
    MPI_Cart_shift(CART_COMM, 1, 1, nbr_SOUTH,  nbr_NORTH);
    MPI_Cart_shift(CART_COMM, 2, 1, nbr_BOTTOM, nbr_TOP);

    // get the (integer) coordinates of this process in the Cartesian topology
    MPI_Cart_get(CART_COMM, 3,  dims, periods, coords);

    MPI_Barrier(CART_COMM);
    if(*myid==0) std::cout << "                                                   W  E  S  N  B  T" << std::endl;
    MPI_Barrier(CART_COMM);

    std::cout << "MPI rank " << *myid << " MPI coordinates: (" << coords[0] << ", " << coords[1] << ", " << coords[2] 
                             << ") - neighbors :" << *nbr_WEST << "  " << *nbr_EAST   << "  " 
                                                  << *nbr_SOUTH << "  " << *nbr_NORTH << "  " 
                                                  << *nbr_BOTTOM << "  " << *nbr_TOP << std::endl;
    MPI_Barrier(CART_COMM);
}
