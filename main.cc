#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>

#include "FillHaloLayers.h"
#include "Initialization.h"
#include "Integral.h"
#include "Kernels.h"
#include "MpiSetup.h"
#include "DomainDecomposition.h"
#include "Definitions.h"
#include "VtkParallelWriter.h"

using namespace std;

int main(int argc, char* argv[])
{
        node_range x_range;
        node_range y_range;
        node_range z_range;

	MpiSetup(argc, argv, &numprocs, &myid, ndims,
                 &dims[0], &coords[0], CART_COMM,
                 &nbr_WEST, &nbr_EAST,
                 &nbr_SOUTH, &nbr_NORTH,
                 &nbr_BOTTOM, &nbr_TOP);	

	DomainDecomposition3D(myid, CART_COMM, dims, coords,
                       	      NX, NY, NZ,
                              delta,
                              x_min, y_min, z_min,
                              x_range, y_range, z_range,
                              local_origin_x, local_origin_y, local_origin_z,
                              LX,                // local nodes along X
                              LY,                // local nodes along Y
                              LZ);               // local nodes along Z

	const int nn = 1;

	const int size1 = (nn+LX+nn) * (nn+LY+nn) * (nn+LZ+nn);

	double *c;
	double *mu;

	if ((c = (double *)malloc((size1)*sizeof(double))) == 0) {fprintf(stderr,"malloc1 Fail \n"); return 1;}
	if ((mu = (double *)malloc((size1)*sizeof(double))) == 0) {fprintf(stderr,"malloc1 Fail \n"); return 1;}

	Initialization(c,LX,LY,LZ,nn,myid);

	MPI_Barrier(CART_COMM);	

	FillHaloLayers(nn,              // ghost layer thickness
                       LX,              // number of nodes along X (local for this MPI process)
                       LY,              // number of nodes along Y (local for this MPI process)
                       LZ,              // number of nodes along Z (local for this MPI process)
                       myid,            // MPI process id or rank
                       CART_COMM,       // Cartesian communicator
                       nbr_WEST,        // neighboring MPI process to my west
                       nbr_EAST,        // neighboring MPI process to my east
                       nbr_SOUTH,       // neighboring MPI process to my south
                       nbr_NORTH,       // neighboring MPI process to my north
                       nbr_BOTTOM,      // neighboring MPI process to my bottom
                       nbr_TOP,         // neighboring MPI process to my top
                       c);

	MPI_Barrier(CART_COMM);

    	int t = 0;

        string name_c = "./out/integral_c.txt";
        ofstream ofile_c (name_c.c_str());

        string name_mu = "./out/integral_mu.txt";
        ofstream ofile_mu (name_mu.c_str());

        double localIntegral_c;

	double globalIntegral_c;

        double localIntegral_mu;

        double globalIntegral_mu;

	localIntegral_c = Integral(c,LX,LY,LZ,nn);

	MPI_Reduce(&localIntegral_c,&globalIntegral_c,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	if (myid == 0) {

	ofile_c << 0 << "," << globalIntegral_c << std::endl;

        }

	std::vector<char*> names;
	names.push_back("c");

	std::vector<double*>colors_initial;
        colors_initial.push_back(c);

    	while (t < t_f) {

		if (myid == 0) {
    		printf("Timestep is: %d\n",t);	
  		}

	 	chemicalPotential(c,mu,delta,delta,delta,kappa,e_AA,e_BB,e_AB,LX,LY,LZ,nn);

		MPI_Barrier(CART_COMM);

        	FillHaloLayers(nn,              // ghost layer thickness
                 	      LX,              // number of nodes along X (local for this MPI process)
                       	      LY,              // number of nodes along Y (local for this MPI process)
                              LZ,              // number of nodes along Z (local for this MPI process)
                              myid,            // MPI process id or rank
                              CART_COMM,       // Cartesian communicator
                              nbr_WEST,        // neighboring MPI process to my west
                              nbr_EAST,        // neighboring MPI process to my east
                              nbr_SOUTH,       // neighboring MPI process to my south
                              nbr_NORTH,       // neighboring MPI process to my north
                              nbr_BOTTOM,      // neighboring MPI process to my bottom
                              nbr_TOP,         // neighboring MPI process to my top
                              mu);

		MPI_Barrier(CART_COMM);

    		cahnHilliard(c,mu,D,dt,delta,delta,delta,LX,LY,LZ,nn);

		MPI_Barrier(CART_COMM);

                FillHaloLayers(nn,              // ghost layer thickness
                              LX,              // number of nodes along X (local for this MPI process)
                              LY,              // number of nodes along Y (local for this MPI process)
                              LZ,              // number of nodes along Z (local for this MPI process)
                              myid,            // MPI process id or rank
                              CART_COMM,       // Cartesian communicator
                              nbr_WEST,        // neighboring MPI process to my west
                              nbr_EAST,        // neighboring MPI process to my east
                              nbr_SOUTH,       // neighboring MPI process to my south
                              nbr_NORTH,       // neighboring MPI process to my north
                              nbr_BOTTOM,      // neighboring MPI process to my bottom
                              nbr_TOP,         // neighboring MPI process to my top
                              c);
		MPI_Barrier(CART_COMM);

	if (t == 0) {

        localIntegral_mu = Integral(mu,LX,LY,LZ,nn);

        MPI_Reduce(&localIntegral_mu,&globalIntegral_mu,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

        if (myid == 0) {

        ofile_mu << t << "," << globalIntegral_mu << std::endl;

        }

	names.push_back("mu");
	
	colors_initial.push_back(mu);

	vtkParallelWriter(argc, argv,colors_initial,names,LX,LY,LZ,x_min,x_max,y_min,y_max,z_min,z_max,local_origin_x,local_origin_y,local_origin_z,nn,t);

        }

	if (t % t_freq == 0 && t > 0) {		

        	localIntegral_c = Integral(c,LX,LY,LZ,nn);

        	MPI_Reduce(&localIntegral_c,&globalIntegral_c,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	        if (myid == 0) {

        	ofile_c << t << "," << globalIntegral_c << std::endl;

        	}

        	localIntegral_mu = Integral(mu,LX,LY,LZ,nn);

        	MPI_Reduce(&localIntegral_mu,&globalIntegral_mu,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

        	if (myid == 0) {
 
        	ofile_mu << t << "," << globalIntegral_mu << std::endl;
 
        	} 

	        std::vector<double*>colors;
		colors.push_back(c);
		colors.push_back(mu);

		vtkParallelWriter(argc, argv,colors,names,LX,LY,LZ,x_min,x_max,y_min,y_max,z_min,z_max,local_origin_x,local_origin_y,local_origin_z,nn,t);		

	}


		t++;

	}

	free(c);
	free(mu);
	
    	MPI_Finalize();

	return 0;
}
