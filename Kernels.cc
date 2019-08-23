#include "Kernels.h"

double Laplacian(double *c, double dx, double dy, double dz, int x, int y, int z, int MX, int MY, int MZ, int nn)
{

  int xp, xn, yp, yn, zp, zn;

  int MXP = nn+MX+nn;     
  int MYP = nn+MY+nn;     

  xp = x+1;
  xn = x-1;
  yp = y+1;
  yn = y-1;
  zp = z+1;
  zn = z-1;

  double cxx = (c[xp+y*MXP+z*MXP*MYP] + c[xn+y*MXP+z*MXP*MYP] - 2.0*c[x+y*MXP+z*MXP*MYP]) / (dx*dx);
  double cyy = (c[x+yp*MXP+z*MXP*MYP] + c[x+yn*MXP+z*MXP*MYP] - 2.0*c[x+y*MXP+z*MXP*MYP]) / (dy*dy);
  double czz = (c[x+y*MXP+zp*MXP*MYP] + c[x+y*MXP+zn*MXP*MYP] - 2.0*c[x+y*MXP+z*MXP*MYP]) / (dz*dz);

  double result = cxx + cyy + czz;

  return result;

}

void chemicalPotential(double *c, double *mu, double dx, double dy, double dz, double kappa, double e_AA, double e_BB, double e_AB, int MX, int MY, int MZ, int nn)
{

    int MXP = nn+MX+nn;
    int MYP = nn+MY+nn;
    int MZP = nn+MZ+nn;

    for (unsigned int idz = 0; idz < MZ; idz++) {
	int IDz = idz + nn;
     for (unsigned int idy = 0; idy < MY; idy++) {
		int IDy = idy + nn;
      for (unsigned int idx = 0; idx < MX; idx++) {
			int IDx = idx + nn;

				int IDflattened = IDx+IDy*MXP+IDz*MXP*MYP;

			        mu[IDflattened] = ( 9.0 / 2.0 )*( ( c[IDflattened] + 1.0 ) * e_AA + ( c[IDflattened] - 1 ) * e_BB - 2.0 * c[IDflattened] * e_AB ) + 3.0 * c[IDflattened] + c[IDflattened] * c[IDflattened] * c[IDflattened] - kappa * Laplacian(c,dx,dy,dz,IDx,IDy,IDz,MX,MY,MZ,nn);

    }
     }
      }

}

void cahnHilliard(double *cnew, double *mu, double D, double dt, double dx, double dy, double dz, int MX, int MY, int MZ, int nn)
{

    int MXP = nn+MX+nn;
    int MYP = nn+MY+nn;
    int MZP = nn+MZ+nn;

    for (unsigned int idz = 0; idz < MZ; idz++) {
        int IDz = idz + nn;
     for (unsigned int idy = 0; idy < MY; idy++) {
                int IDy = idy + nn;
      for (unsigned int idx = 0; idx < MX; idx++) {
                        int IDx = idx + nn;

                                int IDflattened = IDx+IDy*MXP+IDz*MXP*MYP;

                                cnew[IDflattened] = cnew[IDflattened] + dt * D * Laplacian(mu,dx,dy,dz,IDx,IDy,IDz,MX,MY,MZ,nn);

    }
     }
      }

}
