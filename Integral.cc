#include "Integral.h"

double Integral(double *c, int MX, int MY, int MZ, int nn)
{
    double summation = 0.0;    

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

			        summation = summation + c[IDflattened];

    }
     }
      }

    return summation;

}
