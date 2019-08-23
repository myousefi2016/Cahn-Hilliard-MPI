#include "Initialization.h"

void Initialization(double *c, int MX, int MY, int MZ, int nn, int myid)
{

    int MXP = nn+MX+nn;
    int MYP = nn+MY+nn;
    int MZP = nn+MZ+nn;

    srand(time(NULL) + myid);

    for (unsigned int idz = 0; idz < MZ; idz++) {
        int IDz = idz + nn;
     for (unsigned int idy = 0; idy < MY; idy++) {
                int IDy = idy + nn;
      for (unsigned int idx = 0; idx < MX; idx++) {
                        int IDx = idx + nn;

                                int IDflattened = IDx+IDy*MXP+IDz*MXP*MYP;

                                double f = (double)rand() / RAND_MAX;
      				c[IDflattened] = -1.0 + 2.0*f;

    }
     }
      }

}
