double Laplacian(double *c, double dx, double dy, double dz, int x, int y, int z, int MX, int MY, int MZ, int nn);

void chemicalPotential(double *c, double *mu, double dx, double dy, double dz, double kappa, double e_AA, double e_BB, double e_AB, int MX, int MY, int MZ, int nn);

void cahnHilliard(double *cnew, double *mu, double D, double dt, double dx, double dy, double dz, int MX, int MY, int MZ, int nn);
