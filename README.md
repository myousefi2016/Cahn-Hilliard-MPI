# Cahn Hilliard MPI (Phase-Field Simulation of Spinodal Decomposition)

Description: This C++/MPI program does phase-field simulation of spinodal decomposition based on this references:

1. [Cahn-Hilliard Equation](https://en.wikipedia.org/wiki/Cahn%E2%80%93Hilliard_equation)

Dependencies: OpenMPI 4.0.1, Intel Compiler 18.0, VTK 8.9, CMake 3.13.1

Make sure VTK 8.9 is built with MPI enabled.

It might be possible to compile this program with GCC compiler as well. Please create an issue if you encountered any error or unexpected behavior.

In order to compile the program use these commands in UNIX shell terminals:

```
git clone git@github.com:myousefi2016/Cahn-Hilliard-MPI.git
cd Cahn-Hilliard-MPI && mkdir build && cd build
cmake ..
make
```

To run the program after compilation just running this command:

```
mkdir out && mpirun -np #nprocs ./Cahn-Hilliard-MPI #decompX #decompY #decompZ
```

where `#nprocs` is the number of processors available in your machine. And, `#decompX`, `#decompY`, `#decompZ` are the domain decomposition configuration in X, Y, Z directions respectively to create MPI topology. Note that you need to make sure: `#nprocs = #decompX * #decompY * #decompZ`.


It will store the results in `out` directory as vtk files in parallel format. Good luck and if you use this piece of code for your research don't forget to give attribute to this github repository.

![alt text](https://raw.githubusercontent.com/myousefi2016/Cahn-Hilliard-MPI/master/animation/out.gif)
