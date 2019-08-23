// MPI Library
#include <mpi.h>

//VTK Library
#include <vtkXMLPStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkMPIController.h>
#include <vtkProgrammableFilter.h>
#include <vtkInformation.h>

struct Args {
  vtkProgrammableFilter* pf;
  int local_extent[6];
};

void execute (void* arg);
void vtkParallelWriter(int argc, char *argv[],std::vector<double*>colors, std::vector<char*> names, int LX, int LY, int LZ, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, double local_origin_x, double local_origin_y, double local_origin_z, int nn, int timesnapshot);
