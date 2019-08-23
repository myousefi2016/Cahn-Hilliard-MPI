#include <vector>
#include <string>
#include <iostream>
#include "VtkParallelWriter.h"

// function to operate on the point attribute data
void execute (void* arg) {
  Args* args = reinterpret_cast<Args*>(arg);
  auto info = args->pf->GetOutputInformation(0);
  auto output_tmp = args->pf->GetOutput();
  auto input_tmp  = args->pf->GetInput();
  vtkStructuredGrid* output = dynamic_cast<vtkStructuredGrid*>(output_tmp);
  vtkStructuredGrid* input  = dynamic_cast<vtkStructuredGrid*>(input_tmp);
  output->ShallowCopy(input);
  output->SetExtent(args->local_extent);
}

void vtkParallelWriter(int argc, char *argv[],std::vector<double*>colors, std::vector<char*> names, int LX, int LY, int LZ, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, double local_origin_x, double local_origin_y, double local_origin_z, int nn, int timesnapshot) {
  int global_extent[6] = {x_min, x_max, y_min, y_max, z_min, z_max};
  bool flagX, flagY, flagZ;
  int new_local_origin_x;
  int new_local_origin_y;
  int new_local_origin_z;
  flagX = false;
  flagY = false;
  flagZ = false;
  if (local_origin_x == 0) {
	flagX = true;
  }
  if (local_origin_y == 0) {
        flagY = true;
  }
  if (local_origin_z == 0) {
        flagZ = true;
  }
  if (flagX == false) {
    new_local_origin_x = local_origin_x - 1;
  }
  if (flagY == false) {
    new_local_origin_y = local_origin_y - 1;
  }
  if (flagZ == false) {
    new_local_origin_z = local_origin_z - 1;
  }
  if (flagX == true) {
    new_local_origin_x = local_origin_x;
  }
  if (flagY == true) {
    new_local_origin_y = local_origin_y;
  }
  if (flagZ == true) {
    new_local_origin_z = local_origin_z;
  }
  int local_extent[6] = {new_local_origin_x, local_origin_x+LX-1, new_local_origin_y, local_origin_y+LY-1, new_local_origin_z, local_origin_z+LZ-1};
  int dims[3] = {local_origin_x+LX-new_local_origin_x, local_origin_y+LY-new_local_origin_y, local_origin_z+LZ-new_local_origin_z};

  // Create and Initialize vtkMPIController
  auto contr = vtkSmartPointer<vtkMPIController>::New();
  if (timesnapshot == 0) {
  contr->Initialize(&argc, &argv, 1);
  }
  int nranks = contr->GetNumberOfProcesses();
  int rank   = contr->GetLocalProcessId();

  // Create grid points, allocate memory and Insert them
  auto points = vtkSmartPointer<vtkPoints>::New();
  for (int k=0; k<dims[2]; ++k) {
    for (int j=0; j<dims[1]; ++j) {
      for (int i=0; i<dims[0]; ++i) {
        points->InsertNextPoint(new_local_origin_x+i, new_local_origin_y+j, new_local_origin_z+k);

   }
   }
   }

  // Create a density field. Note that the number of cells is always less than
  // number of grid points by an amount of one so we use dims[i]-1
  std::vector<vtkSmartPointer<vtkDoubleArray>> vtkColors;
  for (int iterator = 0; iterator < colors.size(); iterator++) {
  auto tempColor = vtkSmartPointer<vtkDoubleArray>::New();
  tempColor->SetNumberOfComponents(1);
  tempColor->SetName(names[iterator]);
  int LXP = nn+LX+nn;
  int LYP = nn+LY+nn;
  int LZP = nn+LZ+nn;
  for (int k=0; k<dims[2]; ++k) {
        int K;
        if (flagZ == true) {
	K = k + nn;
        }
        if (flagZ == false) {
        K = k;
        }
    for (int j=0; j<dims[1]; ++j) {
                int J;
                if (flagY == true) {
		J = j + nn;
                }
                if (flagY == false) {
                J = j;
                }
      for (int i=0; i<dims[0]; ++i) {
			int I;
			if (flagX == true) {
			I = i + nn;
			}
			if (flagX == false) {
			I = i;
			}
				int IDflattened = I+J*LXP+K*LXP*LYP;	
        tempColor->InsertNextTuple1(colors[iterator][IDflattened]);
      }
      }
      }
  vtkColors.push_back(tempColor);
  }

  // Create a vtkProgrammableFilter
  auto pf = vtkSmartPointer<vtkProgrammableFilter>::New();

  // Initialize an instance of Args
  Args args;
  args.pf = pf;
  for(int i=0; i<6; ++i) args.local_extent[i] = local_extent[i];

  pf->SetExecuteMethod(execute, &args);

  // Create a structured grid and assign point data and cell data to it
  auto structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
  structuredGrid->SetExtent(global_extent);
  pf->SetInputData(structuredGrid);
  structuredGrid->SetPoints(points);
  for (int iterator = 0; iterator < vtkColors.size(); iterator++) {
  structuredGrid->GetPointData()->AddArray(vtkColors[iterator]);
  }

  std::string fileName = std::string("./out/output_") + std::to_string(timesnapshot) + ".pvts";

  // Create the parallel writer and call some functions
  auto parallel_writer = vtkSmartPointer<vtkXMLPStructuredGridWriter>::New();
  parallel_writer->SetInputConnection(pf->GetOutputPort());
  parallel_writer->SetController(contr);
  parallel_writer->SetFileName(fileName.c_str());
  parallel_writer->SetNumberOfPieces(nranks);
  parallel_writer->SetStartPiece(rank);
  parallel_writer->SetEndPiece(rank);
  parallel_writer->SetDataModeToBinary();
  parallel_writer->Update();
  parallel_writer->Write();

}
