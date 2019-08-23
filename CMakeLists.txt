CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(Cahn-Hilliard-MPI)

find_package(MPI)

find_package(VTK REQUIRED)

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

include_directories(SYSTEM ${MPI_INCLUDE_PATH} ${VTK_USE_FILE})
 
FILE(GLOB SOURCES "*.cpp" "*.cc" "*.c" "*.h")
ADD_EXECUTABLE(Cahn-Hilliard-MPI ${SOURCES})

if(VTK_LIBRARIES)
  target_link_libraries(Cahn-Hilliard-MPI ${VTK_LIBRARIES} ${MPI_LIBRARIES})
else()
  target_link_libraries(Cahn-Hilliard-MPI vtkHybrid vtkWidgets ${MPI_LIBRARIES})
endif()
