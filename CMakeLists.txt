# CMakeList.txt : CMake project for PQP_signed, include source and define
# project specific logic here.3.8
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project (PQP_signed)

# Get Eigen3
find_package(Eigen3 REQUIRED)
if (Eigen3_FOUND)
    message(STATUS "${EIGEN3_VERSION_STRING}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif ()

# Get CGAL
find_package(CGAL REQUIRED)
if (CGAL_FOUND)
    include(${CGAL_USE_FILE})
else ()
    message("ERROR: this program requires CGAL and will not be compiled.")
endif ()

# Add source to this project's executable.
file(GLOB Model3D Model3D/*.h Model3D/*.cpp)
file(GLOB PQP PQP/*.h PQP/*.cpp)
file(GLOB mc mc/*.h mc/*.cpp)
add_executable(${PROJECT_NAME} ${Model3D} ${PQP} ${mc})


# TODO: Add tests and install targets if needed.
include_directories(Model3D/ PQP/ mc/)

