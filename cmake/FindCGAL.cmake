# - Try to find the LIBIGL library
# Once done this will define
#
#  LIBIGL_FOUND - system has LIBIGL
#  LIBIGL_INCLUDE_DIR - **the** LIBIGL include directory
if(CGAL_FOUND)
    return()
endif()

find_path(
    PATHS
        F:/xinCode/GraghicLibrary/cgal
)

include(FindPackageHandleStandardArgs)
