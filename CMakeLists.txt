cmake_minimum_required(VERSION 3.7)

project(spoic LANGUAGES C Fortran)

add_compile_options(-fbacktrace -g)

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
message(STATUS "HDF5 libraries found: ${HDF5_FOUND}")
message(STATUS "HDF5 C include directors : ${HDF5_C_INCLUDE_DIRS}")
message(STATUS "HDF5 Fortran include directors : ${HDF5_Fortran_INCLUDE_DIRS}")
message(STATUS "HDF5 Fortran libraries : ${HDF5_Fortran_LIBRARIES}")

add_executable(spoic.x
spoic.f90)
target_link_libraries(spoic.x PUBLIC
    spoic ${HDF5_Fortran_LIBRARIES} ${HDF5_C_LIBRARIES})

add_library(spoic
read_spec.f90
)
target_link_libraries(spoic PUBLIC spec ${HDF5_Fortran_LIBRARIES} ${HDF5_C_LIBRARIES})
target_include_directories(spoic PUBLIC ${HDF5_C_INCLUDE_DIRS} ${HDF5_Fortran_INCLUDE_DIRS})

add_subdirectory("SPEC/")
