cmake_minimum_required(VERSION 3.7)

project(spoic LANGUAGES C Fortran)

add_executable(spoic
spoic.f90
spec)
target_link_libraries(spoic spec)

add_subdirectory("SPEC/")
