# Minimum required Version
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# Projects name
project(time_propagator)

option(DEBUG "Enable debugging" OFF)

enable_language(Fortran)

set(SOURCE_TIME ${SOURCE_TIME}  source/main.F90
                                source/mod_types.F90
                                source/mod_constants.F90
                                source/mod_imRK4_parameters.F90
                                source/mod_math.F90
                                source/mod_matrices.F90
                                source/mod_hamiltonian.F90
                                source/mod_imRK4.F90
                                )

# Execulable names
set(TIME_EXE tp.exe)
set(CMAKE_Fortran_COMPILER ifort)

add_executable(${TIME_EXE} ${SOURCE_TIME})

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

target_link_libraries(${TIME_EXE} -mkl)
target_compile_options(${TIME_EXE} PRIVATE -no-wrap-margin
                                            -fpp
                                            -DVERSION="${GIT_VERSION}"
                                            -mkl
                                            -xHost)


if(DEBUG)
  target_compile_options(${TIME_EXE} PRIVATE  -C
                                              -CB
                                              -check uninit
                                              -debug all
                                              -warn all
                                              -ftrapuv
                                              -traceback
                                              -fpp
                                              -g
                                              -O0)
else()
  target_compile_options(${TIME_EXE} PRIVATE -O3)
endif(DEBUG)
