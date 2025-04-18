# This file is part of monofonIC (MUSIC2)
# A software package to generate ICs for cosmological simulations
# Copyright (C) 2020 by Oliver Hahn
# 
# monofonIC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# monofonIC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.11)
set(PRGNAME monofonIC)
project(monofonIC C CXX)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -fno-omit-frame-pointer -g  -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -Wall -pedantic" CACHE STRING "Flags used by the compiler during Release builds." FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -fno-omit-frame-pointer -Wall -pedantic" CACHE STRING "Flags used by the compiler during RelWithDebInfo builds." FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O1 -march=native -DDEBUG -fno-omit-frame-pointer -Wall -pedantic" CACHE STRING "Flags used by the compiler during Debug builds." FORCE)
set(CMAKE_CXX_FLAGS_DEBUGSANADD "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address " CACHE STRING "Flags used by the compiler during Debug builds with Sanitizer for address." FORCE)
set(CMAKE_CXX_FLAGS_DEBUGSANUNDEF "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined" CACHE STRING "Flags used by the compiler during Debug builds with Sanitizer for undefineds." FORCE)
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "Flags used by the compiler during Release builds." FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "Flags used by the compiler during RelWithDebInfo builds." FORCE)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "Flags used by the compiler during Debug builds." FORCE)
set(CMAKE_C_FLAGS_DEBUGSANADD "${CMAKE_CXX_FLAGS_DEBUGSANADD}" CACHE STRING "Flags used by the compiler during Debug builds with Sanitizer for address." FORCE)
set(CMAKE_C_FLAGS_DEBUGSANUNDEF "${CMAKE_CXX_FLAGS_DEBUGSANUNDEF}" CACHE STRING "Flags used by the compiler during Debug builds with Sanitizer for undefineds." FORCE)

set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "RelWithDebInfo" "DebugSanAdd" "DebugSanUndef")
endif()
mark_as_advanced(CMAKE_CXX_FLAGS_DEBUGSANADD CMAKE_CXX_FLAGS_DEBUGSANUNDEF)
mark_as_advanced(CMAKE_C_FLAGS_DEBUGSANADD CMAKE_C_FLAGS_DEBUGSANUNDEF)
mark_as_advanced(CMAKE_EXECUTABLE_FORMAT CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT)


find_package(PkgConfig REQUIRED)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}")


########################################################################################################################
# Add a custom command that produces version.cc, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.cc
           ${CMAKE_CURRENT_BINARY_DIR}/_version.cc
    COMMAND ${CMAKE_COMMAND} -P
            ${CMAKE_CURRENT_SOURCE_DIR}/version.cmake)

########################################################################################################################
# OpenMP
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

########################################################################################################################
# Pthreads
find_package(Threads REQUIRED)

########################################################################################################################
# mpi support
option(ENABLE_MPI "Enable MPI support" ON)
if(ENABLE_MPI)
  find_package(MPI)
  if(MPI_CXX_FOUND)
    message(STATUS "MPI_CXX_COMPILER = ${MPI_CXX_COMPILER}")
    message(STATUS "MPI_CXX_INCLUDE_PATH = ${MPI_CXX_INCLUDE_PATH}")
    message(STATUS "MPI_CXX_LIBRARIES = ${MPI_CXX_LIBRARIES}")
  endif(MPI_CXX_FOUND)
endif(ENABLE_MPI)

########################################################################################################################
# floating point precision
set (
  CODE_PRECISION "DOUBLE"
  CACHE STRING "Floating point type used for internal computations and FFTs"
)
set_property (
  CACHE CODE_PRECISION
  PROPERTY STRINGS FLOAT DOUBLE LONGDOUBLE
)

########################################################################################################################
# convolver type, right now only orszag or naive
set (
  CONVOLVER_TYPE "ORSZAG"
  CACHE STRING "Convolution algorithm to be used for LPT terms (Naive=no dealiasing, Orszag=dealiased)"
)
set_property (
  CACHE CONVOLVER_TYPE
  PROPERTY STRINGS ORSZAG NAIVE
)

set (
  CONVOLVER_TYPE_PNG "ORSZAG"
  CACHE STRING "Convolution algorithm to be used for PNG terms (Naive=no dealiasing, Orszag=dealiased)"
)
set_property (
  CACHE CONVOLVER_TYPE_PNG
  PROPERTY STRINGS ORSZAG NAIVE
)

########################################################################################################################
# PLT options, right now only on/off
option(ENABLE_PLT "Enable PLT (particle linear theory) corrections" OFF)


########################################################################################################################
# FFTW
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
if(ENABLE_MPI)
  find_package(FFTW3 COMPONENTS SINGLE DOUBLE LONGDOUBLE OPENMP THREADS MPI)
else()
  find_package(FFTW3 COMPONENTS SINGLE DOUBLE LONGDOUBLE OPENMP THREADS)
endif(ENABLE_MPI)
mark_as_advanced(FFTW3_SINGLE_MPI_LIBRARY FFTW3_SINGLE_OPENMP_LIBRARY FFTW3_SINGLE_SERIAL_LIBRARY FFTW3_SINGLE_THREADS_LIBRARY)
mark_as_advanced(FFTW3_DOUBLE_MPI_LIBRARY FFTW3_DOUBLE_OPENMP_LIBRARY FFTW3_DOUBLE_SERIAL_LIBRARY FFTW3_DOUBLE_THREADS_LIBRARY)
mark_as_advanced(FFTW3_LONGDOUBLE_MPI_LIBRARY FFTW3_LONGDOUBLE_OPENMP_LIBRARY FFTW3_LONGDOUBLE_SERIAL_LIBRARY FFTW3_LONGDOUBLE_THREADS_LIBRARY)
mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_MPI_INCLUDE_DIR)
mark_as_advanced(pkgcfg_lib_PC_FFTW_fftw3)

########################################################################################################################
# GSL
find_package(GSL REQUIRED)
mark_as_advanced(pkgcfg_lib_GSL_gsl pkgcfg_lib_GSL_gslcblas pkgcfg_lib_GSL_m)

########################################################################################################################
# HDF5
find_package(HDF5 REQUIRED)
mark_as_advanced(HDF5_C_LIBRARY_dl HDF5_C_LIBRARY_hdf5 HDF5_C_LIBRARY_m HDF5_C_LIBRARY_pthread HDF5_C_LIBRARY_z HDF5_C_LIBRARY_sz)

########################################################################################################################
# PANPHASIA
option(ENABLE_PANPHASIA "Enable PANPHASIA random number generator" ON)
if(ENABLE_PANPHASIA)
enable_language(Fortran)
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132 -implicit-none")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132 -fimplicit-none")
endif()
endif(ENABLE_PANPHASIA)
########################################################################################################################
# PANPHASIA HO (High-Order, new version)
# option(ENABLE_PANPHASIA_HO "Enable PANPHASIA-HO random number generator" ON)
########################################################################################################################
# INCLUDES
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/external)

# SOURCES
# get all the *.cc files in the subfolders
file( GLOB SOURCES
  ${PROJECT_SOURCE_DIR}/src/*.cc
)
# add the auto generated version file
list (APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/version.cc")

# PLUGINS
# get all the *.cc files in the plugin subfolder
file( GLOB PLUGINS
  ${PROJECT_SOURCE_DIR}/src/plugins/*.cc
)

if(ENABLE_PANPHASIA)
if(ENABLE_MPI)
list (APPEND SOURCES 
  ${PROJECT_SOURCE_DIR}/external/panphasia/panphasia_routines.f
  ${PROJECT_SOURCE_DIR}/external/panphasia/generic_lecuyer.f90
  #${PROJECT_SOURCE_DIR}/external/panphasia_ho/main.c
  ${PROJECT_SOURCE_DIR}/external/panphasia_ho/high_order_panphasia_routines.c
  ${PROJECT_SOURCE_DIR}/external/panphasia_ho/pan_mpi_routines.c
  ${PROJECT_SOURCE_DIR}/external/panphasia_ho/uniform_rand_threefry4x64.c
)
# target_include_directories(${PRGNAME} PRIVATE ${PROJECT_SOURCE_DIR}/external/panphasia_ho)
endif()
endif()

########################################################################################################################
# DOXYGEN
# This will generate Doxygen documentation when you run make doc_doxygen in your build directory.
option(BUILD_DOCUMENTATION "Build doxygen documentation" OFF)
if(BUILD_DOCUMENTATION)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

        add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    endif(DOXYGEN_FOUND)
endif(BUILD_DOCUMENTATION)

########################################################################################################################
# project configuration header
configure_file(
  ${PROJECT_SOURCE_DIR}/include/cmake_config.hh.in
  ${PROJECT_SOURCE_DIR}/include/cmake_config.hh
)

########################################################################################################################
# executable and linking
add_executable(${PRGNAME} ${SOURCES} ${PLUGINS})
# target_setup_class(${PRGNAME})

set_target_properties(${PRGNAME} PROPERTIES CXX_STANDARD 14)

########################################################################################################################
# mpi flags and precision set up
if(MPI_CXX_FOUND)
  if(CODE_PRECISION STREQUAL "FLOAT")
    if(FFTW3_SINGLE_MPI_FOUND)
      target_link_libraries(${PRGNAME} PRIVATE FFTW3::FFTW3_SINGLE_MPI)
      target_compile_definitions(${PRGNAME} PRIVATE "USE_FFTW_MPI")
    else()
      message(SEND_ERROR "MPI enabled but FFTW3 library not found with MPI support for single precision!")
    endif()
  elseif(CODE_PRECISION STREQUAL "DOUBLE")
    if(FFTW3_DOUBLE_MPI_FOUND)
      target_link_libraries(${PRGNAME} PRIVATE FFTW3::FFTW3_DOUBLE_MPI)
      target_compile_definitions(${PRGNAME} PRIVATE "USE_FFTW_MPI")
    else()
      message(SEND_ERROR "MPI enabled but FFTW3 library not found with MPI support for double precision!")
    endif()
  elseif(CODE_PRECISION STREQUAL "LONGDOUBLE")
    if(FFTW3_LONGDOUBLE_MPI_FOUND)
      target_link_libraries(${PRGNAME} PRIVATE FFTW3::FFTW3_LONGDOUBLE_MPI)
      target_compile_definitions(${PRGNAME} PRIVATE "USE_FFTW_MPI")
    else()
      message(SEND_ERROR "MPI enabled but FFTW3 library not found with MPI support for long double precision!")
    endif()
  endif()

  target_link_libraries(${PRGNAME} PRIVATE MPI::MPI_CXX)
  target_compile_definitions(${PRGNAME} PRIVATE "USE_MPI")
endif(MPI_CXX_FOUND)
# else()
if(CODE_PRECISION STREQUAL "FLOAT")
  if(FFTW3_SINGLE_THREADS_FOUND)
    target_link_libraries(${PRGNAME} PRIVATE FFTW3::FFTW3_SINGLE_THREADS)
    target_compile_definitions(${PRGNAME} PRIVATE "USE_FFTW_THREADS")
  endif()
  target_link_libraries(${PRGNAME} PRIVATE FFTW3::FFTW3_SINGLE_SERIAL) 
elseif(CODE_PRECISION STREQUAL "DOUBLE")
  if(FFTW3_DOUBLE_THREADS_FOUND) 
    target_link_libraries(${PRGNAME} PRIVATE FFTW3::FFTW3_DOUBLE_THREADS)
    target_compile_definitions(${PRGNAME} PRIVATE "USE_FFTW_THREADS")
  endif()
  target_link_libraries(${PRGNAME} PRIVATE FFTW3::FFTW3_DOUBLE_SERIAL) 
elseif(CODE_PRECISION STREQUAL "LONGDOUBLE")
  if(FFTW3_LONGDOUBLE_THREADS_FOUND) 
    target_link_libraries(${PRGNAME} PRIVATE FFTW3::FFTW3_LONGDOUBLE_THREADS)
    target_compile_definitions(${PRGNAME} PRIVATE "USE_FFTW_THREADS")
  endif()
  target_link_libraries(${PRGNAME} PRIVATE FFTW3::FFTW3_LONGDOUBLE_SERIAL)
endif()


if(HDF5_FOUND)
  # target_link_libraries(${PRGNAME} ${HDF5_C_LIBRARY_DIRS})
  target_link_libraries(${PRGNAME} PRIVATE ${HDF5_LIBRARIES})
  target_include_directories(${PRGNAME} PRIVATE ${HDF5_INCLUDE_DIRS})
  target_compile_definitions(${PRGNAME} PRIVATE "USE_HDF5")
endif(HDF5_FOUND)

if(ENABLE_PANPHASIA)
  target_compile_definitions(${PRGNAME} PRIVATE "USE_PANPHASIA")
endif(ENABLE_PANPHASIA)

# if(ENABLE_PANPHASIA_HO)
#   target_compile_definitions(${PRGNAME} PRIVATE "USE_PANPHASIA_HO")
# endif(ENABLE_PANPHASIA_HO)

if(ENABLE_PLT)
  target_compile_definitions(${PRGNAME} PRIVATE "ENABLE_PLT")
endif(ENABLE_PLT)

target_link_libraries(${PRGNAME} PRIVATE GSL::gsl)


########################################################################################################################
# include CLASS
option(ENABLE_CLASS "Enable CLASS support (as a submodule)." ON)
  
if(ENABLE_CLASS)
  include(${CMAKE_CURRENT_SOURCE_DIR}/external/class.cmake)
  target_link_libraries(${PRGNAME} PRIVATE class::libclass_cpp)
  target_compile_definitions(${PRGNAME} PRIVATE "USE_CLASS")
endif()

########################################################################################################################
# include GenericIO (only if MPI available)
include(CMakeDependentOption)
cmake_dependent_option(ENABLE_GENERICIO "Enable GenericIO (HACC) output support" off "ENABLE_MPI;MPI_CXX_FOUND" off)

if(ENABLE_GENERICIO)
  include(${CMAKE_CURRENT_SOURCE_DIR}/external/genericio.cmake)
  target_link_libraries(${PRGNAME} PRIVATE genericio::genericio_mpi)
  target_compile_definitions(${PRGNAME} PRIVATE "ENABLE_GENERICIO")
endif()

########################################################################################################################
# CPU specifics
if(APPLE)
    # prevent linking problems on recent MacOS, enable exceptions (disabled otherwise)
    target_link_options(${PRGNAME} PRIVATE "-fexceptions" "-Wl,-no_warn_duplicate_libraries,-ld_classic" )
    set(EXE_LINKER_FLAGS "${EXE_LINKER_FLAGS} -Wl,-keep_dwarf_unwind")
    set(SHARED_LINKER_FLAGS "${SHARED_LINKER_FLAGS} -Wl,-keep_dwarf_unwind")
    set(MODULE_LINKER_FLAGS "${MODULE_LINKER_FLAGS} -Wl,-keep_dwarf_unwind")
endif()
