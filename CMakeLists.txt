cmake_minimum_required(VERSION 3.16)
project(MPCNS_Uns LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(UNS_BUILD_TESTS     "Build unit/regression tests"            ON)
option(UNS_BUILD_MERCURY   "Build the MercuryUns solver (needs MPI)" ON)
option(UNS_ENABLE_METIS    "Use METIS for graph partitioning"        OFF)
option(UNS_ENABLE_PETSC    "Use PETSc (implicit solves)"             OFF)

include(cmake/UnsWarnings.cmake)

# Windows (MSVC or MinGW, with MS-MPI): only the MPI runtime (msmpi.dll, from the MS-MPI redistributable) stays external
if(WIN32)
  add_compile_definitions(NOMINMAX _USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS)
  if(MSVC)
    add_compile_options(/utf-8 /EHsc)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>") # static CRT: no vcruntime DLLs to ship
  else()
    add_link_options(-static-libgcc -static-libstdc++)
  endif()
endif()

# RPATH: the binaries (build tree and installed copies) carry the directory of the MPI library, so they run without LD_LIBRARY_PATH.
# Set explicitly (instead of letting CMake rewrite the build RPATH at install time, which is fragile with MPI's own -Wl,-rpath flags).
find_package(MPI COMPONENTS CXX)
if(MPI_CXX_FOUND AND NOT WIN32)
  set(_mpi_dirs "")
  foreach(_l IN LISTS MPI_CXX_LIBRARIES)
    if(IS_ABSOLUTE "${_l}")
      get_filename_component(_d "${_l}" DIRECTORY)
      list(APPEND _mpi_dirs "${_d}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES _mpi_dirs)
  set(CMAKE_INSTALL_RPATH "${_mpi_dirs}")
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()

# ---- libraries: one directory per concern, strictly layered -----------------
#   core -> complex -> geometry -> generate / partition / parallel / refine / io
add_subdirectory(src/uns)

if(UNS_BUILD_TESTS)
  enable_testing()
endif()

# ---- applications (Common = shared MHD blocks, MercuryUns = Mercury, MarsUns = Mars) ------------------------------------------------------------
add_subdirectory(Preprocess)
if(UNS_BUILD_MERCURY)
  find_package(MPI COMPONENTS CXX)
  if(MPI_CXX_FOUND)
    add_subdirectory(Common)
    add_subdirectory(MercuryUns)
    add_subdirectory(MarsUns)
  else()
    message(STATUS "MPI not found: MercuryUns, MarsUns and uns_comm are skipped")
  endif()
endif()

if(UNS_BUILD_TESTS)
  add_subdirectory(tests)
endif()

# `make install` / `cmake --install build --prefix <work-dir>`: ready-to-run CASE folders (see cases/README.md)
include(cmake/InstallCase.cmake)
