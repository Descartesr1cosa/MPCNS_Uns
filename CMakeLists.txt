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

# ---- libraries: one directory per concern, strictly layered -----------------
#   core -> complex -> geometry -> generate / partition / parallel / refine / io
add_subdirectory(src/uns)

# ---- applications ------------------------------------------------------------
add_subdirectory(Preprocess)
if(UNS_BUILD_MERCURY)
  find_package(MPI COMPONENTS CXX)
  if(MPI_CXX_FOUND)
    add_subdirectory(MercuryUns)
  else()
    message(STATUS "MPI not found: MercuryUns and uns_comm are skipped")
  endif()
endif()

if(UNS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
