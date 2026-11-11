# Common : model-independent MHD building blocks shared by MercuryUns and MarsUns (namespace mhd).
#   config/       run configuration (Config) read from setup files
#   physics/      Euler helpers (primitive/conservative, Rusanov), nondimensional scales
#   state/        FieldSet (named entity-numbered arrays)
#   boundary/     face / edge classification (wall, inflow, outflow)
#   fluid/        limited-gradient reconstruction, finite-volume flux right-hand side
#   induction/    constrained transport (Faraday), edge assembly, ideal UCT EMF, cell current, semi-implicit Hall/resistive solve, Hall hook
#   diagnostics/  div B
#   io/           gather of global fields (restart / coarse solution)
#   app/          common main(): setup parsing, mesh read, selftest, task dispatch
add_library(mhd_common STATIC
    config/Config.cpp
    physics/Scales.cpp
    state/FieldSet.cpp
    boundary/BoundaryInfo.cpp
    fluid/Reconstruction.cpp
    fluid/FluxRHS.cpp
    induction/Faraday.cpp
    induction/EdgeEMF.cpp
    induction/IdealEMF.cpp
    induction/Current.cpp
    induction/ImplicitInduction.cpp
    induction/HallEMF.cpp
    diagnostics/DivB.cpp
    io/GlobalFields.cpp
    app/App.cpp)
uns_target_defaults(mhd_common)
target_include_directories(mhd_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mhd_common PUBLIC uns_comm uns_io uns_geometry uns_complex uns_core)

# the MPI launcher this binary was built with: fallback of the self-launching ./MercuryUns / ./MarsUns when mpirun is not in PATH
if(MPIEXEC_EXECUTABLE)
  target_compile_definitions(mhd_common PRIVATE UNS_MPIEXEC="${MPIEXEC_EXECUTABLE}")
endif()
