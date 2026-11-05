# MercuryUns : two-ion Hall-MHD on the unstructured signed cell complex.
# One sub-directory per concern; each only sees the ones listed before it.
add_library(mercury_core STATIC
    state/FieldSet.cpp
    induction/Faraday.cpp
    induction/EdgeEMF.cpp
    induction/IdealEMF.cpp
    induction/Current.cpp
    induction/Resistive.cpp
    induction/ImplicitInduction.cpp
    physics/Exosphere.cpp
    physics/Gyration.cpp
    physics/Ambipolar.cpp
    physics/Scales.cpp
    Config.cpp
    boundary/BoundaryInfo.cpp
    fluid/Reconstruction.cpp
    fluid/FluidRHS.cpp
    initial/Initial.cpp
    io/GlobalFields.cpp
    driver/Solver.cpp
    diagnostics/DivB.cpp
    hall/HallEMF.cpp)
uns_target_defaults(mercury_core)
target_include_directories(mercury_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mercury_core PUBLIC uns_comm uns_io uns_geometry uns_complex uns_core)

add_executable(MercuryUns main.cpp)
uns_target_defaults(MercuryUns)
target_link_libraries(MercuryUns PRIVATE mercury_core)
