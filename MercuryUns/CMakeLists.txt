# MercuryUns : two-ion (H+, Na+) Hall-MHD of Mercury on the unstructured signed cell complex.  Mercury-specific code only;
# everything model independent lives in Common/.
add_library(mercury_core STATIC
    driver/Solver.cpp
    multigrid/Multigrid.cpp
    initial/Initial.cpp
    fluid/MercurySources.cpp
    induction/Resistive.cpp
    physics/Exosphere.cpp
    physics/Gyration.cpp
    physics/Ambipolar.cpp)
uns_target_defaults(mercury_core)
target_include_directories(mercury_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mercury_core PUBLIC mhd_common)

add_executable(MercuryUns main.cpp)
uns_target_defaults(MercuryUns)
target_link_libraries(MercuryUns PRIVATE mercury_core)

add_executable(mg_b_proto bproto.cpp)
uns_target_defaults(mg_b_proto)
target_link_libraries(mg_b_proto PRIVATE mercury_core uns_refine uns_partition uns_parallel)
