# mpcns_pre : generic cell-complex generation / partitioning tool
add_executable(mpcns_pre main.cpp Pipeline.cpp)
uns_target_defaults(mpcns_pre)
target_link_libraries(mpcns_pre PRIVATE uns_io uns_partition uns_parallel uns_generate uns_geometry uns_complex uns_core)
