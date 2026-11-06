# MarsUns : ten-ion Hall-MHD of Mars on the unstructured signed cell complex.  Mars-specific code only (see docs/Mars.md);
# the model-independent parts (flux, CT, UCT EMF, Hall implicit solve, boundaries, IO) come from Common/.
add_library(mars_core STATIC
    physics/MarsChemistry.cpp
    physics/MarsAtmosphereModel.cpp
    physics/MarsDipoles.cpp
    driver/MarsSolver.cpp)
uns_target_defaults(mars_core)
target_include_directories(mars_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mars_core PUBLIC mhd_common)

add_executable(MarsUns main.cpp)
uns_target_defaults(MarsUns)
target_link_libraries(MarsUns PRIVATE mars_core)

if(UNS_BUILD_TESTS)
  add_executable(test_mars tests/test_mars.cpp)
  uns_target_defaults(test_mars)
  target_include_directories(test_mars PRIVATE ${PROJECT_SOURCE_DIR}/tests)
  target_link_libraries(test_mars PRIVATE mars_core uns_refine uns_partition uns_generate)
  foreach(n 1 4)
    add_test(NAME test_mars_np${n}
             COMMAND ${MPIEXEC_EXECUTABLE} --allow-run-as-root --oversubscribe ${MPIEXEC_NUMPROC_FLAG} ${n} $<TARGET_FILE:test_mars>)
  endforeach()
endif()
