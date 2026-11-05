# Each library only depends on the ones listed before it.
function(uns_add_lib name)
  add_library(${name} STATIC ${ARGN})
  uns_target_defaults(${name})
endfunction()

uns_add_lib(uns_complex  complex/Incidence.cpp complex/CellComplex.cpp complex/ComplexBuilder.cpp complex/Validate.cpp)
uns_add_lib(uns_geometry geometry/Geometry.cpp)
target_link_libraries(uns_geometry PUBLIC uns_complex)
uns_add_lib(uns_generate generate/HexMesh.cpp generate/BoxMesh.cpp generate/CubedSphereMesh.cpp)
target_link_libraries(uns_generate PUBLIC uns_complex)
