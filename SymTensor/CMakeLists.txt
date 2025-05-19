
set(mainSym main-sym.cpp TPZMixedSymElasticityND.h TPZMixedSymElasticityND.cpp TPZInterfaceSymTensor.h TPZInterfaceSymTensor.cpp MeshConditioning.h MeshConditioning.cpp)

set(mainHexagon main-hexagon.cpp TPZMixedSymElasticityND.h TPZMixedSymElasticityND.cpp TPZInterfaceSymTensor.h TPZInterfaceSymTensor.cpp MeshConditioning.h MeshConditioning.cpp)


add_executable(SymTensor ${mainSym})
target_link_libraries(SymTensor NeoPZ::pz)
add_executable(SymHexagon ${mainHexagon})
target_link_libraries(SymHexagon NeoPZ::pz)
configure_file(Quad.msh Quad.msh COPYONLY)
configure_file(Triangle.msh Triangle.msh COPYONLY)
configure_file(TriangleBC.msh TriangleBC.msh COPYONLY)
configure_file(honeycomb.1.msh honeycomb.1.msh COPYONLY)
configure_file(honeycomb.2.msh honeycomb.2.msh COPYONLY)
configure_file(honeycomb.4.msh honeycomb.4.msh COPYONLY)
configure_file(honeycomb.8.msh honeycomb.8.msh COPYONLY)

