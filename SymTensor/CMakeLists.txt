
set(mainSym main-sym.cpp TPZMixedSymElasticityND.h TPZMixedSymElasticityND.cpp TPZInterfaceSymTensor.h TPZInterfaceSymTensor.cpp)



add_executable(SymTensor ${mainSym})
target_link_libraries(SymTensor NeoPZ::pz)
configure_file(Quad.msh Quad.msh COPYONLY)
configure_file(Triangle.msh Triangle.msh COPYONLY)
configure_file(TriangleBC.msh TriangleBC.msh COPYONLY)
