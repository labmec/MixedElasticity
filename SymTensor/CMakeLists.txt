
set(mainSym main-sym.cpp)
add_executable(SymTensor ${mainSym})
target_link_libraries(SymTensor NeoPZ::pz)

