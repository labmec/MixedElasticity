
add_executable(Rigid-Body-Rotation main-rigid-body-rotation.cpp)
target_link_libraries(Rigid-Body-Rotation PRIVATE MixedElastoDynamics_lib)

add_executable(MHM-Rigid-Body-Rotation main-MHM-rigid-body-rotation.cpp)
target_link_libraries(MHM-Rigid-Body-Rotation PRIVATE MixedElastoDynamics_lib)

