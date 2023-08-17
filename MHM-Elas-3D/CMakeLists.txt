set(MESHESDIR ${CMAKE_SOURCE_DIR}/MHM-Elas-3D/meshes)
set (CONFIG_H "${CMAKE_SOURCE_DIR}/MHM-Elas-3D/meshpath_config.h")
configure_file ("${CMAKE_SOURCE_DIR}/MHM-Elas-3D/meshpath_config.h.in" ${CONFIG_H})


add_executable(mhm-elas-3d voronoi_mixed_elas.cpp)
target_link_libraries(mhm-elas-3d NeoPZ::pz)
