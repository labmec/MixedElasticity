
set(sources TPZMixedElasticityND.cpp TPZCompelDiscScaled.h TPZCompelDiscScaled.cpp meshgen.h meshgen.cpp pzmixedelasmat.h pzmixedelasmat.cpp)
file(GLOB headers *.h)


set(Mymain main.cpp )
add_executable(MixedElasticitySquare ${Mymain} ${sources})
target_link_libraries(MixedElasticitySquare pz)


set(MymainOsc main_MixedElastOscillat.cpp )
add_executable(MixedElasticityOscillatory ${MymainOsc} ${sources})
target_link_libraries(MixedElasticityOscillatory pz ${OpenCV_LIBS})


set(MymainG mainGirkmann.cpp ${sources} ${headers})
add_executable(MixedElasticityGirkmann ${MymainG})
target_link_libraries(MixedElasticityGirkmann pz)
configure_file(Girkmann.msh Girkmann.msh COPYONLY)


set(mainMixedElasticityFive main-five.cpp ${sources} ${headers})
add_executable(MixedElasticityFive ${mainMixedElasticityFive})
target_link_libraries(MixedElasticityFive pz)

set(mainMixedElasticitySeven main-seven.cpp ${sources} ${headers})
add_executable(MixedElasticitySeven ${mainMixedElasticitySeven})
target_link_libraries(MixedElasticitySeven pz)