
set(sources 
    TPZMixedElasticityND.cpp 
    TPZCompelDiscScaled.h 
    TPZCompelDiscScaled.cpp 
    meshgen.h 
    meshgen.cpp 
    pzmixedelasmat.h 
    pzmixedelasmat.cpp
    TPZMixedElasticityCMeshCreator.cpp
    TPZMixedElasticityUtils.cpp)
file(GLOB headers *.h)


set(Mymain main.cpp )
add_executable(MixedElasticitySquare ${Mymain} ${sources})
target_link_libraries(MixedElasticitySquare NeoPZ::pz)


set(MymainOsc main_MixedElastOscillat.cpp )
add_executable(MixedElasticityOscillatory ${MymainOsc} ${sources})
target_link_libraries(MixedElasticityOscillatory NeoPZ::pz ${OpenCV_LIBS})


set(MymainG mainGirkmann.cpp ${sources} ${headers})
add_executable(MixedElasticityGirkmann ${MymainG})
target_link_libraries(MixedElasticityGirkmann NeoPZ::pz)
configure_file(Girkmann.msh Girkmann.msh COPYONLY)


set(mainMixedElasticityFive main-five.cpp ${sources} ${headers})
add_executable(MixedElasticityFive ${mainMixedElasticityFive})
target_link_libraries(MixedElasticityFive NeoPZ::pz)

set(mainMixedElasticitySeven main-seven.cpp ${sources} ${headers})
add_executable(MixedElasticitySeven ${mainMixedElasticitySeven})
target_link_libraries(MixedElasticitySeven NeoPZ::pz)

set(MHMSimple main_simple.cpp common_files.h common_files.cpp)
add_executable(MHMTest ${MHMSimple})
target_link_libraries(MHMTest NeoPZ::pz )

set(mainH1 main-h1.cpp)
add_executable(ElasticityH1 ${mainH1})
target_link_libraries(ElasticityH1 NeoPZ::pz)

set(mainYotov main-Yotov.cpp ${sources} ${headers})
add_executable(Yotov ${mainYotov})
target_link_libraries(Yotov NeoPZ::pz)

set(mainYotov7spaces main-Yotov7spaces.cpp ${sources} ${headers})
add_executable(Yotov7 ${mainYotov7spaces})
target_link_libraries(Yotov7 NeoPZ::pz)