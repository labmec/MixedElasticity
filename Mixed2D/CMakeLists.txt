
add_library(mixedElast SHARED "")
target_link_libraries(mixedElast PUBLIC NeoPZ::pz)
target_include_directories(mixedElast PUBLIC ${PZ_INCLUDE_DIRS})

target_include_directories(mixedElast PUBLIC
${CMAKE_CURRENT_SOURCE_DIR})

set(sources 
    TPZMixedElasticityND.cpp 
    TPZCompelDiscScaled.h 
    TPZCompelDiscScaled.cpp 
    meshgen.h 
    meshgen.cpp 
    pzmixedelasmat.h 
    pzmixedelasmat.cpp 
    TPZMixedElasticityNDAiry.cpp 
    AiryFunctionHoledPlate.cpp
    AiryFunctionLoadedBeam.cpp
    TPZCompElDiscStress.cpp
    TPZCompElDiscStressBound.cpp
    TPZMixedElasticityCMeshCreator.cpp
    TPZMixedElasticityUtils.cpp
    )
file(GLOB headers *.h)

target_sources(mixedElast PRIVATE
  ${sources})

set(Mymain main.cpp )
add_executable(MixedElasticitySquare ${Mymain})
target_link_libraries(MixedElasticitySquare mixedElast)

set(Main1Fracture main_1fracture.cpp )
add_executable(1Fracture ${Main1Fracture})
target_link_libraries(1Fracture mixedElast)

set(MainTest main_test.cpp )
add_executable(Test ${MainTest})
target_link_libraries(Test mixedElast) 

set(MainBeam main_beam.cpp )
add_executable(Beam ${MainBeam})
target_link_libraries(Beam mixedElast) 

set(MymainOsc main_MixedElastOscillat.cpp )
add_executable(MixedElasticityOscillatory ${MymainOsc})
target_link_libraries(MixedElasticityOscillatory  ${OpenCV_LIBS} mixedElast)


set(MymainG mainGirkmann.cpp ${headers})
add_executable(MixedElasticityGirkmann ${MymainG})
target_link_libraries(MixedElasticityGirkmann mixedElast)
configure_file(Girkmann.msh Girkmann.msh COPYONLY)


set(mainMixedElasticityFive main-five.cpp ${sources} ${headers})
add_executable(MixedElasticityFive ${mainMixedElasticityFive})
target_link_libraries(MixedElasticityFive NeoPZ::pz)

set(mainMixedElasticitySeven main-seven.cpp ${sources} ${headers})
add_executable(MixedElasticitySeven ${mainMixedElasticitySeven})
target_link_libraries(MixedElasticitySeven NeoPZ::pz)

set(MainSimple main_simple.cpp )
add_executable(LinkTest ${MainSimple})
target_link_libraries(LinkTest NeoPZ::pz )

