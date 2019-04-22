
set(sources TPZCompelDiscScaled.h TPZCompelDiscScaled.cpp meshgen.h meshgen.cpp pzmixedelasmat.h pzmixedelasmat.cpp)

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
