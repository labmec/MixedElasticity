file(GLOB headers *.h)
file(GLOB sources meshgen.cpp pzmixedelasmat.cpp)

set(Mymain main.cpp ${sources} ${headers})

add_executable(MixedElasticitySquare ${Mymain})

target_link_libraries(MixedElasticitySquare pz)

set(MymainG mainGirkmann.cpp ${sources} ${headers})

add_executable(MixedElasticityGirkmann ${MymainG})

target_link_libraries(MixedElasticityGirkmann pz)

configure_file(Girkmann.msh Girkmann.msh COPYONLY)
