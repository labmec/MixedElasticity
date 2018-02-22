file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main Main.cpp ${headers})

add_executable(MixedElasticity2D ${main})

target_link_libraries(MixedElasticity2D pz)

