file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main Main.cpp ${headers})

add_executable(Mixed2D ${main})

target_link_libraries(Mixed2D pz)

