# let us create a config file with the directory containing the input files
set(INPUTSDIR ${CMAKE_CURRENT_SOURCE_DIR}/input)

set(CONFIG_H "${CMAKE_CURRENT_SOURCE_DIR}/dirs_config.h")
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/dirs_config.h.in" ${CONFIG_H})

target_include_directories(MixedElastoDynamics_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB ELASTODYNAMICS_SRC CONFIGURE_DEPENDS "*.cpp" "*.c" "*.hpp")
file(GLOB ELASTODYNAMICS_HEADERS CONFIGURE_DEPENDS "*.h")

target_sources(MixedElastoDynamics_lib PRIVATE ${ELASTODYNAMICS_SRC} ${ELASTODYNAMICS_HEADERS})

add_subdirectory(targets)