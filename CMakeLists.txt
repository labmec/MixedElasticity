# @file SBFem/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.10.0)

project (MixedElasticity)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Enables NeoPZ library
find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../install_neopz/ ${CMAKE_SOURCE_DIR}/install_neopz/ "../neopz_install")


IF(APPLE)
   # Mac OS X specific code
   add_definitions(-DMACOSX )
   set(MACOSX MACOSX)
   FIND_LIBRARY(ACCELERATE_LIB Accelerate)
   link_libraries(${ACCELERATE_LIB})
#   SET(OperatingSystem "Mac OS X")

   # Making pz "seeable" in Xcode
   file(GLOB_RECURSE PZ_SRC
        "../neopz/*.h"
        "../neopz/*.cpp"
        "../NeoPzDevelop/*.h"
        "../NeoPzDevelop/*.cpp"
   )
   add_custom_target( GlobalFiles SOURCES ${PZ_SRC})
ENDIF(APPLE)

message("Processor architecture = ${CMAKE_SYSTEM_PROCESSOR}")

add_subdirectory(Mixed2D)
add_subdirectory(MHM-Elas-3D)

#Enabling nightly tests
ENABLE_TESTING()
include(CTest)
