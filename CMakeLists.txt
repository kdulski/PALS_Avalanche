cmake_minimum_required(VERSION 2.6)
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#include(JPetFramework-Utilities)

set(projectName PAv.exe)
set(projectBinary ${projectName}.x)

project(${projectName} CXX) # using only C++

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

add_definitions(-std=c++11 -Wall -Wunused-parameter) # using C++11 standard

################################################################################
# Handling of ROOT 
################################################################################

## Find CERN ROOT
# first try to find ROOT 6 or ROOT 5 compiled with CMake
# if present, such ROOT should be possible to find without explicit setting of any paths
foreach(mode QUIET REQUIRED)
  find_package(ROOT 5 ${mode} COMPONENTS
    Hist
    Physics
    RIO
    Thread
    Tree
    )
  if(ROOT_USE_FILE)
    include(${ROOT_USE_FILE})
  endif()
  if(ROOT_FOUND)
    break()
  endif()
  # if we failed with CMake based find_package, use fallback
  # root-config version from cmake/fallback/
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/fallback)
endforeach()
include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_DEFINITIONS})
message(STATUS "ROOT version: ${ROOT_VERSION}")
if(ROOT_FOUND AND NOT ROOT_VERSION VERSION_LESS "6.0")
  message(FATAL_ERROR "ROOT 6.0 is not compatible")
endif()

add_definitions(${ROOT_DEFINITIONS})

find_package(Boost 1.50.0 REQUIRED COMPONENTS system filesystem)

include_directories(${Boost_INCLUDE_DIRS})
add_definitions(${Boost_DEFINITIONS})

file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_executable(${projectName} ${SOURCES} ${HEADERS} )
target_link_libraries(${projectName} ${ROOT_LIBRARIES} ${Boost_LIBRARIES})