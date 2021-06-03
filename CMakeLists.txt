cmake_minimum_required (VERSION 3.4)
project (PhotonEmission LANGUAGES C CXX)

set(libname "hydro_photonEmission")
set(exename "hydro_photonEmission.e")

# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message(FATAL_ERROR "Insufficient gcc version")
    endif()
endif()

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")
set (CMAKE_INSTALL_PREFIX "${CMAKE_HOME_DIRECTORY}")
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
message("prefix = " ${CMAKE_INSTALL_PREFIX})
message("rpath = " ${CMAKE_INSTALL_RPATH})

option (link_with_lib "Link executable with library" ON)
option (unittest "Build Unit tests" OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -Wall -std=c++11")

if (unittest)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDOCTEST_CONFIG_IMPLEMENT_WITH_MAIN")
else (unittest)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDOCTEST_CONFIG_DISABLE")
endif (unittest)

add_subdirectory (src)
