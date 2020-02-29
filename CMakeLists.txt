cmake_minimum_required(VERSION 3.9)

# set the project name
project(KNN2)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB KNN2_SRC
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/include/*.h)
add_executable(${PROJECT_NAME} ${KNN2_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "kn.exe")




#include(FindPkgConfig)
find_package(armadillo REQUIRED)
find_package(mpi REQUIRED)
find_package(keras2cpp REQUIRED)
