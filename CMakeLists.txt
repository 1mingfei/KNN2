cmake_minimum_required(VERSION 3.15)

# set the project name
project(KNN2)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

set(KNN2_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(KNN2_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

include_directories("${KNN2_INCLUDE_DIR}")
include_directories("${KNN2_SOURCE_DIR}")

add_subdirectory(external/keras2cpp)

file(GLOB KNN2_SRC
        ${KNN2_SOURCE_DIR}/*.cpp
        ${KNN2_INCLUDE_DIR}/*)
add_executable(KNN2 ${KNN2_SRC})

include(FindPkgConfig)

find_package(Armadillo REQUIRED)
find_package(keras2cpp REQUIRED)
include_directories(${LIBARMADILLO_INCLUDE_DIR})
include_directories(${LIBKERAS2CPP_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${LIBARMADILLO_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${LIBKERAS2CPP_LIBRARY})


