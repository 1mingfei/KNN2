cmake_minimum_required(VERSION 3.9)

project(KNN2)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
      -Wall"
            CACHE INTERNAL ""
            )
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
      -Weverything \
      -Wno-class-varargs -Wno-padded \
      -Wc++17-compat \
      -Wno-c++98-compat -Wno-c++98-compat-pedantic"
            # -Wno-c++14-extensions
            # -Wno-c++17-extensions
            # -Wno-switch-enum -Wno-unused-macros
            CACHE INTERNAL ""
            )
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb"
        CACHE INTERNAL ""
        )

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
    -O3 -ffast-math -fno-rtti"
        CACHE INTERNAL ""
        )

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB KNN2_SRC
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/include/*.h)

add_executable(${PROJECT_NAME} ${KNN2_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "kn.exe")



#include(FindPkgConfig)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

find_package(ARMADILLO REQUIRED)
include_directories(${LIBARMADILLO_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${LIBARMADILLO_LIBRARY})

# Lib keras2cpp
add_subdirectory(external/keras2cpp)
include_directories(external/keras2cpp)
target_link_libraries(${PROJECT_NAME} Keras2cpp)
