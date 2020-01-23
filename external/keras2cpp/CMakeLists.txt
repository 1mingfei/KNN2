cmake_minimum_required(VERSION 3.9)

project(keras2cpp)

set(CMAKE_CXX_STANDARD 17)
#
#if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
#      -W -Wall -Wextra -Wpedantic \
#      -Waggressive-loop-optimizations \
#      -Wcast-align -Wcast-qual \
#      -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond \
#      -Wfloat-equal -Wformat=2 -Wformat-signedness -Wframe-larger-than=32768 \
#      -Wlogical-op \
#      -Wnull-dereference \
#      -Wodr -Wold-style-cast \
#      -Wshadow=local -Wshift-overflow=2 -Wstrict-aliasing=2 -Wsuggest-final-methods -Wsuggest-final-types -Wsync-nand \
#      -Wtrampolines \
#      -Wuseless-cast -Wno-unused-but-set-parameter \
#      -Wwrite-strings"
#            CACHE INTERNAL ""
#            )
#endif()
#
#if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
#      -Weverything \
#      -Wno-class-varargs -Wno-padded \
#      -Wc++17-compat \
#      -Wno-c++98-compat -Wno-c++98-compat-pedantic"
#            # -Wno-c++14-extensions
#            # -Wno-c++17-extensions
#            # -Wno-switch-enum -Wno-unused-macros
#            CACHE INTERNAL ""
#            )
#endif()
#
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g"
#        CACHE INTERNAL ""
#        )
#
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} \
#    -g -pg -no-pie \
#    -O3 -ffast-math -fno-rtti"
#        CACHE INTERNAL ""
#        )
#set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO
#        "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} -pg"
#        CACHE INTERNAL ""
#        )
#
## rolled loops works faster here
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
#    -O3 -ffast-math -fno-rtti"
#        CACHE INTERNAL ""
#        )
#
#if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#    set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -Os"
#            CACHE INTERNAL "")
#endif()
#
#if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#    set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -Oz"
#            CACHE INTERNAL "")
#endif()
#
#set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} \
#    -O3 -ffast-math -fno-math-errno -fno-unroll-loops \
#    -fno-rtti -fno-stack-protector -fno-ident \
#    -fomit-frame-pointer -ffunction-sections -fdata-sections \
#    -fmerge-all-constants"
#        CACHE INTERNAL ""
#        )

# project configuration

set (KERAS2CPP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/layers)
set (KERAS2CPP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(KERAS2CPP_SRCS
        ${KERAS2CPP_SOURCE_DIR}/model.cc
        ${KERAS2CPP_SOURCE_DIR}/tensor.cc
        ${KERAS2CPP_SOURCE_DIR}/utils.cc
        ${KERAS2CPP_SOURCE_DIR}/baseLayer.cc
        ${KERAS2CPP_INCLUDE_DIR}/activation.cc
        ${KERAS2CPP_INCLUDE_DIR}/conv1d.cc
        ${KERAS2CPP_INCLUDE_DIR}/conv2d.cc
        ${KERAS2CPP_INCLUDE_DIR}/dense.cc
        ${KERAS2CPP_INCLUDE_DIR}/elu.cc
        ${KERAS2CPP_INCLUDE_DIR}/embedding.cc
        ${KERAS2CPP_INCLUDE_DIR}/flatten.cc
        ${KERAS2CPP_INCLUDE_DIR}/lstm.cc
        ${KERAS2CPP_INCLUDE_DIR}/locally1d.cc
        ${KERAS2CPP_INCLUDE_DIR}/locally2d.cc
        ${KERAS2CPP_INCLUDE_DIR}/maxPooling2d.cc
        ${KERAS2CPP_INCLUDE_DIR}/batchNormalization.cc
        )
include_directories("${KERAS2CPP_INCLUDE_DIR}")
include_directories("${KERAS2CPP_SOURCE_DIR}")

add_library(keras2cpp ${KERAS2CPP_SRCS})