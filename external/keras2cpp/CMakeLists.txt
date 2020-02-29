project(keras2cpp)
set(CMAKE_CXX_STANDARD 17)

set(KERAS2CPP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/layers)
set(KERAS2CPP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(LIBKERAS2CPP_SRCS
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

add_library(keras2cpp STATIC ${LIBKERAS2CPP_SRCS})