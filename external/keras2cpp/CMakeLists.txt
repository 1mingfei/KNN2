set(LIBKERAS2CPP_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/model.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/tensor.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/utils.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/baseLayer.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/activation.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/conv1d.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/conv2d.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/dense.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/elu.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/embedding.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/flatten.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/lstm.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/locally1d.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/locally2d.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/maxPooling2d.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/layers/batchNormalization.cc
        )
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
add_library(Keras2cpp STATIC ${LIBKERAS2CPP_SRCS})
