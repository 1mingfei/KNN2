set(FIND_KERAS2CPP_PATHS ${PROJECT_SOURCE_DIR})


find_path(LIBKERAS2CPP_INCLUDE_DIR model.h tensor.h
        PATH_SUFFIXES keras2cpp
        PATHS "${FIND_KERAS2CPP_PATHS}/external")

find_library(LIBKERAS2CPP_LIBRARY
        NAMES keras2cpp
        PATH_SUFFIXES lib
        PATHS "${FIND_KERAS2CPP_PATHS}")