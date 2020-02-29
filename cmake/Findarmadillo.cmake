set(FIND_ARMADILLO_PATHS /usr/local/)

find_path(LIBARMADILLO_INCLUDE_DIR armadillo
        PATH_SUFFIXES include
        PATHS ${FIND_ARMADILLO_PATHS})
find_library(LIBARMADILLO_LIBRARY
        NAMES armadillo
        PATH_SUFFIXES lib
        PATHS ${FIND_ARMADILLO_PATHS})