if (APPLE)
    message("Working on macOS")
    set(FIND_ARMADILLO_PATHS /usr/local/)
elseif (UNIX)
    message("Working on Linux")
    set(FIND_ARMADILLO_PATHS /usr/)
endif()

find_path(LIBARMADILLO_INCLUDE_DIR armadillo
        PATH_SUFFIXES include
        PATHS ${FIND_ARMADILLO_PATHS})
find_library(LIBARMADILLO_LIBRARY
        NAMES armadillo
        PATH_SUFFIXES lib
        PATHS ${FIND_ARMADILLO_PATHS})