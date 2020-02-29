if (APPLE)
    message("Working on macOS")
    set(FIND_ARMADILLO_PATHS /usr/local/)
elseif (UNIX)
    message("Working on Linux")
    set(FIND_ARMADILLO_PATHS /usr/)
endif()

find_path(LIBMPI_INCLUDE_DIR mpi.h
        PATH_SUFFIXES include
        PATHS ${FIND_MPI_PATHS})

find_library(LIBMPI_LIBRARY
        NAMES mpi
        PATH_SUFFIXES lib
        PATHS ${FIND_MPI_PATHS})