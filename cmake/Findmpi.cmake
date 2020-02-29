set(FIND_MPI_PATHS /usr/local)

find_path(LIBMPI_INCLUDE_DIR mpi.h
        PATH_SUFFIXES include
        PATHS ${FIND_MPI_PATHS})

find_library(LIBMPI_LIBRARY
        NAMES mpi
        PATH_SUFFIXES lib
        PATHS ${FIND_MPI_PATHS})