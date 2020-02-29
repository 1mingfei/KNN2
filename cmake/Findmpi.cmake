if (APPLE)
    message("Looking for mpi on macOS")
    set(FIND_ARMADILLO_PATHS /usr/local/)
elseif (UNIX)
    message("Looking for mpi on Linux")
    set(FIND_ARMADILLO_PATHS /usr/)
endif()

find_path(LIBMPI_INCLUDE_DIR mpi.h
        PATH_SUFFIXES include
        PATHS ${FIND_MPI_PATHS})

find_library(LIBMPI_LIBRARY
        NAMES mpi
        PATH_SUFFIXES lib
        PATHS ${FIND_MPI_PATHS})

include_directories(${LIBMPI_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${LIBMPI_LIBRARY})