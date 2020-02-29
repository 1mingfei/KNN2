if (APPLE)
    message("Looking for mpi on macOS")
    set(FIND_MPI_PATHS /usr/local/)
elseif (UNIX)
    message("Looking for mpi on Linux")
    set(FIND_MPI_PATHS /usr/lib/x86_64-linux-gnu/openmpi/)
endif()

find_path(LIBMPI_INCLUDE_DIR
        NAME mpi.h
        PATHS ${FIND_MPI_PATHS}/include)

find_library(LIBMPI_LIBRARY
        NAMES mpi
        PATHS ${FIND_MPI_PATHS}/lib)

include_directories(${LIBMPI_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${LIBMPI_LIBRARY})