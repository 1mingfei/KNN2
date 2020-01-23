# - Config file for the Armadillo package
# It defines the following variables
#  ARMADILLO_INCLUDE_DIRS - include directories for Armadillo
#  ARMADILLO_LIBRARY_DIRS - library directories for Armadillo (normally not used!)
#  ARMADILLO_LIBRARIES    - libraries to link against

# Tell the user project where to find our headers and libraries
set(FIND_ARMADILLO_PATHS
        /usr/local/Cellar/armadillo/9.800.3_1)

#set(ARMADILLO_INCLUDE_DIRS "/usr/local/Cellar/armadillo/9.800.3_1/include")
#set(ARMADILLO_LIBRARY_DIRS "/usr/local/Cellar/armadillo/9.800.3_1/lib")

find_path(LIBARMADILLO_INCLUDE_DIR armadillo
        PATH_SUFFIXES include
        PATHS ${FIND_ARMADILLO_PATHS})
find_library(LIBARMADILLO_LIBRARY
        NAMES armadillo
        PATH_SUFFIXES lib
        PATHS ${FIND_ARMADILLO_PATHS})
