DIRS = ./
SRC_DIR	=	./src
TGT_DIR	=	./bin
OBJ_DIR	=	./obj
INC_DIR	=	./include
K2C_INC_DIR	=	./external/keras2cpp
LIB_DIR	=	./lib
LMP_LIB     =  ${HOME}/Softwares/lammps/src/

# compiler
CC_SERIAL     =  g++
CC_MPI        =  mpicxx
OMPI_CC       =  mpicxx
OMPI_CLINKER  =  mpicxx
MPI_LIB       =  ${MPI_HOME}
OPTFLAGS	  =  -Wall -O3 -g
CPPFLAGS	  =  -std=c++17 #-DDEBUG_SELECT_TRAP

#include and lib while compile
CINCLUDE 	=  -I${INC_DIR} -I${K2C_INC_DIR}
CDLINK  	=  -L${LIB_DIR} -lm -lmpi -lpthread  -lgmp -larmadillo -lkeras2cpp

MAKETARGET	=	${TGT_DIR}/kn.exe

SRC_FILES   = $(wildcard ${SRC_DIR}/*.cpp)

# parallel
PARALLEL = MPI
CC = ${CC_MPI}

# Rules
# all objects depend on headers
# OBJECTS := $(subst .cpp,.o,${SRC_FILES})
OBJECTS 	:=	$(patsubst ${SRC_DIR}/%.cpp,${OBJ_DIR}/%.o,${SRC_FILES})

# %.o: %.cpp
$(OBJ_DIR)/%.o: ${SRC_DIR}/%.cpp
	${CC} ${CINCLUDE} ${OPTFLAGS} ${CPPFLAGS} -c -o $@ $<

${MAKETARGET}:${OBJECTS}

all: ${MAKETARGET}
	${CC}  -o  ${MAKETARGET}  ${OBJECTS}  ${CDLINK}  ${OPTFLAGS}  ${CPPFLAGS}

clean:
	rm -f  ${OBJ_DIR}/*.o
