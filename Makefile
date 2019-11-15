DIRS = ./
SRC_DIR	=	./src
TGT_DIR	=	./bin
OBJ_DIR	=	./obj
INC_DIR	=	./include

# compiler
CC_SERIAL     =  g++
CC_MPI        =  mpicxx
OMPI_CC       =  mpicxx
OMPI_CLINKER  =  mpicxx
MPI_LIB       =  ${MPI_HOME}
OPTFLAGS	  =  -Wall -O3 -g
CPPFLAGS	  =  -std=c++17 -DDEBUG 

LMP_LIB     =  ${HOME}/Softwares/lammps/src/
CINCLUDE 	=  -I${INC_DIR}
CDLINK  	=  -lm -lmpi -lpthread  -lgmp -larmadillo 

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
