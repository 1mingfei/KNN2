DIRS = ./
ODIR = ./

# compiler
CC_SERIAL     =  g++
CC_MPI        =  mpicxx
OMPI_CC       =  mpicxx
OMPI_CLINKER  =  mpicxx
MPI_LIB       =  ${MPI_HOME}
OPTFLAGS	  =  -Wall -O3 -g
CPPFLAGS	  =  -std=c++11  -DDEBUG

LMP_LIB     =  ${HOME}/Softwares/lammps/src/
CINCLUDE 	=  -I../include  
CDLINK  	=  -lm -lmpi -lpthread  -lgmp -larmadillo 

MAKETARGET 	=  kn.exe

# The source files
KNSRC 	=  ./KNHome.cpp   \
			   ./KNMain.cpp\
                           ./KNSolidSol.cpp\
                           ./KNvasp.cpp\
                           ./KNParam.cpp\
                           ./gbInCnf.cpp\
                           ./gbOutCnf.cpp\
                           ./gbBox.cpp\
                           ./FCCConfig.cpp\
                           ./neighbor.cpp

# parallel
PARALLEL = MPI
CC = ${CC_MPI}

# Rules
# all objects depend on headers
OBJECTS := $(subst .cpp,.o,${KNSRC})

%.o: %.cpp
	${CC} -c  $<  ${CINCLUDE}  ${OPTFLAGS}  ${CPPFLAGS} 

${MAKETARGET}:$(OBJECTS)

all: ${MAKETARGET}
	${CC}  -o  ${MAKETARGET}  ${OBJECTS}  ${CDLINK}  ${OPTFLAGS}  ${CPPFLAGS}

clean:
	rm -f  *.o  *.u  *~ # \#* *.V *.T *.O *.il
