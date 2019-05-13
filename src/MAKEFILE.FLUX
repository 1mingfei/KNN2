DIRS = ./
ODIR = ./

# compiler
CC_SERIAL     =  g++
CC_MPI        =  mpicxx
OMPI_CC       =  mpicxx
OMPI_CLINKER  =  mpicxx
MPI_LIB       =  ${MPI_HOME}
OPTFLAGS	  =  -Wall -O3 -g
CPPFLAGS	  =  -std=c++11 

ARMA_HOME  = ${HOME}/software/armadillo-code-8.300.x

CINCLUDE 	=  -I../include -I${ARMA_HOME}/include
CDLINK  	=  -lm -lmpi -lpthread  -lgmp -DARMA_DONT_USE_WRAPPER -lblas -llapack -L${ARMA_HOME}/lib64

MAKETARGET 	=  em.exe

# The source files
KNSRC 	=  ./KNHome.cpp   \
			   ./KNMain.cpp\
                           ./gbInCnf.cpp\
                           ./gbOutCnf.cpp\
                           ./gbBox.cpp\
                           ./KNParam.cpp\
                           ./KNAlign.cpp\
                           ./KNStat.cpp

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
