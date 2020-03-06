include	makefiles/macos.include

MAKETARGET	=	${TGT_DIR}/kn.exe
# source cpp files
SRC_FILES   = $(wildcard ${SRC_DIR}/*.cpp)

# parallel
PARALLEL = MPI
CC = ${CC_MPI}

# Rules
# all objects depend on headers
OBJECTS 	:=	$(patsubst ${SRC_DIR}/%.cpp,${OBJ_DIR}/%.o,${SRC_FILES})

# %.o: %.cpp
$(OBJ_DIR)/%.o: ${SRC_DIR}/%.cpp
	${CC} ${CINCLUDE} ${OPTFLAGS} ${CPPFLAGS} -c -o $@ $<

${MAKETARGET}:${OBJECTS}

all: ${MAKETARGET}
	${CC}  -o  ${MAKETARGET} ${OBJECTS} ${CINCLUDE} ${CDLINK} ${OPTFLAGS} ${CPPFLAGS}

clean:
	rm -f  ${OBJ_DIR}/*.o
