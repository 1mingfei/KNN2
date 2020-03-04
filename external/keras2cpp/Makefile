
CC_SERIAL = icpc
CPPFLAGS = --std=c++17 -I. -I./layers -Wall -Werror -MMD -O3 -mtune=core2

KRSRC =	./baseLayer.cc \
	./model.cc \
	./tensor.cc \
	./utils.cc \
	./layers/activation.cc \
	./layers/batchNormalization.cc \
	./layers/conv1d.cc \
	./layers/conv2d.cc \
	./layers/dense.cc \
	./layers/elu.cc \
	./layers/embedding.cc \
	./layers/flatten.cc \
	./layers/locally1d.cc \
	./layers/locally2d.cc \
	./layers/lstm.cc \
	./layers/maxPooling2d.cc 


CC = ${CC_SERIAL}
OBJECTS := $(subst .cc,.o,${KRSRC})

%.o: %.cc ./layers/%.cc
	${CC} -c  $<  ${CINCLUDE}  ${OPTFLAGS}  ${CPPFLAGS} 

lib:	$(OBJECTS)
	mkdir -p ../../lib
	ar rcs ../../lib/libkeras2cpp.a model.o tensor.o ${OBJECTS}                                

clean:
	rm -f  *.o *.d ./layers/*.o ./layers/*.d 


-include $(OBJ_FILES:.o=.d)
