LLVM_CONFIG = ${HOME}/.local/opt/llvm/bin/llvm-config

LLVM_BIN := $(shell ${LLVM_CONFIG} --bindir)

#LLVM_CXXFLAGS := $(shell ${LLVM_CONFIG} --cxxflags)
#LLVM_LDFLAGS := $(shell ${LLVM_CONFIG} --ldflags)
#LLVM_LIBS := $(shell ${LLVM_CONFIG} --libs core orcjit native)

#CXX=${LLVM_BIN}/clang++
#CXX=clang++

CXXFLAGS = -std=c++14 -g -Wno-deprecated-register
#CXXFLAGS = -std=c++14 -g -Wno-deprecated-register -fsanitize=address -fno-omit-frame-pointer
#LDFLAGS = -fsanitize=address -g
#LDFLAGS = -fsanitize=address -g

SRCS = main.cpp \
       node.cpp \
       option.cpp \
       lpp.y.cpp \
       context.cpp \
	   util.cpp \
       lpp.l.cpp

OBJS = $(SRCS:.cpp=.o)
TARGET = lcc-prep

all: ${TARGET}

lpp.y.cpp: lpp.y node.h
	bison -d -o $@ $<

lpp.y.hpp: lpp.y.cpp node.h

lpp.l.cpp: lpp.l lpp.y.hpp
	flex -o $@ $<

node.o: node.cpp node.h
main.o: main.cpp node.h option.h
option.o: option.cpp option.h
context.o: context.cpp


${TARGET}: $(OBJS) node.h
	${CXX} ${LDFLAGS}  -o $@ $(OBJS)

test: ${TARGET} sample.lpp
	./${TARGET} sample.lpp

lcc:
	python ../compiler/lcc.py -L ../../data -S output sample.lpp

lcc-run:
	python cell.py

lcc-clean:
	rm -f output/* cell.py cell.*.fa

clean:
	rm -f $(OBJS) lpp.y.cpp lpp.l.cpp lpp.y.hpp ${TARGET}
