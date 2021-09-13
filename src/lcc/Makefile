LLVM_CONFIG = ${HOME}/.local/opt/llvm/bin/llvm-config

LLVM_BIN := $(shell ${LLVM_CONFIG} --bindir)

#LLVM_CXXFLAGS := $(shell ${LLVM_CONFIG} --cxxflags)
#LLVM_LDFLAGS := $(shell ${LLVM_CONFIG} --ldflags)
#LLVM_LIBS := $(shell ${LLVM_CONFIG} --libs core orcjit native)

#CXX=${LLVM_BIN}/clang++

CXXFLAGS = -std=c++11


SRCS = main.cpp \
	   lpp.y.cpp \
	   lpp.l.cpp

OBJS = $(SRCS:.cpp=.o)
TARGET = parser

all: ${TARGET}

lpp.y.cpp: lpp.y
	bison -d -o $@ $^

lpp.y.hpp: lpp.y.cpp

lpp.l.cpp: lpp.l lpp.y.hpp
	flex -o $@ $<

${TARGET}: $(OBJS)
	${CXX} -o $@ $(OBJS)

test: ${TARGET} sample.lpp
	./${TARGET} < sample.lpp

clean:
	rm -f $(OBJS) lpp.y.cpp lpp.l.cpp lpp.y.hpp ${TARGET}
