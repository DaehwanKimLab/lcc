LLVMBASE = ${HOME}/.local/opt/llvm
LLVM_CONFIG = ${LLVMBASE}/bin/llvm-config

#LLVM_CPPFLAGS := $(shell ${LLVM_CONFIG} --cxxflags)
#LLVM_LDFLAGS := $(shell ${LLVM_LDFLAGS} --ldflags)
#LLVM_LIBS := $(shell ${LLVM_CONFIG} --libs core orcjit native)

SRCS = main.cpp \
	   lpp.y.cpp \
	   lpp.l.cpp

OBJS = $(SRCS:.cpp=.o)

all: parser

lpp.y.cpp: lpp.y
	bison -d -o $@ $^

lpp.y.hpp: lpp.y.cpp

lpp.l.cpp: lpp.l lpp.y.hpp
	flex -o $@ $<

%.o: %.cpp
	g++ -std=c++11 -c -o $@ $<

parser: $(OBJS)
	g++ -o $@ $(OBJS)

test: parser sample.lpp
	./parser < sample.lpp

clean:
	rm -f $(OBJS) lpp.y.cpp lpp.l.cpp lpp.y.hpp parser
