all: CSim


PYTHON_CONFIG=python3-config --embed

CFLAGS=$(shell $(PYTHON_CONFIG) --cflags)
LDFLAGS=$(shell $(PYTHON_CONFIG) --ldflags)

CFLAGS += -g

CSim: SimCWrap.cpp
	$(CXX) $(CFLAGS) $< $(LDFLAGS) -o $@


clean:
	rm -f CSim
