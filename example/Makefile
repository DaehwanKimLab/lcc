CC=../src/compiler/lcc.py

LIBS=../data

LDFLAGS= -L$(LIBS)

SRCS=ecoli.lpp

TARGET=mkcell

$(TARGET).py: $(SRCS)
	$(CC) $(LDFLAGS) -o $(TARGET) $(SRCS)

clean:
	rm -f $(TARGET).py $(TARGET).*.fa *.npy

run:
	python $(TARGET).py
