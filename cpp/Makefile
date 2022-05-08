SHELL='bash'
#
# Bridges - PSC
#
# Intel Compilers are loaded by default
# You will need to specifically switch to GNU Modules
# With with `modules.sh`
#

CC = g++
OPENMP = -fopenmp
CFLAGS = -O3
LIBS =


TARGETS = serial

all:	$(TARGETS)

serial: main.o Mdfun.o
	$(CC) -o $@ $(LIBS) main.o Mdfun.o

main.o: main.cpp Mdfun.h
	$(CC) -c $(CFLAGS) main.cpp
Mdfun.o: Mdfun.cpp Mdfun.h
	$(CC) -c $(CFLAGS) Mdfun.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
