CC = g++
LD = g++

PYVERSION = 2.7
PYPATH = /System/Library/Frameworks/Python.framework/Versions/$(PYVERSION)

CFLAGS = -Wall -ansi -pedantic -O3 -pthread -std=c++17 `root-config --cflags` -fPIC
LIBS = -stdlib=libc++ `root-config --libs` -L$(PYPATH)/lib/python$(PYVERSION) -lpython
INCLUDE = -Iinclude -I$(PYPATH)/include/python$(PYVERSION)

OBJS  = $(patsubst src/%.cc,lib/%.o,$(wildcard src/*.cc))
EXECS = $(patsubst exe/%.cc,bin/%,$(wildcard exe/*.cc))
EXEOBJS  = $(patsubst exe/%.cc,lib/%.o,$(wildcard exe/*.cc))

SOLIB =  lib/SRSP.so

all: $(OBJS) $(EXEOBJS) $(EXECS) $(SOLIB)


lib/SRSP.so : lib/SRSP.o lib/SRS.o
	$(LD) -Wall -shared $(LIBS) $(OBJS) -o $@

lib/%.o : src/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : exe/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@


bin/% : $(OBJS) lib/%.o
	$(LD) $(LIBS) $(OBJS) lib/$*.o -o bin/$*





clean:
	rm -f $(EXECS) lib/*.o

