CC = g++
LD = g++
NVCC = nvcc

PYVERSION = 2.7
PYPATH = /System/Library/Frameworks/Python.framework/Versions/$(PYVERSION)

WSTPDIR = /Applications/Mathematica.app/Contents/SystemFiles/Links/WSTP/DeveloperKit/MacOSX-x86-64/CompilerAdditions/

CFLAGS = -DCUDA -Wall -ansi -pedantic -O3 -pthread -std=c++11 -fPIC
CUDACFLAGS = -DCUDA -std=c++11
LIBS = -Llib -L$(PYPATH)/lib/python$(PYVERSION) -lpython -L$(WSTPDIR) -lWSTPi4 -stdlib=libstdc++ -framework Foundation -lc++ -L/usr/local/cuda/lib -lcuda -lcudart_static
INCLUDE = -Iinclude -I$(PYPATH)/include/python$(PYVERSION) -I$(WSTPDIR)

OBJS  = $(patsubst src/%.cc,lib/%.o,$(wildcard src/*.cc))
CUDAOBJS  = $(patsubst src/%.cu,lib/%.o,$(wildcard src/*.cu))
EXECS = $(patsubst exe/%.cc,bin/%,$(wildcard exe/*.cc))
EXEOBJS  = $(patsubst exe/%.cc,lib/%.o,$(wildcard exe/*.cc))

WSTPOBJS  = $(patsubst wstp/%.tm,lib/%_tm.o,$(wildcard wstp/*.tm))
WSTPCCS  = $(patsubst wstp/%.tm,wstp/%_tm.cc,$(wildcard wstp/*.tm))

SOLIB =  lib/SRS.so


WSPREP = $(WSTPDIR)/wsprep


all: $(WSTPCCS) $(WSTPOBJS) $(OBJS) $(CUDAOBJS) $(EXEOBJS) $(EXECS) $(SOLIB)


mma: bin/SRS_MMA


wstp/%_tm.cc : wstp/%.tm
	$(WSPREP) $< -o $@

lib/SRS.so : $(OBJS) $(CUDAOBJS) $(WSTPCCS) $(WSTPOBJS)
	$(LD) -Wall -shared $(LIBS) $(WSTPOBJS) $(OBJS) $(CUDAOBJS) -o $@

lib/%.o : src/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : src/%.cu
	$(NVCC) $(CUDACFLAGS) $(INCLUDE) -c $< -o $@


lib/%.o : exe/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : $(WSTPCCS)
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@


bin/% : $(WSTPOBJS) $(OBJS) $(CUDAOBJS) lib/%.o
	$(LD) $(LIBS) $(OBJS) $(CUDAOBJS) $(WSTPOBJS) lib/$*.o -o bin/$*





clean:
	rm -f $(EXECS) lib/*.o $(SOLIB) $(WSTPOBJS) $(WSTPCCS)

