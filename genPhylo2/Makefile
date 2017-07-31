BINDIR = .
CC = g++
LD = g++
CPPFLAGS = -std=c++0x -m64 -fPIC -g -fexceptions  -DIL_STD -I..
LDFLAGS = $(LIBDIRS) -lm -lpthread -lemon
OBJS = ../newick.o phylo2treecomp.o

phylo2tc: $(OBJS)
	$(LD) $(OBJS) -o $(BINDIR)/phylo2tc $(LDFLAGS)

%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $< -o $@ 

clean:
	rm *.o phylo2tc
