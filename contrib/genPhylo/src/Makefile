include ../MakeVars

OBJS = Newickform.o PhylogeneticTree.o phylo2treecomp.o seqUtil.o Globals.o

phylo2tc: $(OBJS)
	$(LD) $(OBJS) -o $(BINDIR)/phylo2tc $(LDFLAGS)

%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $< -o $@ 

clean:
	rm *.o *~
