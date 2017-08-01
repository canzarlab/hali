CXX = g++
CFLAGS = -Wall -O2 -std=c++11
INCL = -I Eigen
BINARIES = solver phylo2tc filter bgen
GENO_OBJS = main.o geno/augmentedLagrangian.o geno/lbfgsb.o geno/lineSearch.o
TEST_BINARIES = 
HEADERS = 

build: $(BINARIES) 


clean:
	rm -f $(BINARIES) $(TEST_BINARIES)
	rm -f *.o
	rm -f geno/*.o
	rm -f *~
	rm -f core

geno/%.o: geno/%.cpp
	$(CXX) -c $< -o $@ $(CFLAGS) $(INCL)

%.o: %.cpp $(HEADERS)
	$(CXX) -c $< $(CFLAGS) $(INCL)

solver: $(GENO_OBJS) PhylogeneticTree.o newick.o
	$(CXX) -o $@ $^ $(CFLAGS)

phylo2tc: newick.o phylo2treecomp.o
	$(CXX) -o $@ $^ $(CFLAGS)

filter: filter.o
	$(CXX) -o $@ $^ $(CFLAGS)

bgen: generator.o newick.o
	$(CXX) -o $@ $^ $(CFLAGS)
