CXX = g++
CFLAGS = -O2 -std=c++11 -pthread
INCL = -I.
BINARIES = solver filter bgen conflicts
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

solver: $(GENO_OBJS) Graph.o Greedy.o Solver.o LP.o AntichainConstraint.o Constraint.o IndependentSetConstraint.o CrossingConstraint.o newick.o BnB.o Similarity.o LPInt.o
	$(CXX) -o $@ $^ $(CFLAGS)

filter: filter.o
	$(CXX) -o $@ $^ $(CFLAGS)

bgen: generator.o newick.o
	$(CXX) -o $@ $^ $(CFLAGS)

conflicts: conflicts.o newick.o Similarity.o
	$(CXX) -o $@ $^ $(CFLAGS)
