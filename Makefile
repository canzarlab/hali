CXX = g++
CFLAGS = -Wall -O2 -std=c++11
INCL = -I Eigen
BINARIES = genoSolver geno2Solver
COMMON_OBJS = geno/augmentedLagrangian.o geno/lbfgsb.o geno/lineSearch.o
GENO_OBJS = genoSolver.o
GENO2_OBJS = main.o
TEST_BINARIES = 
HEADERS = 

build: $(BINARIES) 


clean:
	rm -f $(BINARIES) $(TEST_BINARIES)
	rm -f *.o
	rm -f geno/*.o
	rm -f geno2/*.o
	rm -f *~
	rm -f core


geno2/%.o: geno/%.cpp
	$(CXX) -c $< -o $@ $(CFLAGS) $(INCL)

geno/%.o: geno/%.cpp
	$(CXX) -c $< -o $@ $(CFLAGS) $(INCL)

%.o: %.cpp $(HEADERS)
	$(CXX) -c $< $(CFLAGS) $(INCL)
	
genoSolver: $(GENO_OBJS) $(COMMON_OBJS) PhylogeneticTree.o newick.o 
	$(CXX) -o $@ $^ $(CFLAGS)

geno2Solver: main.o $(GENO2_OBJS) $(COMMON_OBJS) PhylogeneticTree.o newick.o
	$(CXX) -o $@ $^ $(CFLAGS)
