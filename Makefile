CXX = g++
CFLAGS = -Wall -O2 -std=c++11
INCL = -I Eigen
BINARIES = genoSolver youngSolver
GENO_OBJS = genoSolver.o geno/augmentedLagrangian.o geno/lbfgsb.o geno/lineSearch.o
YOUNG_OBJS = main.o young/matrix.o young/solver.o young/get_time.o
TEST_BINARIES = 
HEADERS = 

build: $(BINARIES) 


clean:
	rm -f $(BINARIES) $(TEST_BINARIES)
	rm -f *.o
	rm -f geno/*.o
	rm -f young/*.o
	rm -f *~
	rm -f core


young/%.o: young/%.cpp 
	$(CXX) -c $< -o $@ $(CFLAGS) $(INCL)

geno/%.o: geno/%.cpp 
	$(CXX) -c $< -o $@ $(CFLAGS) $(INCL)

%.o: %.cpp $(HEADERS)
	$(CXX) -c $< $(CFLAGS) $(INCL)
	
genoSolver: $(GENO_OBJS) PhylogeneticTree.o newick.o 
	$(CXX) -o $@ $^ $(CFLAGS)

youngSolver: main.o $(YOUNG_OBJS) PhylogeneticTree.o newick.o
	$(CXX) -o $@ $^ $(CFLAGS)