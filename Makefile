#MKLROOT = /opt/intel/mkl

CXX = g++
#CFLAGS = -g -std=c++11 -pthread -m64 # for Valgrind
CFLAGS = -O3 -std=c++11 -pthread -m64 #-DEIGEN_USE_MKL_ALL -I$(MKLROOT)/include
INCL = -I.
VPATH = src
BINARIES = hali filter bgen conflicts convert
GENO_OBJS = test_hali.o geno/augmentedLagrangian.o geno/lbfgsb.o geno/lineSearch.o
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

hali: $(GENO_OBJS) Graph.o read_csv.o Greedy.o Solver.o LP.o AntichainConstraint.o Constraint.o IndependentSetConstraint.o CrossingConstraint.o newick.o BnB.o BnG.o Similarity.o LPInt.o LPCP.o LPFInt.o  Parallel.o
	$(CXX) -o $@ $^ $(CFLAGS) #-L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core -ltbb -lstdc++ -lpthread -lm -ldl

filter: filter.o
	$(CXX) -o $@ $^ $(CFLAGS)

bgen: generator.o newick.o
	$(CXX) -o $@ $^ $(CFLAGS)

conflicts: conflicts.o newick.o Similarity.o
	$(CXX) -o $@ $^ $(CFLAGS)

convert: convert.o
	$(CXX) -o $@ $^ $(CFLAGS) -lboost_graph
