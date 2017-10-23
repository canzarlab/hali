#include "Timer.h"
#include "Greedy.h"
#include "LP.h"
#include "BnB.h"
#include "BnG.h"
#include "LPInt.h"
#include <iostream>

int LP::cf;

Solver* MakeSolver(Graph& t1, Graph& t2, int argc, char** argv)
{
    bool dag = (argc == 9);
    LP::cf = stoi(argv[argc - 4]);
    string d = argv[argc - 3];
    double k = stod(argv[argc - 2]);
    int s = stoi(argv[argc - 1]);
    assert(LP::cf >= 0 && LP::cf <= 2);
    assert(d == "j" || d == "s");
    assert(s >= 0 && s <= 4);

    if (s == 0)
        return new Greedy(t1, t2, d, k, dag);
    else if (s == 1)
        return new LP(t1, t2, d, k, dag);
    else if (s == 2)
        return new BnB(t1, t2, d, k, dag);
	else if (s == 3)
    	return new LPInt(t1, t2, d, k, dag);
	return new BnG(t1, t2, d, k, dag);
}

Graph* MakeDAG(const char* f1, const char* f2, int s)
{
    return (s == 0) ? new DAG(f1, f2) : new LDAG(f1, f2);
}

pair<Graph*, Graph*> MakeGraphs(int argc, char** argv)
{
    if (argc == 8)
        return make_pair(new Tree(argv[1]), new Tree(argv[2]));
    else if (argc == 10)
        return make_pair(new Tree(argv[1], argv[2]), new Tree(argv[3], argv[4]));

    int s = stoi(argv[argc - 1]);
    return make_pair(MakeDAG(argv[1], argv[2], s), MakeDAG(argv[3], nullptr, s));
}

int main(int argc, char** argv)
{
    if (argc < 8 || argc > 10)
    {
        cout << "tree usage: " << argv[0] << " <filename.newick> <filename.newick> <align> <0=matching 1=crossing 2=strict> <j=jaccard s=symdif> <k> <0=greedy 1=fractional 2=bnb 3=integral 4=bng>" << endl;
        cout << "dag usage: " << argv[0] << " <yeastnet> <mapping> <go> <align> <0=matching 1=crossing 2=strict> <j=jaccard s=symdif> <k> <0=greedy 1=fractional 2=bnb 3=integral>" << endl;
        cout << "tree usage (2): " << argv[0] << " <tree> <map> <tree> <map> <align> <0=matching 1=crossing 2=strict> <j=jaccard s=symdif> <k> <0=greedy 1=fractional 2=bnb 3=integral>" << endl;
        return EXIT_FAILURE;
    }

    Timer T;
    T.start();
    Graph *t1, *t2;
    tie(t1, t2) = MakeGraphs(argc, argv);
    Solver* solver = MakeSolver(*t1, *t2, argc, argv);
    solver->Solve(argv[argc - 5]);
    T.stop();
    clog << "TOTAL TIME : \t\t" << T.secs() << " secs" << endl;
    delete t1;
    delete t2;
    delete solver;
}
