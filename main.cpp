#include "Timer.h"
#include "Greedy.h"
#include "LP.h"
#include "BnB.h"
#include "BnG.h"
#include "LPInt.h"
#include "LPCP.h"
#include "LPFInt.h"
#include "Similarity.h"
#include <iostream>

int Solver::cf;
bool Solver::tt;

Solver* MakeSolver(Graph& t1, Graph& t2, int argc, char** argv)
{
    int s = stoi(argv[argc - 1]);
    Solver::cf = stoi(argv[4 + (argc == 9) + 2 * (argc == 12)]);
    Solver::tt = argc == 12;
    string d = argv[5 + (argc == 9) + 2 * (argc == 12)];
    double k = stod(argv[6 + (argc == 9) + 2 * (argc == 12)]);
    double c = (s != 2) ? 0 : stod(argv[8 + 2 * (argc == 12)]);
    var_eps = (argc == 9) ? 0 : stod(argv[7 + 2 * (argc == 12)]);

    assert(LP::cf >= 0 && LP::cf <= 2);
    assert(d == "j" || d == "s");
    assert(s >= 0 && s <= 6);

    if (s == 0)
        return new Greedy(t1, t2, d, k, argc == 9);
    else if (s == 1)
        return new LP(t1, t2, d, k, argc == 9);
    else if (s == 2)
        return new BnB(t1, t2, d, k, argc == 9, c);
    else if (s == 3)
        return new LPCP(t1, t2, d, k, argc == 9);
    else if (s == 4)
        return new BnG(t1, t2, d, k, argc == 9);
    else if (s == 5)
        return new LPInt(t1, t2, d, k, argc == 9);
    return new LPFInt(t1, t2, d, k, argc == 9);
}

Graph* MakeDAG(const char* f1, const char* f2, int s)
{
    return (s == 0) ? new DAG(f1, f2) : new LDAG(f1, f2);
}

pair<Graph*, Graph*> MakeGraphs(int argc, char** argv)
{
    if (argc == 10)
        return make_pair(new Tree(argv[1]), new Tree(argv[2]));
    else if (argc == 12)
        return make_pair(new Tree(argv[1], argv[2]), new Tree(argv[3], argv[4]));

    int s = stoi(argv[argc - 1]);
    return make_pair(MakeDAG(argv[1], argv[2], s), MakeDAG(argv[3], nullptr, s));
}

int main(int argc, char** argv)
{
    if (argc < 9 || argc > 12)
    {
        cout << "tree usage: " << argv[0] << " <filename.newick> <filename.newick> <align> <0=matching 1=crossing 2=strict> <j=jaccard s=symdif> <k> <vareps> <coneps> <0=greedy 1=fractional 2=bnb 3=covering-packing 4=bng 5=integral 6=fract/int>" << endl;
        cout << "dag usage: " << argv[0] << " <yeastnet> <mapping> <go> <align> <0=matching 1=crossing 2=strict> <j=jaccard s=symdif> <k> <0=greedy 1=fractional 2=bnb 3=covering-packing 4=bng 5=integral 6=fract/int>" << endl;
        cout << "tree usage (2): " << argv[0] << " <tree> <map> <tree> <map> <align> <0=matching 1=crossing 2=strict> <j=jaccard s=symdif> <k> <vareps> <coneps> <0=greedy 1=fractional 2=bnb 3=covering-packing 4=bng 5=integral 6=fract/int>" << endl;
        return EXIT_FAILURE;
    }

    Timer T;
    T.start();
    Graph *t1, *t2;
    tie(t1, t2) = MakeGraphs(argc, argv);
    Solver* solver = MakeSolver(*t1, *t2, argc, argv);
    solver->Solve(argv[3 + (argc == 9) + 2 * (argc == 12)]);
    T.stop();
    clog << "TOTAL TIME : \t\t" << T.secs() << " secs" << endl;
    delete t1;
    delete t2;
    delete solver;
}
