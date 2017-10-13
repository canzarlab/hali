#include "Timer.h"
#include "Greedy.h"
#include "LP.h"
#include "BnB.h"
#include <iostream>

int LP::cf;

Solver* MakeSolver(Graph& t1, Graph& t2, string d, double k, bool dag, int s)
{
    if (s == 0)
        return new Greedy(t1, t2, d, k, dag);
    else if (s == 1)
        return new LP(t1, t2, d, k, dag);
    return new BnB(t1, t2, d, k, dag);
}

DAG* MakeDAG(const char* f1, const char* f2, int s)
{
    if (s == 0)
        return (new GDAG(f1, f2))->BuildNetwork();
    return (new LDAG(f1, f2))->BuildNetwork();
}

int main(int argc, char** argv)
{
    bool dag = false;
    Graph *t1, *t2;
    if (argc != 7 && argc != 9)
    {
        cout << "tree usage: " << argv[0] << " <filename.newick> <filename.newick> <align> <0=matching 1=crossing 2=strict> <j=jaccard s=symdif> <k>" << endl;
        cout << "dag usage: " << argv[0] << " <yeastnet> <mapping> <go> <align> <0=matching 1=crossing 2=strict> <j=jaccard s=symdif> <k> <0=greedy 1=fractional 2=bnb>" << endl;
        return EXIT_FAILURE;
    }
    else if (argc == 7)
    {
        clog << "Comparing trees " << argv[1] << " " << argv[2] << endl;
        t1 = new Tree(argv[1]);
        t2 = new Tree(argv[2]);
    }
    else
    {
        int s = stoi(argv[argc - 1]);
        clog << "Comparing dags " << argv[1] << " " << argv[3] << endl;
        t1 = MakeDAG(argv[1], argv[2], s);
        t2 = MakeDAG(argv[3], nullptr, s);
        dag = true;
    }

    const char* out = argv[argc - 5];
    LP::cf = stoi(argv[argc - 4]);
    string d = argv[argc - 3];
    double k = stod(argv[argc - 2]);
    assert(LP::cf >= 0 && LP::cf <= 2);
    assert(d == "j" || d == "s");

    Solver* solver = MakeSolver(*t1, *t2, d, k, dag, stoi(argv[argc - 1]));

    Timer T;
    T.start();
    solver->Solve();
    T.stop();
    clog << "TOTAL TIME : \t\t" << T.secs() << " secs" << endl;
    solver->WriteSolution(out);
    delete t1;
    delete t2;
    delete solver;
}
