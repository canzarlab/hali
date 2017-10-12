#include "Timer.h"
#include "Greedy.h"
#include "LP.h"
#include <iostream>

int LP::cf;

Solver* MakeSolver(Graph& t1, Graph& t2, string d, double k, bool dag, bool greedy)
{
    if (greedy)
        return new Greedy(t1, t2, d, k, dag);
    return new LP(t1, t2, d, k, dag);
}

DAG* MakeDAG(const char* f1, const char* f2, bool y, bool greedy)
{
    if (greedy)
        return (new GDAG(f1, f2, y))->BuildNetwork();
    return (new LDAG(f1, f2, y))->BuildNetwork();
}

int main(int argc, char** argv)
{
    bool dag = false, greedy = false;
    Graph *t1, *t2;
    const char* out;
    if (argc == 7)
    {
        clog << "Comparing trees " << argv[1] << " " << argv[2] << endl;
        t1 = new Tree(argv[1]);
        t2 = new Tree(argv[2]);
        out = argv[3];
    }
    else if (argc == 10)
    {
        clog << "Comparing dags " << argv[1] << " " << argv[3] << endl;
        greedy = stoi(argv[argc - 1]);
        t1 = MakeDAG(argv[1], argv[2], true, greedy);
        t2 = MakeDAG(argv[3], argv[4], false, greedy);
        out = argv[5];
        dag = true;
    }
    else
    {
        cout << "tree usage: " << argv[0] << " <filename.newick> <filename.newick> <align> <c> <d> <k>" << endl;
        cout << "dag usage: " << argv[0] << " <yeastnet> <mapping> <precollapse> <mapping> <align> <c> <d> <k> <g>" << endl;
        return EXIT_FAILURE;
    }
    LP::cf = stoi(argv[argc - 4]);
    string d = argv[argc - 3];
    double k = stod(argv[argc - 2]);
    assert(LP::cf >= 0 && LP::cf <= 2);
    assert(d == "j" || d == "s");

    Solver* solver = MakeSolver(*t1, *t2, d, k, dag, greedy);

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
