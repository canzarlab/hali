#include "Greedy.h"
#include <iostream>
#include <fstream>
#include <functional>
using namespace std::placeholders;

Greedy::Greedy(Graph& t1, Graph& t2, string d, double k, bool dag) : Solver(t1, t2, d, k, dag), A(t1.GetNumNodes(), vd(t2.GetNumNodes()))
{
    vb P(t1.GetNumNodes());
    DFSLeft(t1.GetRoot(), P, [&](newick_node* nodel, newick_node* noder, double w)
    {
        if (w != 0 && (dag || nodel->parent) && nodel->child && (dag || noder->parent) && noder->child)
            E.emplace_back(nodel->taxoni, noder->taxoni, w);
    });
    sort(E.begin(), E.end(), [](const iii& a, const iii& b){return get<2>(a) > get<2>(b);});
}

void Greedy::Solve()
{
    for (iii& e : E)
        if (all_of(M.begin(), M.end(), bind(&Greedy::CC, this, _1, cref(e))))
            M.push_back(e), A[get<0>(e)][get<1>(e)] = get<2>(e);
}

void Greedy::WriteSolution(string fileName)
{
    ofstream sol_file(fileName);
    double weight = 0;
    for (vd& v : A)
    {
        for (double d : v)
            sol_file << (d != 0.) << "\t", weight += d;
        sol_file << endl;
    }
    cout << weight << " ";
}


bool Greedy::CC(const iii& a, const iii& b)
{
    GDAG &t1 = static_cast<GDAG&>(this->t1);
    GDAG &t2 = static_cast<GDAG&>(this->t2);
    int i = get<0>(a), j = get<0>(b);
    int x = get<1>(a), y = get<1>(b);
    if (i == j || x == y) return false;
    bool c2 = (t1.D[j][i] || t1.D[i][j]) == (t2.D[x][y] || t2.D[y][x]);
    bool c1 = (t1.D[i][j] == t2.D[x][y]); // assuming c2 is satisfied
    return c2 && c1;
}
