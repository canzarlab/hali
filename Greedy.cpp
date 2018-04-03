#include "Greedy.h"
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
using namespace std::placeholders;

Greedy::Greedy(Graph& t1, Graph& t2, string d, double k, bool dag) : Solver(t1, t2, d, k, dag), A(t1.GetNumNodes(), vd(t2.GetNumNodes()))
{
    vb P(t1.GetNumNodes());
    DFSLeft(t1.GetRoot(), P, [&](newick_node* nodel, newick_node* noder, double w)
    {
        if (w != 0.0 && (dag || nodel->parent) && (dag || nodel->child) && (dag || noder->parent) && (dag || noder->child))
            E.emplace_back(nodel->taxoni, noder->taxoni, w);
    });
    sort(E.begin(), E.end(), [](const iid& a, const iid& b){return get<2>(a) > get<2>(b);});
}

void Greedy::Solve(string filename)
{
    for (const iid& e : E)
        if (all_of(M.begin(), M.end(), bind(&Greedy::CC, this, _1, cref(e))))
            M.push_back(e), A[get<0>(e)][get<1>(e)] = get<2>(e);
    WriteSolution(filename);
}

void Greedy::WriteSolution(string fileName)
{
    ofstream sol_file(fileName);
    double weight = 0.0;
    for (const iid& e : M)
        sol_file << get<0>(e) << " " << get<1>(e) << " 1\n", weight += get<2>(e);
    PrintScore(weight);
}

float Greedy::GetSolution()
{
	double weight = 0.0;
    for (const iid& e : M)
        weight += get<2>(e);
	return weight;
}

bool Greedy::CC(const iid& a, const iid& b) const
{
    int i = get<0>(a), j = get<0>(b);
    int x = get<1>(a), y = get<1>(b);
    return IsNotInConflict(i, j, x, y);
}
