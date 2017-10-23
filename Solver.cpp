#include "Solver.h"
#include <iostream>

Solver::Solver(Graph& t1, Graph& t2, string d, double k, bool dag) : t1(t1), t2(t2), d(d), k(k), dag(dag)
{
}

void Solver::PrintScore(double weight)
{
    if (dag)
        cout << weight << " ";
    else
        cout << ((d == "j") ? JaccardDist(weight) : SymdifDist(weight)) << " ";
}

int Solver::GetMax(newick_node* node, int& hmax) const
{
    int sum = 0;
    for (newick_child* child = node->child; child; child = child->next)
        sum += GetMax(child->node, hmax);
    hmax += sum;
    return node->child ? sum : 1;
}

double Solver::SymdifDist(double weight) const
{
    int max1 = 0, max2 = 0;
    int r1 = GetMax(t1.GetRoot(), max1);
    int r2 = GetMax(t2.GetRoot(), max2);
    return max1 + max2 - r1 - r2 - weight;
}

double Solver::JaccardDist(double weight) const
{
    return t1.GetNumNodes() - t1.L.size() - 1 + t2.GetNumNodes() - 1 - t2.L.size() - weight;
}
