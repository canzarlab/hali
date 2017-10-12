#ifndef SOLVER_H
#define SOLVER_H

#include "Graph.h"
#include <algorithm>
#include <cmath>

class Solver
{
public:
    Solver(Graph& t1, Graph& t2, string d, double k, bool dag);

    virtual ~Solver() {}
    virtual void Solve() = 0;
    virtual void WriteSolution(string fileName) = 0;

protected:
    Graph &t1, &t2;
    string d;
    double k;
    bool dag;

    template <class F>
    void DFSLeft(newick_node* node, vb& P, F f);

private:
    template <class F>
    void DFSRight(newick_node* nodel, newick_node* noder, vb& Q, F f);
    double JaccardSim(const list<string>& L1, const list<string>& L2) const;
    double SymdifSim(const list<string>& L1, const list<string>& L2) const;
};

template <class F>
void Solver::DFSLeft(newick_node* node, vb& P, F f)
{
    P[node->taxoni] = true;
    {
        vb Q(t2.GetNumNodes());
        DFSRight(node, t2.GetRoot(), Q, f);
    }
    for (newick_child* child = node->child; child; child = child->next)
        if (!P[child->node->taxoni])
            DFSLeft(child->node, P, f);
}

template <class F>
void Solver::DFSRight(newick_node* nodel, newick_node* noder, vb& Q, F f)
{
    Q[noder->taxoni] = true;
    if (d == "j")
        f(nodel, noder, JaccardSim(t1.clade[nodel], t2.clade[noder]));
    else
        f(nodel, noder, SymdifSim(t1.clade[nodel], t2.clade[noder]));

    for (newick_child* child = noder->child; child; child = child->next)
        if (!Q[child->node->taxoni])
            DFSRight(nodel, child->node, Q, f);
}

#endif