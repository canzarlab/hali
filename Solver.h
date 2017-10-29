#ifndef SOLVER_H
#define SOLVER_H

#include "Graph.h"
#include "Similarity.h"

class Solver
{
public:
    Solver(Graph& t1, Graph& t2, string d, double k, bool dag);

    virtual ~Solver() {}
    virtual void Solve(string filename) = 0;
    virtual void WriteSolution(string fileName) = 0;

    static int cf;
    static bool tt;
protected:
    Graph &t1, &t2;
    string d;
    double k;
    bool dag;

    template <class F>
    void DFSLeft(newick_node* node, vb& P, F f);
    void PrintScore(double weight);
private:
    template <class F>
    void DFSRight(newick_node* nodel, newick_node* noder, vb& Q, F f);
    int GetMax(newick_node* node, int& hmax) const;
    double GetMax(Graph& t1, Graph& t2, newick_node* root, newick_node* rnode) const;
    double CalcY(Graph& t1, Graph& t2, newick_node* root) const;
    double TumorDist(double weight) const;
    double SymdifDist(double weight) const;
    double JaccardDist(double weight) const;
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
        f(nodel, noder, JaccardSim(t1.clade[nodel], t2.clade[noder], k));
    else
        f(nodel, noder, SymdifSim(t1.clade[nodel], t2.clade[noder]));

    for (newick_child* child = noder->child; child; child = child->next)
        if (!Q[child->node->taxoni])
            DFSRight(nodel, child->node, Q, f);
}

#endif
