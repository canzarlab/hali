#ifndef CROSSING_CONSTRAINT_H
#define CROSSING_CONSTRAINT_H

#include "Constraint.h"

class CrossingConstraint : Constraint
{
public:
    CrossingConstraint(Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp);
    
    int AddTriplets(vector<ET>& Triplets, int nr_rows);
private:
    template <class F>
    pair<newick_node*, double> GetMaxPC(newick_node* nodel, newick_node* noder, F f, bool s);
    pair<newick_node*, double> GetMaxChild(newick_node* nodel, newick_node* noder);
    pair<newick_node*, double> GetMaxParent(newick_node* nodel, newick_node* noder);

    void DFSLeft(newick_node* node, vb& C);
    void KahnLeft(newick_node* node);
    double DFSRight(newick_node* node, newick_node* nodel);
    void Reconstruct(vii& P, newick_node* nodel, newick_node* noder);
    double& GetDP(newick_node* nodel, newick_node* noder, bool s);

    vi PA;
    vvd DP;
};

template <class F>
pair<newick_node*, double> CrossingConstraint::GetMaxPC(newick_node* nodel, newick_node* noder, F f, bool s)
{
    double mx = 0;
    newick_node* mc = nullptr;
    for (newick_child* pc = f(noder); pc; pc = pc->next)
    {
        double cw = GetDP(nodel, pc->node, s);
        if (cw >= mx)
        {
            mx = cw;
            mc = pc->node;
        }
    }
    return make_pair(mc, mx);
}

#endif
