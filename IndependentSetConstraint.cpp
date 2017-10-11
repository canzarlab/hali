#include "IndependentSetConstraint.h"

IndependentSetConstraint::IndependentSetConstraint(Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp) : Constraint(t1, t2, K, x, swp)
{
}

int IndependentSetConstraint::AddTriplets(vector<ET>& Triplets, int nr_rows)
{
    int ncr = 0;
    for (newick_node* node : t1.L)
    {
        dLN L = DFSRight(node, t2.GetRoot());
        if (L.first - EPS <= 1)
            continue;

        vii P;
        for (newick_node* noder : L.second)
            for (newick_node* nodel = node; nodel; nodel = nodel->parent ? nodel->parent->node : nullptr)
                P.emplace_back(nodel->taxoni, noder->taxoni);

        AddConstraint(Triplets, nr_rows + ncr, P);
        ++ncr;
    }
    return ncr;
}

double IndependentSetConstraint::PathSum(newick_node* nodel, newick_node* noder) const
{
    return nodel ? GetWeight(nodel, noder) + PathSum(nodel->parent ? nodel->parent->node : nullptr, noder) : 0;
}

dLN IndependentSetConstraint::DFSRight(newick_node* nodel, newick_node* noder)
{
    double w = PathSum(nodel, noder), sum = 0;
    LN V;

    for (newick_child* child = noder->child; child; child = child->next)
    {
        double ww;
        LN T;
        tie(ww, T) = DFSRight(nodel, child->node);
        sum += ww;
        V.splice(V.begin(), T);
    }

    if (sum > w)
        return make_pair(sum, V);
    return make_pair(w, LN(1, noder));
}
