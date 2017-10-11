#include "CrossingConstraint.h"
#include <queue>

CrossingConstraint::CrossingConstraint(Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp) : Constraint(t1, t2, K, x, swp)
{
    DP.resize(t1.GetNumNodes());
    PA.resize(t1.GetNumNodes());
    for (auto& v : DP)
        v.resize(t2.GetNumNodes());
}

int CrossingConstraint::AddTriplets(vector<ET>& Triplets, int nr_rows)
{
    int ncr = 0;
    KahnLeft(t1.GetRoot());
    for (auto node : t1.L)
    {
        vii P;
        Reconstruct(P, node, t2.GetRoot());

        double sum = 0;
        for (auto k : P)
            sum += GetWeight(k.first, k.second);

        if (sum - EPS > 1)
            AddConstraint(Triplets, nr_rows + ncr, P), ncr++;
    }
    return ncr;
}

double& CrossingConstraint::GetDP(newick_node* nodel, newick_node* noder, bool s = false)
{
    return s ? DP[noder->taxoni][nodel->taxoni] : DP[nodel->taxoni][noder->taxoni];
}

pair<newick_node*, double> CrossingConstraint::GetMaxChild(newick_node* nodel, newick_node* noder)
{
    return GetMaxPC(nodel, noder, [](newick_node* n){return n->child;}, false);
}

pair<newick_node*, double> CrossingConstraint::GetMaxParent(newick_node* nodel, newick_node* noder)
{
    return GetMaxPC(nodel, noder, [](newick_node* n){return n->parent;}, true);
}

void CrossingConstraint::DFSLeft(newick_node* node)
{
    for (newick_child* child = node->child; child; child = child->next)
    {
        PA[child->node->taxoni]++;
        DFSLeft(child->node);
    }
}

void CrossingConstraint::KahnLeft(newick_node* node)
{
    queue<newick_node*> Q;
    Q.push(node);
    DFSLeft(node);
    while (!Q.empty())
    {
        node = Q.front(); Q.pop();
        DFSRight(t2.GetRoot(), node);
        for (newick_child* child = node->child; child; child = child->next)
            if (--PA[child->node->taxoni] == 0)
                Q.push(child->node);
    }
}

double CrossingConstraint::DFSRight(newick_node* node, newick_node* nodel)
{
    double mx = 0;
    for (newick_child* child = node->child; child; child = child->next)
        mx = max(mx, DFSRight(child->node, nodel));
    mx = max(mx, GetMaxParent(node, nodel).second);
    return GetDP(nodel, node) = mx + GetWeight(nodel, node);
}

void CrossingConstraint::Reconstruct(vii& P, newick_node* nodel, newick_node* noder)
{
    double pw, cw;
    newick_node *child, *parent;
    tie(child, cw) = GetMaxChild(nodel, noder);
    tie(parent, pw) = GetMaxParent(noder, nodel);
    P.emplace_back(nodel->taxoni, noder->taxoni);
    if (nodel->parent && (!child || pw > cw))
        Reconstruct(P, parent, noder);
    else if (child && (!parent || cw >= pw))
        Reconstruct(P, nodel, child);
    else
        assert(!parent && !child);
}
