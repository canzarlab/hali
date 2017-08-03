#include "PhylogeneticTree.h"
#include <fstream>

DAG::DAG(const char* f1, const char* f2, bool y)
{
    msvs C, P;
    msn M;
    ifstream ef(f1), lf(f2);
    string n1, n2, s;
    while (ef >> n1 >> n2)
    {
        if (y) ef >> s;
        P[n1].push_back(n2);
        P[n2];
        C[n2].push_back(n1);
    }

    for (auto& i : P)
        if (i.second.empty())
            root = load_dag(i.first, M, C);

    while (lf >> s >> n1)
        clade[M[n1]].push_back(s);

    Init(root);
}

/*
#include <algorithm>

for (string& j : L[i])
    if (find(L[r].begin(), L[r].end(), j) == L[r].end())
        L[r].push_back(j);
*/

newick_node* DAG::load_dag(const string& r, msn& M, msvs& C)
{
    if (newick_node* node = M[r])
        return node;

    newick_child* child = nullptr;
    newick_child** childptr = &child;
    for (const string& i : C[r])
    {
        *childptr = new newick_child(load_dag(i, M, C));
        childptr = &(*childptr)->next;
    }
    return M[r] = new newick_node("", 0, child);
}

void Graph::Init(newick_node* node)
{
    if (!node->child)
    {
        L.push_back(node);
        Leaf(node);
    }

    node->taxon = to_string(node->taxoni = _n++);
    for (newick_child* child = node->child; child; child = child->next)
    {
        newick_node* cnode = child->node;
        newick_parent** parentptr = &cnode->parent;
        if (!cnode->parent)
            Init(cnode);

        while (*parentptr) parentptr = &(*parentptr)->next;
        *parentptr = new newick_parent(node);
        Child(node, child->node);
    }
    clade[node].sort();
}

Tree::Tree(const char* f1)
{
    Init(root = load_tree(f1));
}

void Tree::Leaf(newick_node* node)
{
    clade[node].push_back(node->taxon);
}

void Tree::Child(newick_node* node, newick_node* child)
{
    list<string> &cl = clade[node], &cr = clade[child];
    cl.insert(cl.end(), cr.begin(), cr.end());
}
