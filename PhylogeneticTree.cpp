#include "PhylogeneticTree.h"
#include <fstream>

DAG::DAG(const char* f1, const char* f2, bool y)
{
    msn M;
    root = load_dag(f1, y, M);
    ifstream lf(f2);
    string s, n1;
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
