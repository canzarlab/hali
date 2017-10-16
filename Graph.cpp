#include "Graph.h"
#include <fstream>
#include <limits>

DAG::DAG(const char* f1, const char* f2)
{
    msn M;
    Init(root = load_dag(f1, f2, clade, M));
}

LDAG::LDAG(const char* f1, const char* f2) : DAG(f1, f2)
{
    size_t SZ = _n * 2 + 2;
    G.resize(SZ);
    for (int j = 0; j < NR_THREADS; ++j)
        R[j].resize(SZ, vd(SZ));

    int S = SZ - 2, T = SZ - 1;
    for (int i = 0; i < _n; ++i)
    {
        G[S].push_back(i);
        G[i].push_back(S);
        G[T].push_back(i + _n);
        G[i + _n].push_back(T);
    }
}

DAG* DAG::BuildNetwork()
{
    vvb C(_n, vb(_n));
    for (newick_node* leaf : L)
        BuildNetwork(leaf, leaf, C);
    return this;
}

void DAG::BuildNetwork(newick_node* node, newick_node* rnode, vvb& C)
{
    int l = rnode->taxoni;
    int i = node->taxoni;
    if (l != i)
        Relation(l, i);

    C[i][l] = true;
    for (newick_parent* parent = node->parent; parent; parent = parent->next)
    {
        newick_node* pn = parent->node;
        int pnt = pn->taxoni;
        if (!C[pnt][l])
            BuildNetwork(pn, rnode, C);
        if (!C[pnt][pnt])
            BuildNetwork(pn, pn, C);
    }
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
    ls &cl = clade[node], &cr = clade[child];
    cl.insert(cl.end(), cr.begin(), cr.end());
}

GDAG::GDAG(const char* f1, const char* f2) : DAG(f1, f2), D(_n, vb(_n))
{
}

void LDAG::Relation(int l, int i)
{
    for (int j = 0; j < NR_THREADS; ++j)
        R[j][l][i + _n] = numeric_limits<double>::infinity();
    G[l].push_back(i + _n);
    G[i + _n].push_back(l);
}

void GDAG::Relation(int l, int i)
{
    D[l][i] = true;
}
