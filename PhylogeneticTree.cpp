#include "PhylogeneticTree.h"

PhylogeneticTree::PhylogeneticTree(string filename) : _n(0)
{
    root = load_tree(filename.c_str());
    Init(root);
    N.resize(B.size());
}

PhylogeneticTree::~PhylogeneticTree()
{
    delete root;
}

void PhylogeneticTree::Init(newick_node* node)
{
    unsigned id = node->taxon;
    if (id >= B.size())
        B.resize(id + 1, false);
    B[id] = true;

    if (!node->child)
        L.push_back(node);

    for (newick_child* child = node->child; child; child = child->next)
    {
        child->node->parent = node;
        Init(child->node);
    }
    ++_n;
}
