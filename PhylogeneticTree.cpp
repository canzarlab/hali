#include "PhylogeneticTree.h"

PhylogeneticTree::PhylogeneticTree(string filename) : _n(0)
{
    root = load_tree(filename.c_str());
    Init(root);
    N.resize(B.size());
}

void PhylogeneticTree::Init(newick_node* node)
{
    unsigned id = stoi(node->taxon);
    if (id >= B.size())
        B.resize(id + 1, false);
    B[id] = true;
    
    for (newick_child* child = node->child; child; child = child->next)
    {
        child->node->parent = node;
        Init(child->node);
    }
    ++_n;
}
