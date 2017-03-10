#include "PhylogeneticTree.h"

PhylogeneticTree::PhylogeneticTree(string filename) : _n(0)
{
    root = load_tree(filename.c_str());
    init(root);
}

void PhylogeneticTree::init(newick_node* node)
{
    unsigned id = stoi(node->taxon);
    if (id >= A.size())
        A.resize(id + 1, false);
    A[id] = true;

    if (!node->child)
        L.push_back(node);
    
    for (newick_child* child = node->child; child; child = child->next)
    {
        child->node->parent = node;
        init(child->node);
    }
    ++_n;
}
