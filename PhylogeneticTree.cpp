#include "PhylogeneticTree.h"

PhylogeneticTree::PhylogeneticTree(string filename) : _n(0)
{
    Init(root = load_tree(filename.c_str()));
}

PhylogeneticTree::~PhylogeneticTree()
{
    delete root;
}

void PhylogeneticTree::Init(newick_node* node)
{
    if (!node->child)
    {
        L.push_back(node);
        clade[node].push_back(node->taxon);
    }

    node->taxon = to_string(node->taxoni = _n++);
    for (newick_child* child = node->child; child; child = child->next)
    {
        child->node->parent = node;
        Init(child->node);
        list<string> &cl = clade[node], &cr = clade[child->node];
        cl.insert(cl.end(), cr.begin(), cr.end());
    }
    clade[node].sort();
}
