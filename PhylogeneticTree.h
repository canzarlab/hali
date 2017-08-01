#ifndef PHYLO_TREE_H
#define PHYLO_TREE_H

#include <vector>
#include <map>
#include <list>
#include "newick.h"
using namespace std;

class PhylogeneticTree
{
public:
    PhylogeneticTree(string filename);
    ~PhylogeneticTree();

    int GetNumNodes() const { return _n; }
    newick_node* GetRoot() const { return root; }

    vector<newick_node*> L;
    map<newick_node*, list<string> > clade;
private:
    void Init(newick_node* node);
    newick_node* root;
    long _n;
};

#endif
