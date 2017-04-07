#ifndef PHYLO_TREE_H
#define PHYLO_TREE_H

#include <vector>
#include "newick.h"
using namespace std;

class PhylogeneticTree
{
public:
    PhylogeneticTree(string filename);
    ~PhylogeneticTree();

    int GetNumNodes() const { return _n; }
    bool NodeExists(unsigned n) const { return n < B.size() ? B[n] : false; }
    int GetIndex(newick_node* node) const { return N[node->taxon]; }
    newick_node* GetRoot() const { return root; }

    // TODO: objediniti N i B
    vector<int> N;
    vector<newick_node*> L;
private:
    void Init(newick_node* node);
    
    vector<bool> B;
    newick_node* root;
    long _n;
};

#endif
