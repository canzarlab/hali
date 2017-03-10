#ifndef PHYLO_TREE_H
#define PHYLO_TREE_H

#include <vector>
#include "newick.h"
using namespace std;

class PhylogeneticTree
{
public:
    PhylogeneticTree(string filename);

    long getNumNodes() const { return _n; }
    bool nodeExists(unsigned n) const { return n < A.size() ? A[n] : false; }   
    
private:
    void init(newick_node* node);
    
    vector<bool> A;
    vector<newick_node*> L;
    newick_node* root;
    long _n;
};

#endif
