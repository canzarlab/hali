#ifndef PHYLO_TREE_H
#define PHYLO_TREE_H

#include <vector>
#include <map>
#include <list>
#include "newick.h"
using namespace std;

typedef map<string, newick_node*> msn;
typedef map<string, vector<string> > msvs;

class Graph
{
public:
    Graph() : _n(0) { }
    virtual ~Graph() { delete root; };

    int GetNumNodes() const { return _n; }
    newick_node* GetRoot() const { return root; }
    newick_node* GetNode(int i) const { return N[i]; }

    vector<newick_node*> L, N;
    map<newick_node*, list<string> > clade;

protected:
    void Init(newick_node* node);
    virtual void Leaf(newick_node* node) { }
    virtual void Child(newick_node* node, newick_node* child) { }
    newick_node* root;
    long _n;

};

class DAG : public Graph
{
public:
    DAG(const char* f1, const char* f2, bool y);
    ~DAG() { }
};

class Tree : public Graph
{
public:
    Tree(const char* f1);
    ~Tree() { }

protected:
    void Leaf(newick_node* node);
    void Child(newick_node* node, newick_node* child);
};

#endif
