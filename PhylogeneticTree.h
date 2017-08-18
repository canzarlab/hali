#ifndef PHYLO_TREE_H
#define PHYLO_TREE_H

#include <vector>
#include <map>
#include <list>
#include "newick.h"
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<double> vd;
typedef vector<vd> vvd;
typedef vector<bool> vb;

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

    vvi G;
    vvd R;
private:
    void BuildNetwork(newick_node* node, newick_node* rnode, vb& C);
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
