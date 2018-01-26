#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <map>
#include <list>
#include "newick.h"
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<double> vd;
typedef vector<vd> vvd;
typedef vector<vb> vvb;

// Up to 64 threads supported
const size_t NR_THREADS = 4;

class Graph
{
public:
    Graph() : _n(0) { }
    virtual ~Graph() { dealloc_dag(root, _n); };

    int GetNumNodes() const { return _n; }
    newick_node* GetRoot() const { return root; }
    virtual Graph* Init();

    vn L;
    mnls clade;
    vvb D;
protected:
    void TransitiveClosure();
    void TransitiveClosure(newick_node* node, newick_node* rnode, vvb& C);
    void Init(newick_node* node);
    virtual void Leaf(newick_node* node) { }
    virtual void Child(newick_node* node, newick_node* child) { }
    virtual void Relation(int d, int a) { }
    newick_node* root;
    long _n;
};

class DAG : public Graph
{
public:
    DAG(const char* f1, const char* f2);
    ~DAG() { }

protected:
    virtual void Relation(int d, int a) { }
};

class Tree : public Graph
{
public:
    Tree() { }
    Tree(const char* f1);
    Tree(const char* f1, const char* f2);
    ~Tree() { }

protected:
    virtual void Leaf(newick_node* node) override;
    virtual void Child(newick_node* node, newick_node* child) override;
    bool t;
};

class LDAG : public DAG
{
public:
    LDAG(const char* f1, const char* f2);
    Graph* Init();

    vvi G;
    vvd R[NR_THREADS];
    vector<vn> P;
protected:
    void GenPaths(newick_node* node, vn& P, vector<vn>& Q);
    virtual void Relation(int d, int a) override;
};

#endif
