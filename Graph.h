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

const size_t NR_THREADS = 4;

class Graph
{
public:
    Graph() : _n(0) { }
    virtual ~Graph() { dealloc_dag(root, _n); };

    int GetNumNodes() const { return _n; }
    newick_node* GetRoot() const { return root; }

    vn L;
    mnls clade;

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
    DAG(const char* f1, const char* f2);
    ~DAG() { }

    DAG* BuildNetwork();
protected:
    virtual void Relation(int d, int a) = 0;
private:
    void BuildNetwork(newick_node* node, newick_node* rnode, vvb& C);
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

class LDAG : public DAG
{
public:
    LDAG(const char* f1, const char* f2);

    vvi G;
    vvd R[NR_THREADS];
protected:
    virtual void Relation(int d, int a);
};

class GDAG : public DAG
{
public:
    GDAG(const char* f1, const char* f2);

    vvb D;
protected:
    virtual void Relation(int d, int a);
};


#endif