#ifndef NEWICK_H
#define NEWICK_H

#include <string>
#include <map>
using namespace std;

struct newick_node;

struct newick_child
{
    newick_child(newick_node* node);
    ~newick_child();

    newick_node* node;
    newick_child* next;
};

typedef newick_child newick_parent;

struct newick_node
{
    newick_node(const string& taxon = "", float dist = 0, newick_child* child = nullptr);
    ~newick_node();

    newick_child* child;
    string taxon;
    int taxoni;
    float dist;
    newick_parent* parent;
};

typedef map<string, newick_node*> msn;

newick_node* load_tree(const char* filename);
newick_node* load_dag(const char* f1, bool y, msn& M);
void dealloc_dag(newick_node* node, int n);
void print_tree(newick_node* root, ostream& file);

#endif

