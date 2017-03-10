#ifndef NEWICK_H
#define NEWICK_H

#include <string>
using namespace std;

struct newick_node;

struct newick_child
{
    newick_child(newick_node* node);
    ~newick_child();

    newick_node* node;
    newick_child* next;
};

struct newick_node
{
    newick_node(const string& taxon, float dist, newick_child* child);
    ~newick_node();

    newick_child* child;
    string taxon;
    float dist;
    newick_node* parent;
};

newick_node* load_tree(const char* filename);
void print_tree(newick_node* root);

#endif

