#ifndef PHYLO_TREE_H
#define PHYLO_TREE_H

#include <iostream>
#include <map>
#include <lemon/list_graph.h>
#include "newick.h"

using namespace lemon;
using namespace std;

class PhylogeneticTree
{
public:
    // Constructor.
    //
    // Arguments: 
    //   filename - path to the tree in newick format. 
    //
    // All the preprocessing will be done here, so we will end up having a static data container. 
    PhylogeneticTree(string filename);

    // Get methods.
    const ListDigraph&             getGraph()                            const { return _t; }
    const ListDigraph::Node&       getRoot()                             const { return _r; }
    const list<ListDigraph::Node>& getLeaves()                           const { return _l; }
    long                           getNumNodes()                         const { return _n; }
    string                         getLabel(ListDigraph::Node& x)        const { return _lab[x]; }
    double                         getDistance(ListDigraph::Node& x)     const { return _dst[x]; }
    long                           getNodeId(string s)                   const;
    const ListDigraph::Node        getParent(const ListDigraph::Node& x) const;

    // Legacy.
    const list<string>&            clade(ListDigraph::Node i) const { return _clade[i]; }
    void                           makeClades();

    // Operators.
    friend ostream& operator<<(ostream &o, PhylogeneticTree &t);

private:
    // Data.
    ListDigraph                         _t;   // Tree.
    ListDigraph::Node                   _r;   // Root.
    list<ListDigraph::Node>             _l;   // Leaves.
    ListDigraph::NodeMap<string>        _lab; // Node labels.
    ListDigraph::NodeMap<double>        _dst; // Node distances. I'm not even sure whether the last two containers will be needed.
    map<string, long>                   _map;
    long                                _n;
    ListDigraph::NodeMap<list<string> > _clade;

    void makeTree(newick_node* root, ListDigraph::Node p);
};

#endif
