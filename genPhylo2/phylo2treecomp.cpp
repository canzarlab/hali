/*
 * File:   generalizePhylo.cpp
 * Author: canzar
 *
 * Created on 5 February 2015, 11:57 am
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <math.h>

#include <lemon/list_graph.h>
#include <lemon/arg_parser.h>
#include <lemon/bfs.h>

#include "newick.h"

using namespace std;
using namespace lemon;

class PhylogeneticTree
{
public:
    PhylogeneticTree(string filename);
    const ListDigraph &get_graph() const { return _t; }
    const ListDigraph::NodeMap<string> &getLabels() const { return _label; }
    const list<string> &clade(ListDigraph::Node i) const { return _clade[i]; }
    ListDigraph::Node getRoot() {return _r; }

private:
    ListDigraph _t; // the tree structure as directed lemon graph
    ListDigraph::Node _r; // root node
    ListDigraph::NodeMap<list<string> > _clade;
    ListDigraph::NodeMap<string> _label;
    void _makeTree(newick_node *root, ListDigraph::Node p);
    void _makeClades();
};

PhylogeneticTree::PhylogeneticTree(string filename) : _clade(_t), _label(_t)
{
	newick_node* root = load_tree(filename.c_str());
    _makeTree(root, INVALID);
    _makeClades();
    delete root;
}

void PhylogeneticTree::_makeTree(newick_node *root, ListDigraph::Node p)
{
    ListDigraph::Node r = _t.addNode();
    if (!root->taxon.empty()) _label[r] = root->taxon;
    if (p != INVALID) _t.addArc(p, r); else _r = r;

    for (newick_child *child = root->child; child; child = child->next)
        _makeTree(child->node, r);
}

void PhylogeneticTree::_makeClades()
{
    list<ListDigraph::Node> L;
    Bfs<ListDigraph> bfs(_t);
    bfs.init();
    bfs.addSource(_r);
    while (!bfs.emptyQueue())
        L.push_front(bfs.processNextNode());

    OutDegMap<ListDigraph> out_deg(_t);
    InDegMap<ListDigraph> in_deg(_t);

    for (list<ListDigraph::Node>::const_iterator cit = L.begin(); cit != L.end(); ++cit)
    {
        ListDigraph::Node i = *cit;
        if (out_deg[i] == 0) _clade[i].push_back(_label[i]);
        else for (ListDigraph::OutArcIt a(_t, i); a != INVALID; ++a)
             {
                 for (list<string>::const_iterator cit = _clade[_t.target(a)].begin(); cit != _clade[_t.target(a)].end(); ++cit) _clade[i].push_back(*cit);
             }
        _clade[i].sort();
    }
}

string getNewick(ListDigraph::Node r, PhylogeneticTree& t, ListDigraph::NodeMap<int>& nm)
{
	stringstream ss;
	ListDigraph::OutArcIt a(t.get_graph(), r);
	if (a == INVALID) { //leaf 
		ss << nm[r] << ":";
		return ss.str();
	}
	
	string result = "(";
	for (; a != INVALID; ++a) {
		result += getNewick(t.get_graph().target(a), t, nm);
		result += ",";
	}
	result.back() = ')';
	ss << result << nm[r] << ":";
	return ss.str();
}

double Jaccard_weight(const list<string> &L1, const list<string> &L2, double k)
{
    vector<string> I(min(L1.size(), L2.size()));
    auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
    I.resize(iit - I.begin());
    double u = I.size(), i = L1.size() + L2.size() - I.size();
    return (u == i) ? 0 : 2 - 2 * pow(i / u, k);
}

double symdif_weight(const list<string> &L1, const list<string> &L2)
{
    vector<string> I(min(L1.size(), L2.size()));
    auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
    I.resize(iit - I.begin());
    return L1.size() + L2.size() - 2 * I.size();
}

int main(int argc, char * const argv[])
{
    // Initialize the argument parser
    ArgParser ap(argc, argv);

    string filename_t1, filename_t2;
    double k = 1.0;
    string d = "j";

    // Add a string option with storage reference for file name
    ap.refOption("t1", "Name of file that contains first tree in Newick format []", filename_t1, true);
    ap.refOption("t2", "Name of file that contains second tree in Newick format []", filename_t2, true);
    ap.refOption("k", "Exponent of GRF metric [1.0]", k, false);
    ap.refOption("d", "Distance function [j/s]", d, false);
    // Perform the parsing process
    // (in case of any error it terminates the program)
    ap.parse();

    PhylogeneticTree t1(filename_t1), t2(filename_t2);

    ListDigraph::NodeMap<int> new_label_t1(t1.get_graph());
    ListDigraph::NodeMap<int> new_label_t2(t2.get_graph());
    int count = 1;
    for (ListDigraph::NodeIt i(t1.get_graph()); i != INVALID; ++i)
        new_label_t1[i] = count++;

    count = 1;
    for (ListDigraph::NodeIt i(t2.get_graph()); i != INVALID; ++i)
        new_label_t2[i] = count++;

    ofstream sim_matrix("sim_t1_t2");
    stringstream ss;
    for (ListDigraph::NodeIt i(t1.get_graph()); i != INVALID; ++i)
    {
        ss.str("");
        ListDigraph::InArcIt a1(t1.get_graph(), i);
        ListDigraph::OutArcIt o1(t1.get_graph(), i);
        for (ListDigraph::NodeIt j(t2.get_graph()); j != INVALID; ++j)
        {
            double w = 0.0;
            ListDigraph::InArcIt a2(t2.get_graph(), j);
            ListDigraph::OutArcIt o2(t2.get_graph(), j);
            if (a1 != INVALID && a2 != INVALID && o1 != INVALID && o2 != INVALID)
            {
                if (d == "j")
                    w = 2 - Jaccard_weight(t1.clade(i), t2.clade(j), k);
                else
                    w = t1.clade(i).size() + t2.clade(j).size() - symdif_weight(t1.clade(i), t2.clade(j));
            }
            ss << w << "\t";
        }
        sim_matrix << ss.str() << endl;
    }

    ofstream t1mod("t1mod");
    t1mod << getNewick(t1.getRoot(), t1, new_label_t1) << endl;
    ofstream t2mod("t2mod");
    t2mod << getNewick(t2.getRoot(), t2, new_label_t2) << endl;
}
