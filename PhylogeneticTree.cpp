#include "PhylogeneticTree.h"
#include <lemon/lgf_writer.h>
#include <lemon/bfs.h>

PhylogeneticTree::PhylogeneticTree(string filename) : _lab(_t), _dst(_t), _clade(_t)
{
    newick_node* root = load_tree(filename.c_str());
    makeTree(root, INVALID);
    delete root;
}

void PhylogeneticTree::makeTree(newick_node* root, ListDigraph::Node p)
{
    ListDigraph::Node r = _t.addNode();
    _n++;
    _lab[r] = root->taxon;
    _map[root->taxon] = _t.id(r);
    _dst[r] = root->dist;

    if (p != INVALID)
        _t.addArc(p, r);
    else
        _r = r;

    if (root->child)
        for (newick_child* child = root->child; child; child = child->next)
            makeTree(child->node, r);
    else
        _l.push_front(r);
}

long PhylogeneticTree::getNodeId(string s) const
{
    auto iter = _map.find(s);
    if (iter != _map.end())
        return iter->second;
    return -1;
}

const ListDigraph::Node PhylogeneticTree::getParent(const ListDigraph::Node& x) const
{
    ListDigraph::InArcIt arc(_t, x);
    if (arc == INVALID) return INVALID;
    return _t.source(arc);
}

ostream& operator<<(ostream& o, PhylogeneticTree &t) 
{
    digraphWriter(t._t, o).nodeMap("name", t._lab).nodeMap("distance", t._dst).attribute("caption", "PhylogeneticTree (JRF)").run();
    return o;
}

void PhylogeneticTree::makeClades()
{
  list<ListDigraph::Node> L;
  Bfs<ListDigraph> bfs(_t);
  bfs.init();
  bfs.addSource(_r);

  while (!bfs.emptyQueue()) 
  {
    L.push_front(bfs.processNextNode());
  }
  OutDegMap<ListDigraph> out_deg(_t);
  InDegMap<ListDigraph> in_deg(_t);

  for (list<ListDigraph::Node>::const_iterator cit = L.begin(); cit != L.end(); ++cit)
  {
    ListDigraph::Node i = *cit;
    if (out_deg[i] == 0) _clade[i].push_back(_lab[i]);
    else for (ListDigraph::OutArcIt a(_t, i); a != INVALID; ++a)
    {
      for (list<string>::const_iterator cit = _clade[_t.target(a)].begin(); cit != _clade[_t.target(a)].end(); ++cit) 
        _clade[i].push_back(*cit);
    }
    _clade[i].sort();
  }
  //for (list<ListDigraph::Node>::const_iterator cit = L.begin(); cit != L.end(); ++cit)
  //  if (out_deg[*cit] > 0 && in_deg[*cit] > 0) _internal_nodes.push_front(*cit);
}

