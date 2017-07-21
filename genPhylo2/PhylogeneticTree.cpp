#include "PhylogeneticTree.h"

#include <vector>
#include <stack>
#include "newick.h"
#include <lemon/lgf_writer.h>
#include <lemon/bfs.h>
#include <algorithm> // set_union

using namespace std;

void tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back((str.substr(lastPos, pos - lastPos)).c_str());
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}


void PhylogeneticTree::parseNewick(string filename)
{
//  string line; stack<ListDigraph::Node> S;
//  getline(is, line);
//  vector<string> tokens;
//  tokenize(line, tokens, " \t\",  ;");
//  for (int i = 0; i < tokens.size(); ++i)
//  {
//    cout << tokens[i] << endl;
//    int j = 0; while (tokens[i][j++] == '(')
//    {
//      S.push(_t.addNode());
//    }
//    vector<string> tokentokens;
//    tokenize(tokens[i], tokentokens, "(:"); // should be two parts
//    for (j = 0; j < tokentokens.size(); ++j)
//      cout << "\ttokentoken " << tokentokens[j] << endl;
//    
//  }

	newick_node* root = load_tree(filename.c_str());
    _makeTree(root, INVALID);
    _makeClades();
    delete root;
}

void PhylogeneticTree::_makeTree(newick_node *root, ListDigraph::Node p)
{
  //if (root->taxon) cout << "adding new node (" << root->taxon << " " << root->dist << ")" << endl;
  ListDigraph::Node r = _t.addNode();
  if (!root->taxon.empty()) _label[r] = root->taxon;
  if (root->dist) _dist[r] = root->dist; else root->dist = 0.0;
  if (p != INVALID) _t.addArc(p, r); else _r = r;
  
  if (!root->child) return;
  
  // else
  newick_child *child = root->child;
  while (child != NULL)
  {
      _makeTree(child->node, r);
      child = child->next;
  }
}

void PhylogeneticTree::_makeClades()
{
  list<ListDigraph::Node> L;
  Bfs<ListDigraph> bfs(_t);
  bfs.init();
  bfs.addSource(_r);
  while (!bfs.emptyQueue()) {
    L.push_front(bfs.processNextNode());
  }
  OutDegMap<ListDigraph> out_deg(_t);
  InDegMap<ListDigraph> in_deg(_t);
  
  //cout << "|t| = " << countNodes(_t) << "  |L| = " << L.size() << endl;

  for (list<ListDigraph::Node>::const_iterator cit = L.begin(); cit != L.end(); ++cit)
  {
    ListDigraph::Node i = *cit;
    if (out_deg[i] == 0) _clade[i].push_back(_label[i]);
    else for (ListDigraph::OutArcIt a(_t, i); a != INVALID; ++a)
    {
      for (list<string>::const_iterator cit = _clade[_t.target(a)].begin(); cit != _clade[_t.target(a)].end(); ++cit) _clade[i].push_back(*cit);
    }
    //  cout << _t.id(i) << "\t" << _clade[i].size() << endl;
    _clade[i].sort();
    
   
  }
  
  for (list<ListDigraph::Node>::const_iterator cit = L.begin(); cit != L.end(); ++cit)
    if (out_deg[*cit] > 0 && in_deg[*cit] > 0) _internal_nodes.push_front(*cit);    
  
//  bool binary = false;
//  for (list<ListDigraph::Node>::const_iterator cit = L.begin(); cit != L.end(); ++cit)
//    if (out_deg[*cit] > 2 || in_deg[*cit] > 2) binary = true;
//  if (!binary) cout << "." << flush;
//  
  
  // sort and add clades
//  for (ListDigraph::NodeIt i(_t); i != INVALID; ++i)
//  {
//    _clade[i].sort();
//    list<string> C;
//    for (list<ListDigraph::Node>::const_iterator cit = _clade[i].begin(); cit != _clade[i].end(); ++cit)
//      C.push_back(_label[*cit]);
//    if (out_deg[i] > 0 && in_deg[i] > 0) // non-trivial clades
//      _clades.push_back(C);
//  }
  
  
}

ostream& operator <<(ostream& o, PhylogeneticTree &t) {
  digraphWriter(t.get_graph(), o).nodeMap("labels", t.getLabels())
   .nodeMap("distance", t.getDistance()).attribute("caption", "PhylogeneticTree (JRF)").run();
  
  return o;
}

