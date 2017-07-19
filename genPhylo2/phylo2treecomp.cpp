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

#include "PhylogeneticTree.h"

#include <dirent.h>



using namespace std;
using namespace lemon;

string getNewick(ListDigraph::Node r, PhylogeneticTree& t, ListDigraph::NodeMap<int>& nm)
{
	stringstream ss;
	ListDigraph::OutArcIt a(t.get_graph(), r);
	if (a == INVALID) { //leaf 
		ss << nm[r];
		return ss.str();
	}
	
	string result = "(";
	for (; a != INVALID; ++a) {
		result += getNewick(t.get_graph().target(a), t, nm);
		result += ",";
	}
	result.back() = ')';
        ss << result << nm[r];
	return ss.str();
}

bool equal(const list<string> &L1, const list<string> &L2)
{
  if (L1.size() != L2.size()) return false;
  if (L1 == L2) return true;
  return false;
}

double Jaccard_weight(const list<string> &L1, const list<string> &L2, double k)
{
  if (equal(L1, L2)) return 0.0;
  vector<string>::iterator uit, iit;
  vector<string> U(L1.size() + L2.size()), I(min(L1.size(), L2.size()));
  uit = set_union(L1.begin(), L1.end(), L2.begin(), L2.end(), U.begin());
  iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
  U.resize(uit - U.begin());
  I.resize(iit - I.begin());
  double J = I.size()/(double)U.size();
  return 2 - 2*pow(J, k);
}

int main(int argc, char * const argv[]) {
    // Initialize the argument parser
    ArgParser ap(argc, argv);

    string filename_t1, filename_t2, dirname;
    double k = 1.0;

    // Add a string option with storage reference for file name
    ap.refOption("t1", "Name of file that contains first tree in Newick format []", filename_t1, false);
    ap.refOption("t2", "Name of file that contains second tree in Newick format []", filename_t2, false);
    ap.refOption("d",  "Name of directory with trees in Newick format for all-against-all *.tre comparisons []", dirname, false);
    ap.refOption("k", "Exponent of GRF metric [1.0]", k, false);
    // Perform the parsing process
    // (in case of any error it terminates the program)
    ap.parse();

    if (!filename_t1.empty())
    {
        PhylogeneticTree t1, t2;
        t1.parseNewick(filename_t1);
        t2.parseNewick(filename_t2);
        //cout << "t1: " << t1 << endl << endl;
        //cout << "t2: " << t2 << endl;
        cout << t1.internal_nodes().size() << " " << t2.internal_nodes().size() << endl;

        //const ListDigraph::NodeMap<std::string>& labels_t1 = t1.getLabels();
        //const ListDigraph::NodeMap<std::string>& labels_t2 = t2.getLabels();
        //std::stringstream ss;
        //int n_t1 = 0, n_t2 = 0;
        //for (ListDigraph::NodeIt n(t1.get_graph()); n != INVALID; ++n)
        //{
        //cout << labels_t1[n] << ": " << endl;
        //for (list<string>::const_iterator cit = t1.clade(n).begin(); cit != t1.clade(n).end(); ++cit) cout << *cit << endl;
        //}

        ListDigraph::NodeMap<int> new_label_t1(t1.get_graph());
        ListDigraph::NodeMap<int> new_label_t2(t2.get_graph());
        int count = 1;
        for (ListDigraph::NodeIt i(t1.get_graph()); i != INVALID; ++i) {
            new_label_t1[i] = count;
            count++;
        }
        count = 1;
        for (ListDigraph::NodeIt i(t2.get_graph()); i != INVALID; ++i) {
            new_label_t2[i] = count;
            count++;
        }
        std::ofstream sim_matrix;
        sim_matrix.open("sim_t1_t2");
        stringstream ss;
        //for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i)
        //for (list<ListDigraph::Node>::const_iterator j = t2.internal_nodes().begin(); j != t2.internal_nodes().end(); ++j)
        for (ListDigraph::NodeIt i(t1.get_graph()); i != INVALID; ++i) {
            ss.str("");
            ListDigraph::InArcIt a1(t1.get_graph(), i);
            ListDigraph::OutArcIt o1(t1.get_graph(), i);
            for (ListDigraph::NodeIt j(t2.get_graph()); j != INVALID; ++j)
            {
                double w = 0.0;
                ListDigraph::InArcIt a2(t2.get_graph(), j);
                ListDigraph::OutArcIt o2(t2.get_graph(), j);
                if (a1 != INVALID && a2 != INVALID && o1 != INVALID && o2 != INVALID) //t1.clade(i).size() > 1 && t2.clade(j).size() > 1)
                    w = 2 - Jaccard_weight(t1.clade(i), t2.clade(j), k);
                ss << w << "\t";

                //source_in_t1[e] = *i; target_in_t2[e] = *j;
                //            //cout << "RF_weight: " << Jaccard_RF_weight(t1.clade(*i), t2.clade(*j)) << ", JRF_weight: " << Jaccard_weight(t1.clade(*i), t2.clade(*j)) << endl;
                //                    }
            }
            sim_matrix << ss.str() << endl;
        }
        sim_matrix.close();

        std::ofstream tmod;
        tmod.open("t1mod");
        string newick_out = "";
        ListDigraph::Node r = t1.getRoot();
        newick_out = getNewick(r,t1,new_label_t1);
        tmod << newick_out << endl;
        tmod.close();
        tmod.open("t2mod");
        newick_out = "";
        r = t2.getRoot();
        newick_out = getNewick(r,t2,new_label_t2);
        tmod << newick_out << endl;
        tmod.close();


        //    cout << "Robinson-Foulds(" << filename_t1 << ", " << filename_t2 << "):\t" << RobinsonFouldsDistance(t1, t2) << endl;
        //    cout << "Jaccard-Robinson-Foulds(" << filename_t1 << ", " << filename_t2 << "):\t" << JaccardRobinsonFouldsDistance(t1, t2, k) << endl;
        //    std::pair<std::pair<double, double>, double> d = JaccardRobinsonFouldsDistanceCPLEX(t1, t2, k);
        //    cout << "Jaccard-Robinson-Foulds_CPLEX(" << filename_t1 << ", " << filename_t2 << "):\t(" << d.first.first << ", " << d.first.second << ") " << d.second << " %" << endl;
        //    cout << "|t1|: " << t1.internal_nodes().size() << endl;
        //    cout << "|t2|: " << t2.internal_nodes().size() << endl;
        /*
    cout <<filename_t1 << "|" << filename_t2 << "\t" << k << "\t" << flush;
    lemon::Timer clk_RF; double d1 = RobinsonFouldsDistance(t1, t2);
    if (verbosity >= 3) cout << d1 << "\t" << (t1.internal_nodes().size() + t2.internal_nodes().size() - d1)/2 << "\t" << clk_RF.userTime() << "\t" << flush;
    lemon::Timer clk_JRF; std::pair<double, int> d2 = JaccardRobinsonFouldsDistance(t1, t2, k);
    if (verbosity >= 3) cout << d2.first << "\t" << d2.second << "\t" << clk_JRF.userTime() << "\t" << flush;
    lemon::Timer clk_JRF_CPLEX; std::pair<std::pair<double, double>, std::pair<double, double> > d3 = JaccardRobinsonFouldsDistanceCPLEX(t1, t2, k);
      if (verbosity >= 3) cout << d3.first.first << "\t" << d3.first.second << "\t" << d3.second.first << "\t" << d3.second.second << "\t" << clk_JRF_CPLEX.userTime() << flush;
    cout << endl;
*/

        //std::pair<std::pair<double, double>, std::pair<double, double> > d = GeneralTreeDistanceGreedy(t1, t2, filename_sim);
        //std::pair<std::pair<double, double>, std::pair<double, double> > d = JaccardRobinsonFouldsDistanceCPLEX(t1, t2, k);
        //cout << d.first.first << " " << d.first.second << " " << d.second.first << " " << d.second.second << endl;

    }

    if (!dirname.empty())
    {
        DIR*    dir;
        dirent* pdir;
        std::vector<std::string> treefiles;

        dir = opendir(dirname.c_str());

        while ((pdir = readdir(dir))) {
            string fn = pdir->d_name;
            if (fn.substr(fn.find_last_of(".") + 1) == "tre") treefiles.push_back(fn);
        }

        int n = treefiles.size();
        std::vector<PhylogeneticTree*> T(n);
        for (int i = 0; i < n; ++i)
        {
            PhylogeneticTree *t = new PhylogeneticTree();
            t->parseNewick(dirname + "/" + treefiles[i]);
            T[i] = t;
        }
        
        //    for (int i = 0; i < n; ++i) for (int j = i; j < n; ++j)
        //    {
        //      //if (verbosity >= 6) cout << "\r" << (i*n+j+1)/(double)(n*(n+1)/2)*100 << " %                      " << flush;
        //      cout << treefiles[i] << "|" << treefiles[j] << "\t" << flush;
        //      lemon::Timer clk_RF; double d1 = RobinsonFouldsDistance(*T[i], *T[j]);
        //      if (verbosity >= 3) cout << d1 << "\t" << clk_RF.userTime() << "\t" << flush;
        //      lemon::Timer clk_JRF; double d2 = JaccardRobinsonFouldsDistance(*T[i], *T[j], k);
        //      if (verbosity >= 3) cout << d2 << "\t" << clk_JRF.userTime() << "\t" << flush;
        //      lemon::Timer clk_JRF_CPLEX; std::pair<std::pair<double, double>, std::pair<double, double> > d3 = JaccardRobinsonFouldsDistanceCPLEX(*T[i], *T[j], k);
        //      if (verbosity >= 3) cout << d3.first.first << " " << d3.first.second << " " << d3.second.first << "\t" << clk_JRF_CPLEX.userTime() << endl;
        //    }
    }



    return 0;
}
