#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <list>
#include <map>
#include "newick.h"
using namespace std;

class tree
{
    int counter;
public:
    map<newick_node*, list<string> > clade;
    newick_node* root;

    tree(const char* filename) : counter(1)
    {
        Init(root = load_tree(filename));
    }

    ~tree()
    {
        delete root;
    }

private:
    void Init(newick_node* node)
    {
        if (!node->child)
            clade[node].push_back(node->taxon);

        node->taxon = to_string(counter++);
        for (newick_child* child = node->child; child; child = child->next)
        {
            child->node->parent = node;
            Init(child->node);
            list<string> &cl = clade[node], &cr = clade[child->node];
            cl.insert(cl.end(), cr.begin(), cr.end());
        }
        clade[node].sort();
    }
};

class generator
{
    string d;
    double k;
    tree t1, t2;
    ofstream sim;
public:
    generator(const char* f1, const char* f2, string d, double k) :
        d(d), k(k), t1(f1), t2(f2), sim("sim_t1_t2")
    {
    }

    void generate()
    {
        DFSLeft(t1.root);
        ofstream t1m("t1mod"), t2m("t2mod");
        print_tree(t1.root, t1m);
        print_tree(t2.root, t2m);
    }

private:
    void DFSLeft(newick_node* node)
    {
        DFSRight(node, t2.root);
        sim << '\n';
        for (newick_child* child = node->child; child; child = child->next)
            DFSLeft(child->node);
    }

    void DFSRight(newick_node* nodel, newick_node* noder)
    {
        double w = 0;
        if (nodel->parent && nodel->child && noder->parent && noder->child)
        {
            if (d == "j")
                w = 2 - jaccard(t1.clade[nodel], t2.clade[noder], k);
            else
                w = t1.clade[nodel].size() + t2.clade[noder].size() - symdif(t1.clade[nodel], t2.clade[noder]);
        }
        sim << w << "\t";

        for (newick_child* child = noder->child; child; child = child->next)
            DFSRight(nodel, noder->child->node);
    }

    double jaccard(const list<string> &L1, const list<string> &L2, double k)
    {
        vector<string> I(min(L1.size(), L2.size()));
        auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
        I.resize(iit - I.begin());
        double i = I.size(), u = L1.size() + L2.size() - I.size();
        return (u == i) ? 0 : 2 - 2 * pow(i / u, k);
    }

    double symdif(const list<string> &L1, const list<string> &L2)
    {
        vector<string> I(min(L1.size(), L2.size()));
        auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
        I.resize(iit - I.begin());
        return L1.size() + L2.size() - 2 * I.size();
    }
};

int main(int argc, char** argv)
{
    if (argc < 4 || argc > 5)
    {
        cout << "usage: " << argv[0] << " <t1> <t2> <d> [k]" << endl;
        return EXIT_FAILURE;
    }
    generator(argv[1], argv[2], argv[3], (argc == 5) ? stod(argv[4]) : 1).generate();
}
