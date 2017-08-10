#include "newick.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

newick_child::newick_child(newick_node* node) : node(node), next(nullptr)
{
}

newick_child::~newick_child()
{
    delete node;
    delete next;
}

newick_node::newick_node(const string& taxon, float dist, newick_child* child) : child(child), taxon(taxon), dist(dist), parent(nullptr)
{
    try { taxoni = stoi(taxon); } catch(...) { taxoni = -1; }
}

newick_node::~newick_node()
{
    delete child;
}

static newick_node* parse_node(string& str, size_t& itr, newick_child* child)
{
    size_t index = str.find_first_of(",);\r\n\0", itr);
    string data = str.substr(itr, index - itr);
    itr = index;

    size_t index2 = data.find_first_of(':');
    string taxon = data.substr(0, index2);
    float dist = 0;
    if (index2 != string::npos && data.back() != ':')
        dist = stof(data.substr(index2 + 1));

    return new newick_node(taxon, dist, child);
}

static newick_node* parse_tree(string& str, size_t& itr)
{
    newick_child* child = nullptr;
    newick_child** childptr = &child;
    for (itr += 1; itr < str.size() && str[itr - 1] != ')'; ++itr)
    {
        newick_node* new_child = nullptr;
        if (str[itr] == '(')
            new_child = parse_tree(str, itr);
        else
            new_child = parse_node(str, itr, nullptr);

        *childptr = new newick_child(new_child);
        childptr = &(*childptr)->next;
    }
    return parse_node(str, itr, child);
}

newick_node* load_tree(const char* filename)
{
    ifstream File(filename);
    if (!File)
        return nullptr;

    int Size;
    File.seekg(0, ios::end);
    Size = File.tellg();
    File.seekg(0, ios::beg);

    string str(Size, 0);
    File.read(&str[0], Size);

    size_t itr = 0;
    return parse_tree(str, itr);
}

typedef map<string, vector<string> > msvs;

newick_node* load_dag_internal(const string& r, msn& M, msvs& C)
{
    if (newick_node* node = M[r])
        return node;

    newick_child* child = nullptr;
    newick_child** childptr = &child;
    for (const string& i : C[r])
    {
        *childptr = new newick_child(load_dag_internal(i, M, C));
        childptr = &(*childptr)->next;
    }
    return M[r] = new newick_node("", 0, child);
}

newick_node* load_dag(const char* f1, bool y, msn& M)
{
    msvs C, P;
    ifstream ef(f1);
    string n1, n2, s;
    while (ef >> n1 >> n2)
    {
        if (y) ef >> s;
        P[n1].push_back(n2);
        P[n2];
        C[n2].push_back(n1);
    }
    for (auto& i : P)
        if (i.second.empty())
            return load_dag_internal(i.first, M, C);
    return nullptr;
}

void print_tree(newick_node* root, ostream& file)
{
    file << fixed << setprecision(6);
    if (root->child)
    {
        file << "(";
        for (newick_child* child = root->child; child; child = child->next)
        {
            print_tree(child->node, file);
            if (child->next)
                file << ",";
        }
        file << ")" << root->taxon << ":" /*<< root->dist*/;
    }
    else
        file << root->taxon << ":" /*<< root->dist*/;
}
