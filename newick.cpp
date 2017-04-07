#include "newick.h"
#include <fstream>
#include <iostream>
#include <iomanip>

newick_child::newick_child(newick_node* node) : node(node), next(nullptr)
{
}

newick_child::~newick_child()
{
    delete node;
    delete next;
}

newick_node::newick_node(const string& taxon, float dist, newick_child* child) : child(child), taxon(stoi(taxon)), dist(dist), parent(nullptr)
{
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

void print_tree(newick_node* root)
{
    cout << fixed << setprecision(6);
    if (root->child)
    {
        cout << "(";
        for (newick_child* child = root->child; child; child = child->next)
        {
            print_tree(child->node);
            if (child->next)
                cout << ",";
        }
        cout << ")" << root->taxon << ":" << root->dist;
    }
    else
        cout << root->taxon << ":" << root->dist;
}
