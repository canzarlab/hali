#include "newick.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <algorithm>
#include <cassert>
using namespace std;

newick_node* yule(int n)
{
	unsigned random[2 * n - 1], size = 1;
	vector<newick_node*> L;
	
	for(int i = 0; i < (2 * n - 1); i++)
		random[i] = i + 1;
		
	random_shuffle(random, random + 2 * n - 1);
	
	newick_node* root = new newick_node(to_string(random[0]), 0, nullptr);
	L.push_back(root);
	
	while(L.size() < n)
	{
		int i = rand() % L.size();
		newick_node* leaf1 = new newick_node(to_string(random[size++]), 0, nullptr);
		newick_node* leaf2 = new newick_node(to_string(random[size++]), 0, nullptr);
		newick_child* child = new newick_child(leaf1);
		child->next = new newick_child(leaf2);
		
		L[i]->child = child;
		L[i] = leaf1;
		L.push_back(leaf2);
	}
    return root;
}

newick_node* uniform(int n)
{
	unsigned random[2 * n - 1], size = 3;
	vector<pair<newick_node*, newick_node*>> E;
	
	for(int i = 0; i < (2 * n - 1); i++)
		random[i] = i + 1;
		
	random_shuffle(random, random + 2 * n - 2);
	
	newick_node* leaf1 = new newick_node(to_string(random[0]), 0, nullptr);
	newick_node* leaf2 = new newick_node(to_string(random[1]), 0, nullptr);
	newick_child* child = new newick_child(leaf1);
	child->next = new newick_child(leaf2);
	
	newick_node* root = new newick_node(to_string(random[2]), 0, child);
	E.emplace_back(root, leaf1);
	E.emplace_back(root, leaf2);
	
	for(int k = 0; k < n - 2; k++)
	{
		int i = rand() % E.size();
		
		leaf1 = new newick_node(to_string(random[size++]), 0, nullptr);
		leaf2 = new newick_node(to_string(random[size++]), 0, new newick_child(E[i].second));
		leaf2->child->next = new newick_child(leaf1);

		for(newick_child* c = E[i].first->child; c; c = c->next)
        {
			if(c->node == E[i].second)
			{
				c->node = leaf2;
				break;
			}
		}	
		E.emplace_back(leaf2, leaf1);
		E.emplace_back(leaf2, E[i].second);
		E[i] = make_pair(E[i].first, leaf2); 
	}
    return root;
}

newick_node* (*F[])(int) = { uniform, yule };

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        cout << "usage: " << argv[0] << " <n-leaves> <n-trees>\n";
        return 1;
    }
    srand(unsigned(time(0)));
    int n_leaves = stoi(argv[1]);
    int n_trees = stoi(argv[2]);
    for (int i = 0; i < n_trees; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            newick_node* root = F[j](n_leaves);
            ofstream file(string("data") + to_string(j) + string("/a") + to_string(i));
            assert(file);
            print_tree(root, file);
            file << '\n';
            delete root;
        }
    }
}
