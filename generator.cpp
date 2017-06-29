#include "newick.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <algorithm>
using namespace std;

void print(newick_node* root, ofstream& file)
{
	if(root->child)
	{
		file << "(";
		for(newick_child* child = root->child; child; child = child->next)
		{
			print(child->node, file);
			if(child->next)
			file << ",";
		}
		file << ")" << root->taxon << ":";
	}
	else
       file << root->taxon << ":";
}

void yule(int n, const char* filename)
{
	srand(unsigned(time(0)));
	ofstream file(filename);
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
	
	if(file.is_open())
	{
		print(root, file);
		file.close();
	}
	
	delete root;
}

void uniform(int n, const char* filename)
{
	srand(unsigned(time(0)));
	ofstream file(filename);
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
			if(c->node == E[i].second)
			{
				c->node = leaf2;
				break;
			}
			
		E.emplace_back(leaf2, leaf1);
		E.emplace_back(leaf2, E[i].second);
		E[i] = make_pair(E[i].first, leaf2); 
	}

	if(file.is_open())
	{
		print(root, file);
		file.close();
	}
	
	delete root;
}
/*
int main()
{
	int n;
	cin >> n;
	uniform(n, "data/a1");
}
*/