#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "Geno.h"
#include "Graph.h"
#include <utility>
using namespace std;

#define EPS 1e-2

typedef vector<pair<int, int> > vii;

class Constraint
{
public:
    Constraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp);

    virtual int AddTriplets(int nr_rows) = 0;

protected:
    int GetCol(int i, int j) const;
    int GetCol(newick_node* nodel, newick_node* noder) const;
    double GetWeight(newick_node* nodel, newick_node* noder) const;
    double GetWeight(int i, int j) const;
    void AddConstraint(int row, vii& P);

    vector<ET>& Triplets;
    Graph &t1, &t2;
    vvi& K;
    Vector& x;
    bool swp;
};

#endif
