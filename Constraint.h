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
    inline int GetCol(int i, int j) const;
    inline int GetCol(newick_node* nodel, newick_node* noder) const;
    inline double GetWeight(newick_node* nodel, newick_node* noder) const;
    inline double GetWeight(int i, int j) const;
    void AddConstraint(int row, vii& P);

    vector<ET>& Triplets;
    Graph &t1, &t2;
    vvi& K;
    Vector& x;
    bool swp;
};

inline int Constraint::GetCol(int i, int j) const
{
    return swp ? K[j][i] : K[i][j];
}

inline int Constraint::GetCol(newick_node* nodel, newick_node* noder) const
{
    return GetCol(nodel->taxoni, noder->taxoni);
}

inline double Constraint::GetWeight(newick_node* nodel, newick_node* noder) const
{
    return GetWeight(nodel->taxoni, noder->taxoni);
}

inline double Constraint::GetWeight(int i, int j) const
{
    int in = GetCol(i, j);
    return in == -1 ? 0 : x(in);
}

#endif
