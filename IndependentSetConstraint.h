#ifndef INDEPENDENT_SET_CONSTRAINT_H
#define INDEPENDENT_SET_CONSTRAINT_H

#include "Constraint.h"

typedef list<newick_node*> LN;
typedef pair<double, LN> dLN;

class IndependentSetConstraint : Constraint
{
public:
    IndependentSetConstraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp);

    int AddTriplets(int nr_rows);
private:
    double PathSum(newick_node* nodel, newick_node* noder) const;
    void DFSRight(newick_node* noder);
    void DFSLeft(newick_node* nodel, newick_node* noder, double w);
    dLN DFSRight(newick_node* nodel, newick_node* noder);

    vvd D;
};

#endif
