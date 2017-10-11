#ifndef INDEPENDENT_SET_CONSTRAINT_H
#define INDEPENDENT_SET_CONSTRAINT_H

#include "Constraint.h"

typedef list<newick_node*> LN;
typedef pair<double, LN> dLN;

class IndependentSetConstraint : Constraint
{
public:
    IndependentSetConstraint(Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp);

    int AddTriplets(vector<ET>& Triplets, int nr_rows);
private:
    double PathSum(newick_node* nodel, newick_node* noder) const;
    dLN DFSRight(newick_node* nodel, newick_node* noder);
};

#endif
