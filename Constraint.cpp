#include "Constraint.h"

Constraint::Constraint(Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp) : t1(t1), t2(t2), K(K), x(x), swp(swp)
{
}


int Constraint::GetCol(int i, int j) const
{
    return swp ? K[j][i] : K[i][j];
}

int Constraint::GetCol(newick_node* nodel, newick_node* noder) const
{
    return GetCol(nodel->taxoni, noder->taxoni);
}

double Constraint::GetWeight(newick_node* nodel, newick_node* noder) const
{
    return GetWeight(nodel->taxoni, noder->taxoni);
}

double Constraint::GetWeight(int i, int j) const
{
    int in = GetCol(i, j);
    return in == -1 ? 0 : x(in);
}

void Constraint::AddConstraint(vector<ET>& Triplets, int row, vii& P)
{
    for (auto k : P)
    {
        int col = GetCol(k.first, k.second);
        if (col != -1)
            Triplets.emplace_back(row, col, 1.);
    }
}
