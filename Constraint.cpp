#include "Constraint.h"

Constraint::Constraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp) : Triplets(Triplets), t1(t1), t2(t2), K(K), x(x), swp(swp)
{
}

void Constraint::AddConstraint(int row, vii& P)
{
    for (auto k : P)
    {
        int col = GetCol(k.first, k.second);
        if (col != -1)
            Triplets.emplace_back(row, col, 1.);
    }
}
