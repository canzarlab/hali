#include "Solver.h"

Solver::Solver(Graph& t1, Graph& t2, string d, double k, bool dag) : t1(t1), t2(t2), d(d), k(k), dag(dag)
{
}

double Solver::JaccardSim(const list<string>& L1, const list<string>& L2) const
{
    vector<string> I(min(L1.size(), L2.size()));
    auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
    I.resize(iit - I.begin());
    double i = I.size(), u = L1.size() + L2.size() - I.size();
    return 2 * pow(i / u, k);
}

double Solver::SymdifSim(const list<string>& L1, const list<string>& L2) const
{
    vector<string> I(min(L1.size(), L2.size()));
    auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
    I.resize(iit - I.begin());
    return 2 * I.size();
}
