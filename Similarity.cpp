#include "Similarity.h"
#include <algorithm>
#include <cmath>
#include <vector>

double JaccardSim(const ls& L1, const ls& L2, double k)
{
    vector<string> I(min(L1.size(), L2.size()));
    auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
    I.resize(iit - I.begin());
    double i = I.size(), u = L1.size() + L2.size() - I.size();
    return 2.0 * pow(i / u, k);
}

double SymdifSim(const ls& L1, const ls& L2)
{
    vector<string> I(min(L1.size(), L2.size()));
    auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
    I.resize(iit - I.begin());
    return 2.0 * I.size();
}
