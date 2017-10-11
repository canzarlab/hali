#ifndef GREEDY_H
#define GREEDY_H

#include "Solver.h"
#include <tuple>

typedef tuple<int, int, int> iii;
typedef vector<iii> viii;

class Greedy : public Solver
{
public:
    Greedy(Graph& t1, Graph& t2, string d, double k, bool dag);

    void Solve();
    void WriteSolution(string fileName);
private:
    bool CC(const iii& a, const iii& b);

    vvd A;
    viii E, M;
};

#endif
