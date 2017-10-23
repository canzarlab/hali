#ifndef GREEDY_H
#define GREEDY_H

#include "Solver.h"
#include <tuple>

typedef tuple<int, int, double> iid;
typedef vector<iid> viid;

class Greedy : public Solver
{
public:
    Greedy(Graph& t1, Graph& t2, string d, double k, bool dag);

    virtual void Solve(string filename) override;
    void WriteSolution(string fileName) override;
	float GetSolution();
private:
    bool CC(const iid& a, const iid& b);

    vvd A;
    viid E, M;
};

#endif
