#ifndef LPINT_H
#define LPINT_H

#include "LP.h"

typedef tuple<int, int> ii;

class LPInt : public LP
{
public:
    LPInt(Graph& t1, Graph& t2, string d, double k, bool dag);

    virtual void Solve(string filename) override;
private:
    bool CC(const ii& a, const ii& b);
    void AddConstraint(const ii& a, const ii& b);
    virtual bool SolveLP() override;
};

#endif
