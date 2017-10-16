#ifndef LPINT_H
#define LPINT_H

#include "LP.h"

class LPInt : public LP
{
public:
    LPInt(Graph& t1, Graph& t2, string d, double k, bool dag);

    virtual void Solve(string filename) override;
private:
    bool SolveLP() override;
};

#endif
