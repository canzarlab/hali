#ifndef LPINT_H
#define LPINT_H

#include "LP.h"

class LPInt : public LP
{
public:
    LPInt(Graph& t1, Graph& t2, string d, double k, bool dag);

    void Solve();
private:
    bool SolveLP();
    // backup x->warm_x and y->warm_y for two consecutive iterations
    Vector warm_x, warm_y;
};

#endif
