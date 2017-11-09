#ifndef LPFINT_H
#define LPFINT_H

#include "LPInt.h"

class LPFInt : public LPInt
{
public:
    LPFInt(Graph& t1, Graph& t2, string d, double k, bool dag);

    virtual void Solve(string filename) override;
};

#endif
