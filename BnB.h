#ifndef BNB_H
#define BNB_H

#include "LP.h"

class BnB : public LP
{
public:
    void Solve();

private:
    void Cleanup(size_t nr_t, size_t nr_r);
    bool SolveLP();
    bool SolveRec(size_t pos, bool b);

    double       sys_lb;
    Vector       sys_b;
    vector<bool> sys_x;
};

#endif
