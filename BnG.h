#ifndef BNG_H
#define BNG_H

#include "LP.h"

class BnG : public LP
{
public:
    BnG(Graph& t1, Graph& t2, string d, double k, bool dag);
    virtual void Solve(string filename) override;

private:
    void  Cleanup(size_t nr_t, size_t nr_r);
    bool  SolveLP() override;
    float Geno();

    vector<bool> sys_x;
	Vector		 sys_lo;
	Vector		 sys_hi;
};

#endif
