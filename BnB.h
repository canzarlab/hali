#ifndef BNB_H
#define BNB_H

#include "LP.h"
#include "Greedy.h"

class BnB : public LP
{
public:
    BnB(Graph& t1, Graph& t2, string d, double k, bool dag, double c);
    virtual void Solve(string filename) override;

private:
    void Cleanup(size_t nr_t, size_t nr_r);
    bool SolveLP(Vector xp, int depth); // override je bilo
    bool SolveRec(unsigned int k, vector<pair<int, int>>& p, Vector xp, int depth);

	Greedy       G;
    double       sys_lb;
    vector<bool> sys_x;
	Vector       sys_s;
	Vector		 sys_lo;
	Vector		 sys_hi;
	
	double       con_eps;
};

#endif
