#ifndef ANTICHAIN_CONSTRAINT_H
#define ANTICHAIN_CONSTRAINT_H

#include "Constraint.h"
#include <mutex>
#include <queue>

class AntichainConstraint : Constraint
{
public:
    AntichainConstraint(Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp);
    int AddTriplets(vector<ET>& Triplets, int nr_rows);

private:
    void RunParallel();
    void AntichainJob(int id);
    bool AugmentingPath(vi& Q, vvd& R);
    double MaxFlow(vi& Q, vvd& R);
    double Reset(vn& P, vvd& R);
    void Antichain(vn& P, vvd& R);

    vvi& G;
    int ncr, nr_rows, S, T, Z, SZ;
    vector<ET>* Triplets;
    mutex cmutex, qmutex;
    int pi;
};

#endif
