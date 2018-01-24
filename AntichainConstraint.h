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
    void Antichain(int ci, vn& P, vvd& R);

    vvi& G;
    int S, T, Z, SZ;
    mutex cmutex, qmutex;
    vector<vii> C;
    int pi;
};

#endif
