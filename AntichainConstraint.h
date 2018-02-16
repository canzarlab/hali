#ifndef ANTICHAIN_CONSTRAINT_H
#define ANTICHAIN_CONSTRAINT_H

#include "Constraint.h"
#include <mutex>
#include <queue>

class AntichainConstraint : Constraint
{
public:
    AntichainConstraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp);
    int AddTriplets(int nr_rows);

private:
    void RunParallel();
    void AntichainJob(int id);
    double MaxFlow(vi& D, vvd& R);
    double Push(int x, double flow, vvd& R, vi& D);
    double Reset(vn& P, vvd& R);
    void Antichain(vn& P, vvd& R);

    vvi& G;
    int ncr, nr_rows, S, T, Z, SZ;
    mutex qmutex;
    vb B;
    int pi;
};

#endif
