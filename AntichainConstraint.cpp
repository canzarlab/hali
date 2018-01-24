#include "Graph.h"
#include "AntichainConstraint.h"
#include <thread>

AntichainConstraint::AntichainConstraint(Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp) : Constraint(t1, t2, K, x, swp), G(((LDAG*)&t2)->G), B(t1.GetNumNodes()), pi(0)
{
    Z = t2.GetNumNodes();
    SZ = Z * 2 + 2;
    S = SZ - 2;
    T = SZ - 1;
}

int AntichainConstraint::AddTriplets(vector<ET>& Triplets, int nr_rows)
{
    int ncr = 0;
    RunParallel();
    for (int i = 0; i < C.size(); ++i)
        AddConstraint(Triplets, nr_rows + ncr++, C[i]);
    return ncr;
}

void AntichainConstraint::RunParallel()
{
    vector<thread> vt;
    for (int i = 0; i < NR_THREADS; ++i)
        vt.emplace_back(&AntichainConstraint::AntichainJob, this, i);

    for (int i = 0; i < NR_THREADS; ++i)
        vt[i].join();
}

void AntichainConstraint::AntichainJob(int id)
{
    LDAG &g1 = (LDAG&)t1, &g2 = (LDAG&)t2;
    while (true)
    {
        int i;
        {
            lock_guard<mutex> g(qmutex);
            if (pi == g1.P.size())
                break;
            i = pi++;
            if (all_of(g1.P[i].begin(), g1.P[i].end(), [&](newick_node* node){return B[node->taxoni];}))
                continue;
        }
        Antichain(i, g1.P[i], g2.R[id]);
    }
}

bool AntichainConstraint::AugmentingPath(vi& Q, vvd& R)
{
    queue<int> W;
    W.push(S);
    fill(Q.begin(), Q.end(), -1);
    while (!W.empty())
    {
        int k = W.front(); W.pop();
        if (k == T)
            return true;

        for (int x : G[k])
            if (Q[x] == -1 && R[k][x] > 0)
                Q[x] = k, W.push(x);
    }
    return false;
}

double AntichainConstraint::MaxFlow(vi& Q, vvd& R)
{
    double flow = 0;
    while (AugmentingPath(Q, R))
    {
        double aug = numeric_limits<double>::infinity();
        for (int x = T; x != S; x = Q[x])
            aug = min(aug, R[Q[x]][x]);

        for (int x = T; x != S; x = Q[x])
            R[Q[x]][x] -= aug, R[x][Q[x]] += aug;

        flow += aug;
    }
    return flow;
}

double AntichainConstraint::Reset(vn& P, vvd& R)
{
    double sum = 0;
    for (int i = 0; i < Z; ++i)
    {
        R[S][i] = R[i][S] = R[T][i + Z] = 0;
        for (int j : G[i + Z])
            R[i + Z][j] = 0;

        for (newick_node* nodel : P)
        {
            R[S][i] += GetWeight(nodel->taxoni, i);
            R[i + Z][T] += GetWeight(nodel->taxoni, i);
        }
        sum += R[S][i];
    }
    return sum;
}

void AntichainConstraint::Antichain(int ci, vn& P, vvd& R)
{
    vi Q(SZ, -1);
    if (Reset(P, R) - MaxFlow(Q, R) <= 1 + EPS)
        return;

    vii PN;
    for (int i = 0; i < Z; ++i)
        if (Q[i] != -1 && Q[i + Z] == -1)
            for (newick_node* nodel : P)
                PN.emplace_back(nodel->taxoni, i);

    lock_guard<mutex> g(qmutex);
    if (!all_of(P.begin(), P.end(), [&](newick_node* node){return B[node->taxoni];}))
        C.push_back(move(PN));
    for (newick_node* node : P)
        B[node->taxoni] = true;
}
