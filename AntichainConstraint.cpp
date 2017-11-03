#include "Graph.h"
#include "AntichainConstraint.h"
#include <thread>

AntichainConstraint::AntichainConstraint(Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp) : Constraint(t1, t2, K, x, swp), G(((LDAG*)&t2)->G)
{
    Z = t2.GetNumNodes();
    SZ = Z * 2 + 2;
    S = SZ - 2;
    T = SZ - 1;
}

int AntichainConstraint::AddTriplets(vector<ET>& Triplets, int nr_rows)
{
    ncr = 0;
    this->nr_rows = nr_rows;
    this->Triplets = &Triplets;
    vn P;
    for (newick_node* leaf : t1.L)
        DFSLeft(leaf, P);
    RunParallel();
    return ncr;
}

void AntichainConstraint::DFSLeft(newick_node* node, vn& P)
{
    P.push_back(node);
    if (!node->parent)
        PQ.push(P);
    for (newick_parent* parent = node->parent; parent; parent = parent->next)
        DFSLeft(parent->node, P);
    P.pop_back();
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
    while (true)
    {
        vn P;
        {
            lock_guard<mutex> g(qmutex);
            if (PQ.empty())
                break;
            P = move(PQ.front());
            PQ.pop();
        }
        Antichain(P, ((LDAG&)t2).R[id]);
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

void AntichainConstraint::Antichain(vn& P, vvd& R)
{
    vi Q(SZ, -1);
    if (Reset(P, R) - MaxFlow(Q, R) <= 1 + EPS)
        return;

    vii PN;
    for (int i = 0; i < Z; ++i)
        if (Q[i] != -1 && Q[i + Z] == -1)
            for (newick_node* nodel : P)
                PN.emplace_back(nodel->taxoni, i);

    lock_guard<mutex> g(cmutex);
    AddConstraint(*Triplets, nr_rows + ncr, PN);
    ++ncr;
}
