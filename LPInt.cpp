#include "LPInt.h"
#include "Timer.h"
#include <iostream>

LPInt::LPInt(Graph& t1, Graph& t2, string d, double k, bool dag) : LP(t1, t2, d, k, dag)
{
}

void LPInt::Solve(string filename)
{
    MatchingConstraints();
    int cnt = 1;
    for (int i = 0; cnt; i++)
    {
        Timer T_lp, T_cross, T_indep;
        T_lp.start();
        SolveLP();
        WriteSolution(filename);
        T_lp.stop();
        clog << ">>> Time for solve: \t\t" << T_lp.secs() << " secs" << endl;
        if (cf == 0)
            break;

        vector<ii> E;
        for (size_t i = 0; i < K.size(); i++)
            for (size_t j = 0; j < K[i].size(); j++)
                if (K[i][j] != -1 && x(K[i][j]) > 1-1e-1)
                    E.emplace_back(i, j);

        cnt = 0;
        for (int i = 0; i < E.size(); ++i)
            for (int j = i + 1; j < E.size(); ++j)
                if (!CC(E[i], E[j]))
                    AddConstraint(E[i], E[j]), cnt++;
        clog << "Added " << cnt << " rows." << endl;
    }
}

bool LPInt::CC(const ii& a, const ii& b)
{
    int i = get<0>(a), j = get<0>(b);
    int x = get<1>(a), y = get<1>(b);
    if (i == j || x == y) return false;
    if (cf == 0) return true;
    bool c2 = (t1.D[j][i] || t1.D[i][j]) == (t2.D[x][y] || t2.D[y][x]);
    bool c1 = (t1.D[i][j] == t2.D[x][y]); // assuming c2 is satisfied
    if (cf == 1) return !c2 || c1;
    return c2 && c1;
}

void LPInt::AddConstraint(const ii& a, const ii& b)
{
    int i = get<0>(a), j = get<1>(a);
    int k = get<0>(b), l = get<1>(b);
    Triplets.emplace_back(nr_rows, K[i][j], 1.);
    Triplets.emplace_back(nr_rows++, K[k][l], 1.);
}

bool LPInt::SolveLP()
{
    clog << "nr_rows = " << nr_rows << " and nr_cols = " << nr_cols << endl;

    SpMat A(nr_rows, nr_cols);
    A.setFromTriplets(Triplets.begin(), Triplets.end());
    SpMat A_t = A.transpose();
    Vector b = Vector::Ones(nr_rows);

    x = warm_x;
    //x = Vector::Zero(nr_cols);
    y = Vector::Zero(nr_rows + nr_cols);

    Vector c1 = -c;
    IntegerPackingJRF simpleJRF(A, b, c1, warm_x, y);
    AugmentedLagrangian solver(simpleJRF, 15);
    solver.setParameter("verbose", false);
    solver.setParameter("pgtol", 1e-1); // should influence running time a lot
    solver.setParameter("constraintsTol", 1e-3);
    Timer timeGeno;
    timeGeno.start();
    solver.solve();
    timeGeno.stop();

    clog << "f = " << solver.f() << " computed in time: " << timeGeno.secs() << " secs" << endl;

    warm_x = x = Vector::ConstMapType(solver.x(), nr_cols);
    y = Vector::ConstMapType(solver.y(), nr_rows);
    return true;
}
