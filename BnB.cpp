#include "BnB.h"
#include "IndependentSetConstraint.h"
#include "AntichainConstraint.h"
#include "CrossingConstraint.h"
#include <iostream>

void BnB::Solve()
{
    MatchingConstraints();
    sys_lb = (double)INF;
    sys_x.resize(nr_cols);
    if (SolveLP()) clog << "TMP: Success." << endl;

    // TODO delete me...
    double sum = 0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        x(i) = round(x(i));
        sum += c(i) * x(i);
    }
    clog << -sys_lb << ' ' << sum << endl;
}

void BnB::Cleanup(size_t nr_t, size_t nr_r)
{
    Triplets.resize(nr_t);
    sys_b.resize(nr_t);
    nr_rows = nr_r;
}

bool BnB::SolveLP()
{
    int nr_t = Triplets.size();
    int nr_r = nr_rows;

    while (1)
    {
        SpMat A(nr_rows, nr_cols);
        A.setFromTriplets(Triplets.begin(), Triplets.end());
        sys_b.conservativeResizeLike(Vector::Ones(nr_rows));
        Vector d = -c;

        x = Vector::Zero(nr_cols);
        y = Vector::Zero(nr_rows);

        PackingJRF simpleJRF(A, sys_b, d, x, y);
        AugmentedLagrangian solver(simpleJRF, 15);
        solver.setParameter("verbose", false);
        solver.setParameter("pgtol", 1e-1);
        solver.setParameter("constraintsTol", 1e-4);
        solver.solve();

        x = Vector::ConstMapType(solver.x(), nr_cols);

        if (cf == 0)
            return 1;
        else if (cf == 1 && Add<CrossingConstraint>())
            continue;
        else if (cf == 2 && (Add<CrossingConstraint>() + Add<IndependentSetConstraint>()))
            continue;


        if (sys_lb == (double)INF)
            sys_lb = solver.f();
        else if (solver.f() <= sys_lb * 1.01)
        {
            clog << "!!!!!!!!!!!!!!! odbacio !!!!!!!!!!!!!!" << endl;
            Cleanup(nr_t, nr_r);
            return 0;
        }

        x = Vector::ConstMapType(solver.x(), nr_cols);

        break;
    }

    size_t pos = x.size();
    double val = 0;

    for (size_t i = 0; i < x.size(); ++i)
        if (x(i) > 1e-3 && x(i) < 1 - 1e-3 && !sys_x[i])
            if (c(i) > val)
            {
                pos = i;
                val = c(i);
            }

    while (pos < x.size())
    {
        if (SolveRec(pos, 1)) break;
        if (SolveRec(pos, 0)) break;
        Cleanup(nr_t, nr_r);
        return 0;
    }

    Cleanup(nr_t, nr_r);
    return 1;
}

bool BnB::SolveRec(size_t pos, bool b)
{
    Triplets.emplace_back(nr_rows++, pos, 1.0);
    if (b)
    {
        Triplets.emplace_back(nr_rows++, pos, -1.0);
        sys_b.conservativeResizeLike(Vector::Ones(nr_rows));
        sys_b(nr_rows - 1) = -1.0;
    }
    else
        sys_b.conservativeResizeLike(Vector::Zero(nr_rows));

    sys_x[pos] = 1;
    bool f = SolveLP();
    sys_x[pos] = 0;

    Cleanup(Triplets.size() - 1 - b, nr_rows - 1 - b);
    return f;
}
