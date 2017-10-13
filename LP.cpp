#include "LP.h"
#include "Timer.h"
#include <iostream>
#include <fstream>

LP::LP(Graph& t1, Graph& t2, string d, double k, bool dag) : Solver(t1, t2, d, k, dag), c(t1.GetNumNodes() * t2.GetNumNodes()), nr_rows(0), nr_cols(0)
{
    K.resize(t1.GetNumNodes(), vi(t2.GetNumNodes(), -1));
}

LP::~LP()
{
}

void LP::MatchingConstraints()
{
    cnt = 0;
    vb P(t1.GetNumNodes());
    int n = t1.GetNumNodes(), m = t2.GetNumNodes();
    DFSLeft(t1.GetRoot(), P, [&](newick_node* nodel, newick_node* noder, double w)
    {
        if (w != 0 && (dag || nodel->parent) && nodel->child && (dag || noder->parent) && noder->child)
        {
            int i = nodel->taxoni, j = noder->taxoni;
            int col = i * m + j - cnt;
            K[i][j] = col;
            Triplets.emplace_back(i, col, 1.);
            Triplets.emplace_back(n + j, col, 1.);
            c(col) = w;
        }
        else
            ++cnt;
    });
    nr_rows = n + m;
    nr_cols = n * m - cnt;
    c.conservativeResize(nr_cols);
}

void LP::Solve()
{
    MatchingConstraints();
    int cnt = 1;
    for (int i = 0; cnt; i++)
    {
        Timer T_lp, T_cross, T_indep;
        T_lp.start();
        SolveLP();
        T_lp.stop();
        clog << ">>> Time for solve: \t\t" << T_lp.secs() << " secs" << endl;
        if (cf == 0)
            break;

        T_cross.start();
        cnt = Add<1>();
        T_cross.stop();
        clog << ">>> Time for crossing constraints: \t\t" << T_cross.secs() << " secs" << endl;

        if (cf == 2)
        {
            T_indep.start();
            cnt += Add<2>();
            T_indep.stop();
            clog << ">>> Time for independent set constraints: \t\t" << T_indep.secs() << " secs" << endl;
        }
        clog << "Added " << cnt << " rows." << endl;
    }
}

bool LP::SolveLP()
{
    clog << "nr_rows = " << nr_rows << " and nr_cols = " << nr_cols << endl;

    SpMat A(nr_rows, nr_cols);
    A.setFromTriplets(Triplets.begin(), Triplets.end());
    SpMat A_t = A.transpose();
    Vector b = Vector::Ones(nr_rows);

    x = Vector::Zero(nr_cols);
    y = Vector::Zero(nr_rows);

    Vector c1 = -c;
    PackingJRF simpleJRF(A, b, c1, x, y);
    AugmentedLagrangian solver(simpleJRF, 15);
    solver.setParameter("verbose", false);
    solver.setParameter("pgtol", 1e-1); // should influence running time a lot
    solver.setParameter("constraintsTol", 1e-4);
    Timer timeGeno;
    timeGeno.start();
    solver.solve();
    timeGeno.stop();

    clog << "f = " << solver.f() << " computed in time: " << timeGeno.secs() << " secs" << endl;

    x = Vector::ConstMapType(solver.x(), nr_cols);
    y = Vector::ConstMapType(solver.y(), nr_rows);
    return true;
}

int LP::GetMax(newick_node* node, int& hmax) const
{
    int sum = 0;
    for (newick_child* child = node->child; child; child = child->next)
        sum += GetMax(child->node, hmax);
    hmax += sum;
    return node->child ? sum : 1;
}

void LP::WriteSolution(string fileName)
{
    ofstream sol_file(fileName);
    float weight = 0;
    for (size_t i = 0; i < K.size(); i++)
    {
        for (size_t j = 0; j < K[i].size(); j++)
        {
            if (K[i][j] != -1)
            {
                weight += x(K[i][j]) * c(K[i][j]);
                sol_file << x(K[i][j]) << "\t" ;
            }
            else
                sol_file << 0 << "\t";
        }
        sol_file << endl;
    }
    if (dag)
        cout << weight << " ";
    else
        cout << ((d == "j") ? JaccardDist(weight) : SymdifDist(weight)) << " ";
}

float LP::SymdifDist(float weight) const
{
    int max1 = 0, max2 = 0;
    int r1 = GetMax(t1.GetRoot(), max1);
    int r2 = GetMax(t2.GetRoot(), max2);
    return max1 + max2 - r1 - r2 - weight;
}

float LP::JaccardDist(float weight) const
{
    return t1.GetNumNodes() - t1.L.size() - 1 + t2.GetNumNodes() - 1 - t2.L.size() - weight;
}
