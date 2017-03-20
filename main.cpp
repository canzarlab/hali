#include "PhylogeneticTree.h"
#include "Timer.h" 
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "geno/genoNLP.hpp"
#include "geno/augmentedLagrangian.hpp"

Scalar const INF = numeric_limits<Scalar>::infinity();
#define EPS 1e-5

// declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<Scalar> SpMat;
typedef PhylogeneticTree Tree;
typedef Eigen::Triplet<double> ET;
typedef vector<vector<double> > DPT;


class SimpleJRF : public GenoNLP
{
 public:
    SimpleJRF(const SpMat& A,
              const Vector& b,
              const Vector& c,
              Vector& x,
              Vector& y)
        : _A(A), _b(b), _c(c), _x(x), _y(y), _n(A.cols()), _m(A.rows())
     {
     }
                                                

    virtual bool getInfo(Index& n, Index& m)
    {
        // number of variables
        n = _n;

        // number of constraints (only real constraints, bound constraints do not count)
        m = _m;
        return true;
    }

    // bounds on the variables
    virtual bool getBounds(Scalar* lb, Scalar* ub)
    {
        Vector::MapType(lb, _n) = Vector::Constant(_n, 0.0);
        Vector::MapType(ub, _n) = Vector::Constant(_n, INF);
        return true;
    }
 
    virtual bool getBoundsConstraints(Scalar* cl, Scalar* cu) 
    {
        // we have equality constraints here
        Vector::MapType(cl, _m) = _b;
        Vector::MapType(cu, _m) = Vector::Constant(_m, INF);
        return true;
    };
  
    virtual bool getStartingPoint(Scalar* x)
    {
//        Vector::MapType(x, _n) = Vector::Zero(_n);
        Vector::MapType(x, _n) = _x;
        return true;
    }

    virtual bool getStartingPointDual(Scalar* y) 
    {
//        Vector::MapType(y, _m) = Vector::Zero(_m);
        Vector::MapType(y, _m) = _y;
        return true;
    };

    virtual bool functionValueAndGradient(const Scalar* variablesPtr,
	                                      Scalar& functionValue,
	                                      Scalar* gradientPtr)
    {
        Vector::ConstMapType x = Vector::ConstMapType(variablesPtr, _n);
        Vector::MapType g = Vector::MapType(gradientPtr, _n);
        functionValue = _c.dot(x);
        g = _c;
        return true;
    }

    virtual bool functionValueConstraints(const Scalar* variablesPtr,
                                          Scalar* constraintValuesPtr)
    {
        Vector::ConstMapType x = Vector::ConstMapType(variablesPtr, _n);
        Vector::MapType constraintValues = Vector::MapType(constraintValuesPtr, _m);
        constraintValues = _A * x;
        return true;
    }
  
    virtual bool gradientConstraintsTimesVector(const Scalar* variablesPtr,
                                                const Scalar* dualVariablesPtr,
                                                Scalar* gradientPtr)
    {
        Vector::ConstMapType x = Vector::ConstMapType(variablesPtr, _n);
        Vector::ConstMapType y = Vector::ConstMapType(dualVariablesPtr, _m);
        Vector::MapType gradient = Vector::MapType(gradientPtr, _n);
        gradient = _A.transpose() * y;
        return true;
    }
  
private:
    const SpMat& _A;
    const Vector& _b;
    const Vector& _c;
    Vector& _x;
    Vector& _y;
    Index _n;
    Index _m;
};

class CrossingConstraint
{
public:
    CrossingConstraint(Tree& t1, Tree& t2, vector<vector<int> >& K, Vector& x, bool swp) : t1(t1), t2(t2), K(K), x(x), swp(swp)
    {
        DP.resize(t1.GetNumNodes());
        for (auto& v : DP)
            v.resize(t2.GetNumNodes());
    }

    int AddTriplets(vector<ET>& Triplets, int nr_rows)
    {
        int ncr = 0;
        DFSLeft(t1.GetRoot());
        for (auto node : t1.L)
        {
            vector<pair<newick_node*, newick_node*> > P;
            Reconstruct(P, node, t2.GetRoot());

            double sum = 0;
            for (auto k : P)
                sum += GetWeight(k.first, k.second);

            if (sum - EPS > 1.0)
            {
                for (auto k : P)
                {
                    int col = GetCol(k.first, k.second);
                    if (col != -1)
                        Triplets.push_back(ET(nr_rows + ncr, col, 1.));
                }
                ncr++;
            }
        }
        return ncr;
    }

    int GetCol(int i, int j)
    {
        return swp ? K[j][i] : K[i][j];
    }

    int GetCol(newick_node* nodel, newick_node* noder)
    {
        return GetCol(t1.GetIndex(nodel), t2.GetIndex(noder));
    }

    double& GetDP(newick_node* nodel, newick_node* noder)
    {
        int i1 = t1.GetIndex(nodel);
        int i2 = t2.GetIndex(noder);
        return DP[i1][i2];
    }

    double GetWeight(newick_node* nodel, newick_node* noder)
    {
        int in = GetCol(nodel, noder);
        return in == -1 ? 0 : -x(in);
    }

    double GetWeightParent(newick_node* nodel, newick_node* noder)
    {
        if (newick_node* parent = nodel->parent)
            return GetDP(parent, noder);
        return 0;
    }

    pair<newick_node*, double> GetMaxChild(newick_node* nodel, newick_node* noder)
    {
        double mx = 0;
        newick_node* mc = nullptr;
        for (newick_child* child = noder->child; child; child = child->next)
        {
            double cw = GetDP(nodel, child->node);
            if (cw > mx)
            {
                mx = cw;
                mc = child->node;
            }
        }
        return make_pair(mc, mx);
    }

    void DFSLeft(newick_node* node)
    {
        DFSRight(t2.GetRoot(), node);
        for (newick_child* child = node->child; child; child = child->next)
            DFSLeft(child->node);
    }

    double DFSRight(newick_node* node, newick_node* nodel)
    {
        double mx = 0;
        for (newick_child* child = node->child; child; child = child->next)
            mx = max(mx, DFSRight(child->node, nodel));
        mx = max(mx, GetWeightParent(nodel, node));
        return GetDP(nodel, node) = mx + GetWeight(nodel, node);
    }

    void Reconstruct(vector<pair<newick_node*, newick_node*> >& P, newick_node* nodel, newick_node* noder)
    {
        double pw = GetWeightParent(nodel, noder), cw;
        newick_node* child;
        tie(child, cw) = GetMaxChild(nodel, noder);
        P.emplace_back(nodel, noder);
        if (nodel->parent && (!child || pw > cw))
            Reconstruct(P, nodel->parent, noder);
        else if (child && (!nodel->parent || cw >= pw))
            Reconstruct(P, nodel, child);
        else
            assert(!nodel->parent && !child);
    }

private:
    DPT DP;
    Tree& t1;
    Tree& t2;
    vector<vector<int> >& K;
    Vector& x;
    bool swp;
};

class LP
{
public:
    LP(const char* p1, const char* p2) : t1(p1), t2(p2), c(t1.GetNumNodes() * t2.GetNumNodes()), nr_rows(0), nr_cols(0)
    {
        K.resize(t1.GetNumNodes());
        for (auto& v : K)
            v.resize(t2.GetNumNodes());
    }

    ~LP()
    {
    }

    void MatchingConstraints(const char* path)
    {
        ifstream SimFile(path);
        if (!SimFile)
        {
            cout << "Failed to open " << path << endl;
            exit(EXIT_FAILURE);
        }

        int n = t1.GetNumNodes(), m = t2.GetNumNodes(), cnt = 0;
        string line;
        for (int i = 0, k = 0; getline(SimFile, line) && !line.empty(); ++i)
        {
            istringstream ss(line);
            if (!t1.NodeExists(i + 1))
                continue;

            t1.N[i + 1] = k;
            double w;
            for (int j = 0, l = 0; ss >> w; ++j)
            {
                if (!t2.NodeExists(j + 1))
                    continue;

                t2.N[j + 1] = l;
                if (w != 0)
                {
                    int col = k * m + l - cnt;
                    K[k][l] = col;
                    Triplets.push_back(ET(k, col, 1.));
                    Triplets.push_back(ET(n + l, col, 1.));
                    c(col) = w;
                }
                else
                {
                    K[k][l] = -1;
                    ++cnt;
                }
                ++l;
            }
            ++k;
        }
        nr_rows = n + m;
        nr_cols = n * m - cnt;
        c.conservativeResize(nr_cols);
        x = Vector::Zero(nr_cols);
        y = Vector::Zero(nr_rows);
        Solve();
    }

    void CrossingConstraints()
    {
//        cout << nr_rows << endl;
        CrossingConstraint cc12(t1, t2, K, x, false);
        nr_rows += cc12.AddTriplets(Triplets, nr_rows);
//        cout << nr_rows << endl;
        CrossingConstraint cc21(t2, t1, K, x, true);
        nr_rows += cc21.AddTriplets(Triplets, nr_rows);
//        cout << "Total number of rows: " << nr_rows << endl;
        
        y.conservativeResizeLike(Vector::Zero(nr_rows)); // resizes y with 0's, but keeping old values intact.         
        Solve();
    }

    void Solve()
    {
        SpMat A(nr_rows, nr_cols);
        A.setFromTriplets(Triplets.begin(), Triplets.end());
        SpMat A_t = A.transpose();
        Vector b = Vector::Ones(nr_rows);
//        x = Vector::Zero(nr_cols);
//        y = Vector::Zero(nr_cols);
        SimpleJRF simpleJRF(A_t, c, b, y, x);
        AugmentedLagrangian solver(simpleJRF, 15);
        solver.setParameter("verbose", false);
        solver.setParameter("pgtol", 1e-1); // should influence running time a lot
        solver.setParameter("constraintsTol", 1e-3);
        Timer timeGeno;
        timeGeno.start();
        solver.solve();
        timeGeno.stop();
        cout << "nr_rows = " << nr_rows << " and nr_cols = " << nr_cols << endl;
        cout << "f = " << solver.f() << " computed in time: " << timeGeno.secs() << " secs" << endl;
        x = Vector::ConstMapType(solver.y(), nr_cols);
        y = Vector::ConstMapType(solver.x(), nr_rows);
    }

private:
    vector<ET> Triplets;
    Vector x;
    Vector y;
    vector<vector<int> > K;
    Tree t1, t2;
    Vector c;
    int nr_rows, nr_cols;
};

int main(int argc, char** argv)
{
    if (argc != 4)
    {
        cout << "usage: " << argv[0] << " <filename.newick> <filename.newick> <filename.sim>" << endl;
        return EXIT_FAILURE;
    }

    LP lp(argv[1], argv[2]);
    lp.MatchingConstraints(argv[3]);
    for (int i = 0; i < 1000; ++i)
        lp.CrossingConstraints();
}
