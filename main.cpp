#include "PhylogeneticTree.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "geno/genoNLP.hpp"
#include "geno/augmentedLagrangian.hpp"

Scalar const INF = numeric_limits<Scalar>::infinity();
#define EPS 1e-9

// declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<Scalar> SpMat;

class SimpleJRF : public GenoNLP
{
public:
    SimpleJRF(const SpMat& A, 
              const Vector& b,
              const Vector& c) 
        : _A(A), _b(b), _c(c), _n(A.cols()), _m(A.rows())
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
        Vector::MapType(x, _n) = Vector::Zero(_n);
        return true;
    }

    virtual bool getStartingPointDual(Scalar* y) 
    {
        Vector::MapType(y, _m) = Vector::Zero(_m);
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
    Index _n;
    Index _m;
};

class LP
{
    typedef PhylogeneticTree Tree;
    typedef Eigen::Triplet<double> ET;
    typedef vector<vector<double> > DPT;
public:
    LP(const char* p1, const char* p2) : t1(p1), t2(p2), ncr(0), c(t1.GetNumNodes() * t2.GetNumNodes())
    {
        DP.resize(t1.GetNumNodes());
        for (auto& v : DP)
            v.resize(t2.GetNumNodes());
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

            t1.N[i] = k;
            double w;
            for (int j = 0, l = 0; ss >> w; ++j)
            {
                if (!t2.NodeExists(j + 1))
                    continue;

                t2.N[j] = l;
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
        int nr_rows = n + m;
        int nr_cols = n * m - cnt;
        c.conservativeResize(nr_cols);
        
        SpMat A(nr_rows, nr_cols);
        A.setFromTriplets(Triplets.begin(), Triplets.end());
        SpMat A_t = A.transpose();
        Vector b = Vector::Ones(nr_rows);        
        SimpleJRF simpleJRF(A_t, c, b);
        AugmentedLagrangian solver(simpleJRF, 15);
        solver.setParameter("verbose", false);
        solver.setParameter("pgtol", 1e-1); // should influence running time a lot
        solver.setParameter("constraintsTol", 1e-3); 
        solver.solve();
        cout << solver.f() << endl;
        x = Vector::ConstMapType(solver.y(), nr_cols);
    }

    void CrossingConstraints()
    {
        dfs1(t1.GetRoot());
        for (auto node : t1.L)
        {
            vector<int> P;
            f(P, node, t2.GetRoot());

            int n = t1.GetNumNodes(), m = t2.GetNumNodes();
            double sum = 0;
            for (auto k : P)
                sum += DP[k / m][k % m];

            if (sum - EPS > 1)
            {
                for (auto k : P)
                    Triplets.push_back(ET(n + m + ncr + 1, k, 1.));
                ncr++;
            }
        }
        cout << ncr << endl;
    }

private:
    int GetCol(newick_node* nodel, newick_node* noder)
    {
        int i1 = t1.GetIndex(nodel);
        int i2 = t2.GetIndex(noder);
        return K[i1][i2];
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

    void dfs1(newick_node* node)
    {
        dfs2(t2.GetRoot(), node);
        for (newick_child* child = node->child; child; child = child->next)
            dfs1(child->node);
    }

    double dfs2(newick_node* node, newick_node* nodel)
    {
        double mx = 0;
        for (newick_child* child = node->child; child; child = child->next)
            mx = max(mx, dfs2(child->node, nodel));
        mx = max(mx, GetWeightParent(nodel, node));
        return GetDP(nodel, node) = mx + GetWeight(nodel, node);
    }

    void f(vector<int>& P, newick_node* nodel, newick_node* noder)
    {        
        double pw = GetWeightParent(nodel, noder), cw;
        newick_node* child;
        tie(child, cw) = GetMaxChild(nodel, noder);
        P.push_back(GetCol(nodel, noder));

        if (nodel->parent && (!child || pw > cw))
            f(P, nodel->parent, noder);
        else if (child && (!nodel->parent || cw >= pw))
            f(P, nodel, child);
        else
            assert(!nodel->parent && !child);
    }

    DPT DP;
    vector<ET> Triplets;
    Vector x;
    vector<vector<int> > K;
    PhylogeneticTree t1, t2;
    int ncr;
    Vector c;
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
    lp.CrossingConstraints();    
}
