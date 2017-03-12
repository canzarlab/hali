#include "PhylogeneticTree.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "geno/genoNLP.hpp"
#include "geno/augmentedLagrangian.hpp"

Scalar const INF = std::numeric_limits<Scalar>::infinity();

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
        Vector::MapType(cu, _m) = _b;
        Vector::MapType(cl, _m) = Vector::Constant(_m, -INF);
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
    typedef vector<vector<double> > DPT;
public:
    LP(const char* p1, const char* p2) : t1(p1), t2(p2)
    {
        DP.resize(t1.GetNumNodes());
        for (auto& v : DP)
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

        vector<vector<double> > C;
        string line;
        for (int i = 0, k = 0; getline(SimFile, line) && !line.empty(); ++i)
        {
            istringstream ss(line);
            if (!t1.NodeExists(i + 1))
                continue;

            t1.N[i] = k++;        
            C.push_back(vector<double>());
            double w;
            for (int j = 0, l = 0; ss >> w; ++j)
            {
                if (!t2.NodeExists(j + 1))
                    continue;

                t2.N[j] = l++;
                C.back().push_back(w);
            }
        }

        typedef Eigen::Triplet<double> ET;
        vector<ET> Triplets;
        int n = t1.GetNumNodes(), m = t2.GetNumNodes(), k = 0, nr_rows = 0, nr_cols = 0;
        K.resize(n * m);
        Vector c(n * m);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                if (C[i][j] != 0)
                {
                    int col = i * m + j - k;
                    K[i * m + j] = col;
                    Triplets.push_back(ET(i, col, 1.));
                    Triplets.push_back(ET(n + j, col, 1.));
                    
                    if (col > nr_cols) nr_cols = col;
                    if (n + j > nr_rows) nr_rows = n + j;
                    
                    c(col) = -C[i][j];
                }
                else
                    ++k;
            }
        }
        nr_rows++;
        nr_cols++;
        c.conservativeResize(nr_cols);
        
        SpMat A(nr_rows, nr_cols);
        A.setFromTriplets(Triplets.begin(), Triplets.end());
        Vector b = Vector::Ones(nr_rows);        
        SimpleJRF simpleJRF(A, b, c);
        AugmentedLagrangian solver(simpleJRF, 15);
        solver.setParameter("verbose", false);
        solver.setParameter("pgtol", 1e-1); // should influence running time a lot
        solver.setParameter("constraintsTol", 1e-3); 
        solver.solve();
//        x = solver.x();
        x = Vector::ConstMapType(solver.x(), nr_rows);
    }

    void CrossingConstraints()
    {
        dfs1(t1.GetRoot());
    }
        
private:
    double GetWeight(newick_node* nodel, newick_node* noder)
    {
        int i1 = t1.GetIndex(nodel);
        int i2 = t2.GetIndex(noder);
        return x(K[i1 * t2.GetNumNodes() + i2]);
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
        return DP[t1.GetIndex(nodel)][t2.GetIndex(node)] = mx + GetWeight(nodel, node);
    }

    DPT DP;
    Vector x;
    vector<int> K;
    PhylogeneticTree t1, t2;
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
