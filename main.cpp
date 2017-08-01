#include "PhylogeneticTree.h"
#include "Timer.h" 
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "geno/genoNLP.hpp"
#include "geno/augmentedLagrangian.hpp"

Scalar const INF = numeric_limits<Scalar>::infinity();
#define EPS 1e-3

// declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<Scalar> SpMat;
typedef PhylogeneticTree Tree;
typedef Eigen::Triplet<double> ET;
typedef vector<vector<double> > DPT;
typedef vector<pair<newick_node*, newick_node*> > VPN;

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
 
/*    virtual bool getBoundsConstraints(Scalar* cl, Scalar* cu) 
    {
        // we have equality constraints here
        Vector::MapType(cl, _m) = _b;
        Vector::MapType(cu, _m) = Vector::Constant(_m, INF);
        return true;
    };
*/  
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
  
protected:
    const SpMat& _A;
    const Vector& _b;
    const Vector& _c;
    Vector& _x;
    Vector& _y;
    Index _n;
    Index _m;
};

class CoveringJRF : public SimpleJRF
{
 public:
    CoveringJRF(const SpMat& A, 
                const Vector& b, 
                const Vector& c, 
                Vector& x, 
                Vector& y): SimpleJRF(A, b, c, x, y) { }
    bool getBoundsConstraints(Scalar* cl, Scalar* cu) override
    {
    // we have equality constraints here
        Vector::MapType(cl, _m) = _b;
        Vector::MapType(cu, _m) = Vector::Constant(_m, INF);
        return true;
    } 
     
};


class PackingJRF : public SimpleJRF
{
 public:
    PackingJRF(const SpMat& A, 
                const Vector& b, 
                const Vector& c, 
                Vector& x, 
                Vector& y): SimpleJRF(A, b, c, x, y) { }
    bool getBoundsConstraints(Scalar* cl, Scalar* cu) override
    {
    // we have equality constraints here
        Vector::MapType(cu, _m) = _b;
        Vector::MapType(cl, _m) = Vector::Constant(_m, -INF);
        return true;
    } 
};


/*
 min c^tx
     x_i(1-x_i) = 0, forall i=1,...,n
     Ax <= b
     x > 0
*/
class IntegerPackingJRF : public SimpleJRF {
 public:
     IntegerPackingJRF(const SpMat& A, 
                       const Vector& b,
                       const Vector& c,
                       Vector& x, 
                       Vector& y) : SimpleJRF(A,b,c,x,y){ _m += _n; }
                       
     bool getBoundsConstraints(Scalar* cl, Scalar* cu) override
     {
         Vector::MapType cLower = Vector::MapType(cl, _m);
         Vector::MapType cUpper = Vector::MapType(cu, _m);
         cLower.head(_n) = Vector::Zero(_n);
         cUpper.head(_n) = Vector::Zero(_n);
         cUpper.tail(_A.rows()) = _b;
         cLower.tail(_A.rows()) = Vector::Constant(_A.rows(), -INF);    
         return true;
     }
     bool functionValueConstraints(const Scalar* variablesPtr,
                                   Scalar* constraintValuesPtr) override
     {
          Vector::ConstMapType x = Vector::ConstMapType(variablesPtr, _n);
          Vector::MapType constraintValues = Vector::MapType(constraintValuesPtr, _m);
          Vector t = Vector::Ones(_n);
          constraintValues.head(_n) = x.cwiseProduct(x-t);
          constraintValues.tail(_A.rows()) = _A * x;
          return true;
      }
      bool gradientConstraintsTimesVector(const Scalar* variablesPtr,
                                          const Scalar* dualVariablesPtr,
                                          Scalar* gradientPtr) override
      {
           Vector::ConstMapType x = Vector::ConstMapType(variablesPtr, _n);
           Vector::ConstMapType y = Vector::ConstMapType(dualVariablesPtr, _m);
           Vector::MapType gradient = Vector::MapType(gradientPtr, _n);
           gradient = (2*x).cwiseProduct(y.head(_n)) - y.head(_n);
           gradient += _A.transpose() * y.tail(_A.rows());
           return true;
       }
      
};                                   


class Constraint
{
public:
    Constraint(Tree& t1, Tree& t2, vector<vector<int> >& K, Vector& x, bool swp) : t1(t1), t2(t2), K(K), x(x), swp(swp)
    {
    }

    virtual int AddTriplets(vector<ET>& Triplets, int nr_rows) = 0;

protected:
    int GetCol(int i, int j)
    {
        return swp ? K[j][i] : K[i][j];
    }

    int GetCol(newick_node* nodel, newick_node* noder)
    {
        return GetCol(nodel->taxoni, noder->taxoni);
    }

    double GetWeight(newick_node* nodel, newick_node* noder)
    {
        int in = GetCol(nodel, noder);
        return in == -1 ? 0 : x(in);
    }

    bool AddConstraint(vector<ET>& Triplets, int row, VPN& P, double sum)
    {
        if (sum - EPS <= 1.0)
            return false;
        
        for (auto k : P)
        {
            int col = GetCol(k.first, k.second);
            if (col != -1)
                Triplets.push_back(ET(row, col, 1.));
        }
        return true;
    }

    Tree& t1;
    Tree& t2;
    vector<vector<int> >& K;
    Vector& x;
    bool swp;
};

class CrossingConstraint : Constraint
{
public:
    CrossingConstraint(Tree& t1, Tree& t2, vector<vector<int> >& K, Vector& x, bool swp) : Constraint(t1, t2, K, x, swp)
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
            VPN P;
            Reconstruct(P, node, t2.GetRoot());

            double sum = 0;
            for (auto k : P)
                sum += GetWeight(k.first, k.second);

            if (AddConstraint(Triplets, nr_rows + ncr, P, sum))
                ncr++;
        }
        return ncr;
    }

private:
    double& GetDP(newick_node* nodel, newick_node* noder)
    {
        return DP[nodel->taxoni][noder->taxoni];
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
            if (cw >= mx)
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

    void Reconstruct(VPN& P, newick_node* nodel, newick_node* noder)
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

    DPT DP;
};

class IndependentSetConstraint : Constraint
{
    typedef list<newick_node*> LN;
    typedef pair<double, LN> dLN;
public:
    IndependentSetConstraint(Tree& t1, Tree& t2, vector<vector<int> >& K, Vector& x, bool swp) : Constraint(t1, t2, K, x, swp)
    {
    }

    int AddTriplets(vector<ET>& Triplets, int nr_rows)
    {
        int ncr = 0;
        for (auto node : t1.L)
        {
            dLN L = DFSRight(node, t2.GetRoot());
            VPN P;
            for (newick_node* noder : L.second)
                for (newick_node* nodel = node; nodel; nodel = nodel->parent)
                    P.emplace_back(nodel, noder);

            if (AddConstraint(Triplets, nr_rows + ncr, P, L.first))
                ++ncr;
        }
        return ncr;
    }

private:
    double PathSum(newick_node* nodel, newick_node* noder)
    {
        return nodel ? GetWeight(nodel, noder) + PathSum(nodel->parent, noder) : 0;
    }
    
    dLN DFSRight(newick_node* nodel, newick_node* noder)
    {
        double w = PathSum(nodel, noder), sum = 0;
        LN V;
        
        for (newick_child* child = noder->child; child; child = child->next)
        {
            double ww;
            LN T;
            tie(ww, T) = DFSRight(nodel, child->node);
            sum += ww;
            V.splice(V.begin(), T);
        }
        
        if (sum > w)
            return make_pair(sum, V);
        return make_pair(w, LN(1, noder));
    }
};

class LP
{
public:
    LP(const char* p1, const char* p2, string d, double k) : d(d), t1(p1), t2(p2), c(t1.GetNumNodes() * t2.GetNumNodes()), nr_rows(0), nr_cols(0), k(k)
    {
        K.resize(t1.GetNumNodes());
        for (auto& v : K)
            v.resize(t2.GetNumNodes());
    }

    ~LP()
    {
    }
    
    void MatchingConstraints()
    {
        cnt = 0;
        DFSLeft(t1.GetRoot());
        int n = t1.GetNumNodes(), m = t2.GetNumNodes();
        nr_rows = n + m;
        nr_cols = n * m - cnt;
        c.conservativeResize(nr_cols);
        // set warm_x and warm_y initialy to zeroes
        warm_x = Vector::Zero(nr_cols);
        warm_y = Vector::Zero(nr_rows);
    }

    int CrossingConstraints()
    {
        int row_old = nr_rows;
        CrossingConstraint cc12(t1, t2, K, x, false);
        nr_rows += cc12.AddTriplets(Triplets, nr_rows);
        CrossingConstraint cc21(t2, t1, K, x, true);
        nr_rows += cc21.AddTriplets(Triplets, nr_rows);
        return nr_rows - row_old;
    }

    int IndependentSetConstraints()
    {
        int row_old = nr_rows;
        IndependentSetConstraint isc12(t1, t2, K, x, false);
        nr_rows += isc12.AddTriplets(Triplets, nr_rows);
        IndependentSetConstraint isc21(t2, t1, K, x, true);
        nr_rows += isc21.AddTriplets(Triplets, nr_rows);
        return nr_rows - row_old;
    }
    
    void Solve()
    {
        clog << "nr_rows = " << nr_rows << " and nr_cols = " << nr_cols << endl;        
        warm_y.conservativeResizeLike(Vector::Zero(nr_rows)); // resizes y with 0's, but keeping old values intact.
        
        SpMat A(nr_rows, nr_cols);
        A.setFromTriplets(Triplets.begin(), Triplets.end());
        SpMat A_t = A.transpose();
        Vector b = Vector::Ones(nr_rows);

        x = Vector::Zero(nr_cols);
        y = Vector::Zero(nr_rows);

//        CoveringJRF simpleJRF(A_t, c, b, warm_y, warm_x);
        
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
       
        /* when Packing: x->x, y->y, when Covering: -y -> x, x -> y */
        //x = -Vector::ConstMapType(solver.y(), nr_cols);
        //y = Vector::ConstMapType(solver.x(), nr_rows);

        x = Vector::ConstMapType(solver.x(), nr_cols);
        y = Vector::ConstMapType(solver.y(), nr_rows);
/*
        warm_x = Vector::ConstMapType(solver.y(), nr_cols);
        warm_y = Vector::ConstMapType(solver.x(), nr_rows);
        
        
        Vector t = A_t*y-c;
        int nr_tight_constr =  nr_cols - (t.array() > 0.1).count();
        clog << "Number of tight constraints in the dual: " << nr_tight_constr << endl;
        //idea, truncate matrix A for columns that correspond to non-tight constraints in dual
        SpMat truncA(nr_rows, nr_tight_constr);
        Vector truncc(nr_tight_constr);
        Vector truncx(nr_tight_constr);
        int truncA_col = 0;
        for (int i=0; i<nr_cols; i++)
            if (t(i)> 0.1)
                x(i) = 0.0;
            else {
                for (SpMat::InnerIterator it(A,i); it; ++it)
                {
                    truncA.coeffRef(it.row(), truncA_col) = it.value();
                }
                truncc(truncA_col) = - c(i);
                truncx(truncA_col) = x(i);
                truncA_col ++;
            }
            
        assert(truncA_col == nr_tight_constr);

        clog << "Truncated matrix formed ... resolve" << endl;
        Vector expy = Vector::Zero(truncA.rows() + truncA.cols());
        IntegerPackingJRF simpleJRF1(truncA, b, truncc, truncx, expy);
        //PackingJRF simpleJRF1(truncA, b, truncc, truncx, y);
        AugmentedLagrangian solver1(simpleJRF1, 15);
        solver1.setParameter("verbose", false);
        solver1.setParameter("pgtol", 1e-1); // should influence running time a lot
        solver1.setParameter("constraintsTol", 1e-5);
        Timer timeGeno1;
        timeGeno1.start();
        solver1.solve();
        timeGeno1.stop();
        clog << "trunc f = " << solver1.f() << " computed in time: " << timeGeno1.secs() << " secs" << endl;
        
        // map the solution back to vector x       
        truncx = Vector::ConstMapType(solver1.x(), nr_tight_constr);
        clog <<"MAX INTEGER PACKING VALUE: " <<  truncx.maxCoeff() << endl;
        truncA_col = 0;
        for (int i=0; i<nr_cols; i++)
            if (t(i)<= 0.1){
                x(i) = truncx(truncA_col);
                truncA_col++;
            }
            
            assert(truncA_col == nr_tight_constr);
*/
    }

    int GetMax(newick_node* node, int& hmax)
    {
        int sum = 0;
        for (newick_child* child = node->child; child; child = child->next)
            sum += GetMax(child->node, hmax);
        hmax += sum;
        return node->child ? sum : 1;
    }

    void WriteSolution(string fileName)
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
        sol_file.close();
        cout << ((d == "j") ? JaccardDist(weight) : SymdifDist(weight)) << " ";
    }

private:
    void DFSLeft(newick_node* node)
    {
        DFSRight(node, t2.GetRoot());
        for (newick_child* child = node->child; child; child = child->next)
            DFSLeft(child->node);
    }

    void DFSRight(newick_node* nodel, newick_node* noder)
    {
        int i = nodel->taxoni, j = noder->taxoni;
        if (nodel->parent && nodel->child && noder->parent && noder->child)
        {
            double w = 0;
            if (d == "j")
                w = JaccardSim(t1.clade[nodel], t2.clade[noder]);
            else
                w = SymdifSim(t1.clade[nodel], t2.clade[noder]);

            if (w != 0)
            {
                int n = t1.GetNumNodes(), m = t2.GetNumNodes();
                int col = i * m + j - cnt;
                K[i][j] = col;
                Triplets.push_back(ET(i, col, 1.));
                Triplets.push_back(ET(n + j, col, 1.));
                c(col) = w;
            }
            else
            {
                K[i][j] = -1;
                ++cnt;
            }
        }
        else
        {
            K[i][j] = -1;
            ++cnt;
        }

        for (newick_child* child = noder->child; child; child = child->next)
            DFSRight(nodel, child->node);
    }

    int cnt;

    double JaccardSim(const list<string>& L1, const list<string>& L2)
    {
        vector<string> I(min(L1.size(), L2.size()));
        auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
        I.resize(iit - I.begin());
        double i = I.size(), u = L1.size() + L2.size() - I.size();
        return 2 * pow(i / u, k);
    }

    double SymdifSim(const list<string>& L1, const list<string>& L2)
    {
        vector<string> I(min(L1.size(), L2.size()));
        auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
        I.resize(iit - I.begin());
        return 2 * I.size();
    }

    float SymdifDist(float weight)
    {
        int max1 = 0, max2 = 0;
        int r1 = GetMax(t1.GetRoot(), max1);
        int r2 = GetMax(t2.GetRoot(), max2);
        return max1 + max2 - r1 - r2 - weight;
    }

    float JaccardDist(float weight)
    {
        return t1.GetNumNodes() - t1.L.size() - 1 + t2.GetNumNodes() - 1 - t2.L.size() - weight;
    }

    double k;
    string d;
    vector<ET> Triplets;
    Vector x;
    Vector y;
    // backup x->warm_x and y->warm_y for two consecutive iterations 
    Vector warm_x;
    Vector warm_y;
    vector<vector<int> > K;
    Tree t1, t2;
    Vector c;
    int nr_rows, nr_cols;
};

int main(int argc, char** argv)
{
    if (argc < 3 || argc > 6)
    {
        clog << "usage: " << argv[0] << " <filename.newick> <filename.newick> [c] [d] [k]" << endl;
        return EXIT_FAILURE;
    }
    int c = (argc < 4) ? 2 : stoi(argv[3]);
    assert(c >= 0 && c <= 2);
    string d = (argc < 5) ? "j" : argv[4];
    assert(d == "j" || d == "s");
    double k = (argc < 6) ? 1 : stod(argv[5]);

    clog << "Comparing trees " << argv[1] << " " << argv[2] << endl;

    LP lp(argv[1], argv[2], d, k);
    lp.MatchingConstraints();
    int cnt = 1, i;
    Timer T;
    T.start();
    for (i = 0; cnt; i++)
    {        
        Timer T_lp, T_cross, T_indep;
        T_lp.start();
        lp.Solve();
        T_lp.stop();
        clog << ">>> Time for solve: \t\t" << T_lp.secs() << " secs" << endl;
        if (c == 0)
            break;

        T_cross.start();
        cnt = lp.CrossingConstraints();
        T_cross.stop();
        clog << ">>> Time for crossing constraints: \t\t" << T_cross.secs() << " secs" << endl;

        if (c == 2)
        {
            T_indep.start();
            cnt += lp.IndependentSetConstraints();
            T_indep.stop();
            clog << ">>> Time for independent set constraints: \t\t" << T_indep.secs() << " secs" << endl;
        }
        clog << "Added " << cnt << " rows." << endl;        
    }
    T.stop();
    clog << "Done comparing trees " << argv[1] << " " << argv[2] << endl;
    clog << "TOTAL TIME : \t\t" << T.secs() << " secs" << endl;
    clog << "Total number of iterations: " <<  i + 1 << endl;
    
    lp.WriteSolution("yeastnet_precollapse.solution");
}
