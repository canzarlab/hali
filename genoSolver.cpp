#include <iostream>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "PhylogeneticTree.h"


#include "geno/genoNLP.hpp"
#include "geno/augmentedLagrangian.hpp"

/*
  solves problem
  
  min_x c'*x
  s.t.  A*x >= b
        x   >= 0

 Hence, A here is transposed to our original maximization system matrix. 
*/

Scalar const INF = std::numeric_limits<Scalar>::infinity();

// declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<Scalar> SpMat; // probably boolean matrix would also work.

class SimpleJRF : public GenoNLP {
public:
  SimpleJRF(const SpMat& A, 
	    const Vector& b,
	    const Vector& c) 
    : _A(A), _b(b), _c(c), _n(A.cols()), _m(A.rows())
  {    
  }

  virtual bool getInfo(Index& n, Index& m) {
    // number of variables
    n = _n;

    // number of constraints (only real constraints, bound constraints do not count)
    m = _m;
    return true;
  }

  // bounds on the variables
  virtual bool getBounds(Scalar* lb, Scalar* ub) {
    Vector::MapType(lb, _n) = Vector::Constant(_n, 0.0);
    Vector::MapType(ub, _n) = Vector::Constant(_n, INF);
    return true;
  }
 
  virtual bool getBoundsConstraints(Scalar* cl, Scalar* cu) 
  {
    // we have equality constraints here
    Vector::MapType(cl, _m) = Vector::Constant(_m, -INF);
    Vector::MapType(cu, _m) = _b;
    return true;
  };
  
  virtual bool getStartingPoint(Scalar* x) {
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
	                                      Scalar* gradientPtr) {
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


/*void runExample(int m, int n) {
  // generate some random data
  Matrix A = Matrix::Random(m, n);
  Vector xOrig = Vector::Random(n);
  for (int i = 0; i < n; ++i)
    if (xOrig[i] < 0)
      xOrig[i] = 0;
		  
  Vector b = A*xOrig;
  Vector c = Vector::Ones(n);

  SimpleNLP simpleNLP(A, b, c);
  AugmentedLagrangian solver(simpleNLP, 15);

  solver.solve();
  Vector x = Vector::ConstMapType(solver.x(), n);
  double f = solver.f();

  //  std::cout << "xOrig = \n" << xOrig.transpose() << std::endl;
  //  std::cout << "x = \n" << x.transpose() << std::endl;

  std::cout << "norm(A*xOrig - b)^2 = " << (A*xOrig - b).squaredNorm() << std::endl;
  std::cout << "norm(A*x - b)^2     = " << (A*x - b).squaredNorm() << std::endl;
  std::cout << "minimal x = " << x.minCoeff() << std::endl;
  //  std::cout << "x = " << x.transpose() << std::endl;
  std::cout << "function value = " << f << std::endl;
}


void printUsage() {
  std::cout << "usage: example3 <m> <n>" << std::endl;
  std::cout << "\n\tsolves min_x c'*x" << std::endl;
  std::cout << "\t       s.t.  A*x <= b" << std::endl;
  std::cout << "\t               x >= 0\n" << std::endl;
  std::cout << "<m> - number of rows of random A" << std::endl;
  std::cout << "<n> - number of cols of random A" << std::endl;
}
*/

int main(int argc, char** argv)
{
    if (argc != 4)
    {
        cout << "usage: " << argv[0] << " <filename.newick> <filename.newick> <filename.sim> " << endl;
        return EXIT_FAILURE;
    }

    PhylogeneticTree t1(argv[1]);
    PhylogeneticTree t2(argv[2]);
    int n = t1.getNumNodes(), m = t2.getNumNodes();

    ifstream SimFile(argv[3]);
    double* C = new double[n * m];
    for (int i = 0; i < n * m; ++i)
        SimFile >> C[i];

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    
    int k = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            if (C[i * n + j] != 0)
            {
                tripletList.push_back(T(m * i + j - k, i, 1));
                tripletList.push_back(T(m * i + j - k, n + j, 1));
//                solver->add_entry(i, m * i + j - k, 1. / C[i * n + j]);
//                solver->add_entry(n + j, i * m + j - k, 1. / C[i * n + j]);
            }
            else
                ++k;
        }
    }
    

    delete[] C;

    SpMat A(n, m); // TODO: set nr_rows and nr_cols
    // TODO: set vectors c and b (in dual). 

    A.setFromTriplets(tripletList.begin(), tripletList.end()); 
    
    

    
}


