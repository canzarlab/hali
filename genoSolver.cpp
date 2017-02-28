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

 Hence, A here is transposed to our original maximization system matrix, b is c and c is b. 
*/

Scalar const INF = std::numeric_limits<Scalar>::infinity();

// declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<Scalar> SpMat;

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
    Vector::MapType(cu, _m) = Vector::Constant(_m, INF);
    Vector::MapType(cl, _m) = _b;
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

/*
void runExample(int m, int n) {
  // generate some random data
  Matrix A = Matrix::Random(m, n).cwiseAbs();
  Vector xOrig = Vector::Random(n);
  for (int i = 0; i < n; ++i)
    if (xOrig[i] < 0)
      xOrig[i] = 0;
		  
//  Vector b = A*xOrig;
  Vector b = Vector::Random(m).cwiseAbs();
  Vector c = Vector::Ones(n);

  SimpleJRF simpleNLP(A, b, c);
  AugmentedLagrangian solver(simpleNLP, 15);
  //solver.setParameter("constraintsTol", .5);

  solver.solve();
  Vector x = Vector::ConstMapType(solver.x(), n);
  double f = solver.f();

  //  std::cout << "xOrig = \n" << xOrig.transpose() << std::endl;
  //  std::cout << "x = \n" << x.transpose() << std::endl;

  std::cout << "norm(A*xOrig - b)^2 = " << (A*xOrig - b).squaredNorm() << std::endl;
  std::cout << "norm(A*x - b)^2     = " << (A*x - b).squaredNorm() << std::endl;
//  cout << "TRALALAL" << A*x -b << endl;
  std::cout << "minimal x = " << x.minCoeff() << std::endl;
  //  std::cout << "x = " << x.transpose() << std::endl;
  std::cout << "function value = " << f << std::endl;
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
    int nr_rows = 0;
    int nr_cols = 0;
    Vector c(n*m); // n*m is pessimistic. c should be resized at the end. 
    
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            if (C[i * n + j] != 0)
            {
//              solver->add_entry(i, m * i + j - k, 1. / C[i * n + j]);
                tripletList.push_back(T(i, m * i + j - k, 1.));
//              solver->add_entry(n + j, i * m + j - k, 1. / C[i * n + j]);
                tripletList.push_back(T(n + j, i * m + j - k, 1.));
                
                if (m*i+j-k > nr_cols)	nr_cols = m*i+j-k;
                if (n+j>nr_rows ) nr_rows = n+j;
                
                c(m*i+j-k) = C[i*n+j];
            }
            else
                ++k;
        }
    }
    
    nr_rows++;
    nr_cols++;
//  cout << "Vector c = " << c << endl;  
    c.conservativeResize(nr_cols);
    
//    cout << "nr_rows = " << nr_rows << " and nr_cols = " << nr_cols << endl;
//    cout << "c = " << c << endl;
    
    

    delete[] C;

    SpMat A(nr_rows, nr_cols); 
    
    cout << "number of rows and columns: " << A.rows() << "   " << A.cols()<< endl;
    
    A.setFromTriplets(tripletList.begin(), tripletList.end()); 
    
//    Matrix A_dense = Matrix(A);

    
    Vector b = Vector::Ones(nr_rows);
    assert (A.cols() == c.rows());
    assert (A.rows() == b.rows());
   
   
  //  cout << "Transposed Matrix A = " << A_dense.transpose() << endl;
  //  cout << "Vector c = " << c << endl;
  //  cout << "Vector b = " << b << endl;


//    runExample(5000,500);
    SpMat A_t = A.transpose(); 
    SimpleJRF simpleJRF(A_t, c, b);  
    AugmentedLagrangian solver(simpleJRF, 15);
       
    solver.solve();

    Vector x = Vector::ConstMapType(solver.x(), nr_rows);
    double f = solver.f();
             
//  std::cout << "x = \n" << x.transpose() << std::endl;
                 
    std::cout << "norm(A.transpose()*x - b)^2     = " << (A_t*x - c).squaredNorm() << std::endl;
    std::cout << "minimal x = " << x.minCoeff() << std::endl;
//  std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "function value = " << f << std::endl;    
   
 
}


