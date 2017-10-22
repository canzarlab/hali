#include "BnG.h"
#include <iostream>

BnG::BnG(Graph& t1, Graph& t2, string d, double k, bool dag) : LP(t1, t2, d, k, dag)
{
}

void BnG::Solve(string filename)
{
	MatchingConstraints();

	sys_lo.conservativeResizeLike(Vector::Zero(nr_cols));
	sys_hi.conservativeResizeLike(Vector::Ones(nr_cols));
	sys_x.resize(nr_cols);	
	sys_x.assign(nr_cols, false);
	
	Geno();
	SolveLP();

	for (size_t i = 0; i < x.size(); ++i)
		x(i) = round(x(i));
	WriteSolution(filename);
}

void BnG::Cleanup(size_t nr_t, size_t nr_r)
{
    Triplets.resize(nr_t);	
	nr_rows = nr_r;
}

bool BnG::SolveLP()
{
    size_t pos = x.size();
	double val = 0;

	int qwe = 0;

	for (size_t i = 0; i < x.size(); ++i)
		if (x(i) > 0.05 && x(i) < 0.95 && !sys_x[i])
		{    
			if (c(i) > val) 
			{
				pos = i;
				val = c(i);				
			} 
			qwe++;
		}

	if (pos < x.size())
	{	
		sys_x[pos] = true;

		size_t nr_t = Triplets.size();
		size_t nr_r = nr_rows;

		float f1, f2;
		Vector y;

		sys_lo(pos) = 1;
		sys_hi(pos) = 1;
		f1 = Geno();
		y = x;		
		Cleanup(nr_t, nr_r);			

		sys_lo(pos) = 0;
		sys_hi(pos) = 0;
		f2 = Geno();

		cout << max(f1, f2) << ' ' << qwe << endl;

		if (f1 < f2)
		{
			sys_lo(pos) = 1;
			sys_hi(pos) = 1;
			Cleanup(nr_t, nr_r); 
			x = y;
		}

		SolveLP();					
	}

	return 1;
}

float BnG::Geno()
{
	while(1)
	{
	    SpMat A(nr_rows, nr_cols);
	    A.setFromTriplets(Triplets.begin(), Triplets.end());
		Vector b = Vector::Ones(nr_rows);   	
		x = Vector::Zero(nr_cols);
		y = Vector::Zero(nr_rows);
		Vector d = -c;	

	    BranchingJRF simpleJRF(A, b, d, x, y);
		simpleJRF.lo = sys_lo;
		simpleJRF.hi = sys_hi;
	    AugmentedLagrangian solver(simpleJRF, 15);
	    solver.setParameter("verbose", false);
	    solver.setParameter("pgtol", 1e-1); 
	    solver.setParameter("constraintsTol", 1e-4);
	    solver.solve();	 

		x = Vector::ConstMapType(solver.x(), nr_cols);

		if (LP::cf == 1 && Add<CrossingConstraint>())
			continue;
		else if (LP::cf == 2 && (Add<CrossingConstraint>() + Add<IndependentSetConstraint>()))
			continue;			

		return solver.f();	
	}
}
