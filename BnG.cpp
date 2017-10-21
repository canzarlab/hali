#include "BnG.h"
#include <iostream>

BnG::BnG(Graph& t1, Graph& t2, string d, double k, bool dag) : LP(t1, t2, d, k, dag)
{
}

void BnG::Solve(string filename)
{
	MatchingConstraints();
	sys_x.resize(nr_cols);				
	Geno();
	SolveLP();
	for (size_t i = 0; i < x.size(); ++i)
		x(i) = round(x(i));
	WriteSolution(filename);
}

void BnG::Cleanup(size_t nr_t, size_t nr_r)
{
    Triplets.resize(nr_t);
	sys_b.conservativeResizeLike(Vector::Ones(nr_r));		
	nr_rows = nr_r;
}

bool BnG::SolveLP()
{
    size_t pos = x.size();
	double val = 0;

	for (size_t i = 0; i < x.size(); ++i)
		if (x(i) > 1e-3 && x(i) < 1 - 1e-3 && !sys_x[i])
		    if (c(i) > val) 
			{
				pos = i;
				val = c(i);				
			}

	if (pos < x.size())
	{	
		sys_x[pos] = 1;

		size_t nr_t = Triplets.size();
		size_t nr_r = nr_rows;

		Vector y;
		float f1, f2;

		Setup(pos, 1); 
		f1 = Geno();
		y = x;		
		Cleanup(nr_t, nr_r);			

		Setup(pos, 0);
		f2 = Geno();

		if (f1 < f2)
		{
			Cleanup(nr_t, nr_r);
			Setup(pos, 1); 
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
		sys_b.conservativeResizeLike(Vector::Ones(nr_rows));
		//y.conservativeResizeLike(Vector::Zero(nr_rows));	    	
		x = Vector::Zero(nr_cols);
		y = Vector::Zero(nr_rows);
		Vector d = -c;	

	    PackingJRF simpleJRF(A, sys_b, d, x, y);
	    AugmentedLagrangian solver(simpleJRF, 15);
	    solver.setParameter("verbose", false);
	    solver.setParameter("pgtol", 1e-1); 
	    solver.setParameter("constraintsTol", 1e-4);
	    solver.solve();	 

		x = Vector::ConstMapType(solver.x(), nr_cols);

		if (LP::cf == 0) 
			return solver.f();
		else if (LP::cf == 1 && Add<CrossingConstraint>())
			continue;
		else if (LP::cf == 2 && (Add<CrossingConstraint>() + Add<IndependentSetConstraint>()))
			continue;			

		return solver.f();	
	}
}

void BnG::Setup(size_t pos, bool b)
{
	Triplets.push_back(ET(nr_rows++, pos, 1.0));
	if (b)
	{
		Triplets.push_back(ET(nr_rows++, pos, -1.0));
		sys_b.conservativeResizeLike(Vector::Ones(nr_rows));						
		sys_b(nr_rows - 1) = -1.0;
	}
	else
		sys_b.conservativeResizeLike(Vector::Zero(nr_rows));
}
