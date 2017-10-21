#include "BnB.h"
#include <iostream>

BnB::BnB(Graph& t1, Graph& t2, string d, double k, bool dag) : LP(t1, t2, d, k, dag)
{
}

void BnB::Solve(string filename)
{
	MatchingConstraints();	
	x = Vector::Zero(nr_cols);
	sys_lb = (double)INF;
	sys_x.resize(nr_cols);
	SolveLP();
	x = sys_s;

	for (size_t i = 0; i < x.size(); ++i)
		x(i) = round(x(i));

    WriteSolution(filename);
}

void BnB::Cleanup(size_t nr_t, size_t nr_r)
{
    Triplets.resize(nr_t);
	sys_b.conservativeResizeLike(Vector::Ones(nr_r));		
	nr_rows = nr_r;
}

bool BnB::SolveLP()
{
    int nr_t = Triplets.size();
	int nr_r = nr_rows;
	double f;

	while(1)
	{
	    SpMat A(nr_rows, nr_cols);
	    A.setFromTriplets(Triplets.begin(), Triplets.end());
		sys_b.conservativeResizeLike(Vector::Ones(nr_rows));		    
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
		{
			sys_s = x;
			return 1;
		}
		else if (LP::cf == 1 && Add<CrossingConstraint>())
			continue;
		else if (LP::cf == 2 && (Add<CrossingConstraint>() + Add<IndependentSetConstraint>()))
			continue;

		f = solver.f();			

		if (sys_lb != double(INF) && f >= sys_lb * 1.0) 
		{				
			Cleanup(nr_t, nr_r);			
			return 0; 
		}			

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

	if (pos < x.size())
		while(1)
	    {	
			bool b1 = (SolveRec(pos, 1));
			bool b2 = (SolveRec(pos, 0));								
			if (b1 || b2) break;			
			Cleanup(nr_t, nr_r);			
			return 0;
		}
	else if (sys_lb > f)
	{
		sys_s = x;
		sys_lb = f; 
	}

	Cleanup(nr_t, nr_r);
	return 1;
}

bool BnB::SolveRec(size_t pos, bool b)
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

	sys_x[pos] = 1;
	bool f = SolveLP();		
	sys_x[pos] = 0;

	Cleanup(Triplets.size() - 1 - b, nr_rows - 1 - b);
	return f;
}
