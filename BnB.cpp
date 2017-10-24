#include "BnB.h"
#include "Timer.h"
#include <iostream>

BnB::BnB(Graph& t1, Graph& t2, string d, double k, bool dag) : LP(t1, t2, d, k, dag), G(Greedy(t1, t2, d, k, dag))
{
}

void BnB::Solve(string filename)
{
	MatchingConstraints();	

	G.Solve("");
	sys_lb = -G.GetSolution();
	float g = sys_lb;

	geno_calls = 0;
	geno_time  = 0.0;

	x = Vector::Zero(nr_cols);
	sys_lo.conservativeResizeLike(Vector::Zero(nr_cols));
	sys_hi.conservativeResizeLike(Vector::Ones(nr_cols));
	sys_x.resize(nr_cols);	
	sys_x.assign(nr_cols, false);

	SolveLP();
	if (sys_s.size() == x.size())	
	{
		x = sys_s;
		for (size_t i = 0; i < x.size(); ++i)
			x(i) = round(x(i));
		WriteSolution(filename);
	}	
	else
		G.WriteSolution(filename);

	//cout << "total geno calls: " << geno_calls << endl;
	//cout << "total geno time: " << geno_time << endl;

	/*float weight = 0.0;
	for (size_t i = 0; i < x.size(); ++i)
	{
		x(i) = round(x(i));		
		weight += x(i) * c(i);	
	}*/

	//cout << endl;
	//cout << "Greedy: " << -g << endl;
	//cout << "Optimal: " << weight << endl;
    //WriteSolution(filename);
}

// ./solver inputs/T2_s200_25.ploidyless.dag inputs/T2_s200_25.ploidyless.map inputs/T9_s200_25.ploidyless.dag inputs/T9_s200_25.ploidyless.map align 2 j 1 2

void BnB::Cleanup(size_t nr_t, size_t nr_r)
{
    Triplets.resize(nr_t);	
	nr_rows = nr_r;
}

bool BnB::SolveLP()
{
    int nr_t = Triplets.size();
	int nr_r = nr_rows;
	double f;

	int qwe = 0;
	double qwer = 0.0;
	Timer T;

	while(1)
	{
	    SpMat A(nr_rows, nr_cols);
	    A.setFromTriplets(Triplets.begin(), Triplets.end());
		Vector b = Vector::Ones(nr_rows);   	
		//x = Vector::Zero(nr_cols);
		y = Vector::Zero(nr_rows);
		Vector d = -c;			

	    BranchingJRF simpleJRF(A, b, d, x, y);
		simpleJRF.lo = sys_lo;
		simpleJRF.hi = sys_hi;
	    AugmentedLagrangian solver(simpleJRF, 15);
	    solver.setParameter("verbose", false);
	    solver.setParameter("pgtol", 1e-1); 
	    solver.setParameter("constraintsTol", 1e-4);

		T.start();
	    SolverStatus status = solver.solve();	 
		T.stop();

		qwe++;
		qwer += T.secs();
		geno_calls++;
		geno_time += T.secs();

		if (status == INFEASIBLE) 
		{	
			clog << "ERROR: infeasible solution" << endl;			
			Cleanup(nr_t, nr_r);			
			return 0; 
		}

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

		//cout << "geno calls: " << qwe << endl;
		//cout << "geno time: " << qwer << endl;
		//cout << "total geno calls: " << geno_calls << endl;
		//cout << "total geno time: " << geno_time << endl;
		//cout << "upper: " << f << endl;
		//cout << "lower: " << sys_lb << endl;
		
		if (sys_lb != double(INF) && f >= sys_lb * 1.0) 
		{				
			//cout << "The branch was cut." << endl << endl;
			Cleanup(nr_t, nr_r);			
			return 0; 
		}			

		//cout << endl;

		// ./solver T9_s200_25.ploidyless.dag T9_s200_25.ploidyless.map T2_s200_25.ploidyless.dag T2_s200_25.ploidyless.map align 2 j 1 2
		// ./solver inputs/a28 inputs/a10028 outputs/out 2 j 1 2 2>/dev/null > stats2.txt

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
	sys_lo(pos) = b;
	sys_hi(pos) = b;

	sys_x[pos] = 1;
	bool f = SolveLP();		
	sys_x[pos] = 0;

	sys_lo(pos) = 0;
	sys_hi(pos) = 1;

	return f;
}
