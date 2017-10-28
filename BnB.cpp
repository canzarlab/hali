#include "BnB.h"
#include "Timer.h"
#include <iostream>

// ./solver inputs/C5.new.dag inputs/C5.new.map inputs/C32.new.dag inputs/C32.new.map align 0 j 1 0.025 0 2 2>/dev/null
// ./solver inputs/C5.new.dag inputs/C5.new.map inputs/C32.new.dag inputs/C32.new.map align 1 j 1 0.025 0 2 2>/dev/null
// ./solver inputs/C5.new.dag inputs/C5.new.map inputs/C32.new.dag inputs/C32.new.map align 2 j 1 0.025 0 2 2>/dev/null 185.395

// 182.445
// 182.469
// 185.334

BnB::BnB(Graph& t1, Graph& t2, string d, double k, bool dag, double c) : LP(t1, t2, d, k, dag), G(Greedy(t1, t2, d, k, dag)), con_eps(c)
{
}

int geno_calls = 0;
double geno_time = 0;

void BnB::Solve(string filename)
{
	MatchingConstraints();	

	G.Solve("");
	sys_lb = -G.GetSolution();

	x = Vector::Zero(nr_cols);
	sys_lo.conservativeResizeLike(Vector::Zero(nr_cols));
	sys_hi.conservativeResizeLike(Vector::Ones(nr_cols));
	sys_x.resize(nr_cols);	
	sys_x.assign(nr_cols, false);

	if (SolveLP())	
	{
		x = sys_s;
		for (size_t i = 0; i < x.size(); ++i)
			x(i) = round(x(i));
		WriteSolution(filename);
	}	
	else
		G.WriteSolution(filename);

	//cout << endl;
	//cout << "total geno calls: " << geno_calls << endl;
	//cout << "total geno time: " << geno_time << endl;
}

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
		y = Vector::Zero(nr_rows);
		Vector d = -c;			

	    BranchingJRF simpleJRF(A, b, d, x, y);
		simpleJRF.lo = sys_lo;
		simpleJRF.hi = sys_hi;
	    AugmentedLagrangian solver(simpleJRF, 15);
	    solver.setParameter("verbose", false);
	    solver.setParameter("pgtol", 1e-1); 
	    solver.setParameter("constraintsTol", 1e-4);

		//T.start();
	    SolverStatus status = solver.solve();	 
		//T.stop();

		/*qwe++;
		qwer += T.secs();
		geno_calls++;
		geno_time += T.secs();*/

		if (status == INFEASIBLE) 
		{				
			Cleanup(nr_t, nr_r);			
			return 0; 
		}

		x = Vector::ConstMapType(solver.x(), nr_cols); 

		if (LP::cf == 1 && Add<1>())
			continue;
        else if (LP::cf == 2 && (Add<1>() + Add<2>()))
			continue;

		f = solver.f();	

		/*cout << "rows: " << nr_rows << endl;
		cout << "cols: " << nr_cols << endl;
		cout << "geno calls: " << qwe << endl;
		cout << "geno time: " << qwer << endl;
		cout << "total geno calls: " << geno_calls << endl;
		cout << "total geno time: " << geno_time << endl;
		cout << "upper: " << f << endl;
		cout << "lower: " << sys_lb << endl;
		int qwert = 0;
		for (size_t i = 0; i < x.size(); ++i)
			if (x(i) > 1e-3 && x(i) < 1 - 1e-3 && !sys_x[i])
			 	qwert++;
		cout << "fracs: " << qwert << endl;*/

		if (sys_lb != double(INF) && f >= sys_lb * (1.0 + con_eps)) 
		{				
			//cout << "The branch was cut." << endl << endl;
			Cleanup(nr_t, nr_r);			
			return 0; 
		}			

		//cout << endl;

		break;
	}
	
	size_t pos = x.size();
	double val = 1;

	for (size_t i = 0; i < x.size(); ++i)
		if (x(i) > 1e-3 && x(i) < 1 - 1e-3 && !sys_x[i])
		 	if (abs(0.5 - x(i)) < val) pos = i, val = abs(0.5 - x(i));				

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
