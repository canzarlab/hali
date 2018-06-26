/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. 
*/

#ifndef Parallel_H
#define Parallel_H

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "BnG.h"

#define MAX_THREADS 9

class ParallelSolver
{

	public:
  
	ParallelSolver(Graph& t1, Graph& t2, string d, double k, bool dag, int nthreads);
  
	// Solves the model and writes it down to "filename".
  void Solve(string filename);

	// Updates the best upper bound and solution with regards to locking.
	void PushUB(Vector& var, double  val);
	void PullUB(Vector& var, double& val);

	double  GetUBVal(); 
	Vector& GetUBVar(); 

	bool    Finished(); 

	private:

	void Callback(string filename, GenericBnBSolver* solver);
	
	// threading locks and return values
	mutex							 thr_block; // Locks the best upper bound.
	mutex							 thr_slock; // Locks the solution.
	condition_variable thr_cond;  // Finishes all threads once one solver finds an optimal solution.
	atomic<thread::id> thr_val;   // Value of the finished thread.
	int                thr_num;   // Number of running threads.

	// Best upper bound and solution.
	Vector sys_sol;
	double sys_ub;

	// Solver related data.
	Graph& t1, t2;
	string d;
	double k;
	bool dag;

	bool sol;
};

// Different BnB solvers

#define PAR_CLASS(X, P, VS)                                                                 \
class X : public P																																					\
{																																														\
	public:																																										\
                                                                                            \
	X(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) :           \
		P(t1, t2, dist, k, dag), par(par)                                                       \
	{                                                                                         \
	}                                                                                         \
                                                                                            \
	protected:                                                                                \
                                                                                            \
	double VarScore       (int i) { return VS; }                                              \
	void   OnNodeStart    ()      { finished = par.Finished(); par.PullUB(sys_sol, sys_ub); } \
	void   OnUpdateUB     ()      { par.PushUB(sys_sol, sys_ub); }                            \
	void   OnSolverFinish ()      { par.PullUB(sys_sol, sys_ub); }                            \
                                                                                            \
	ParallelSolver& par;                                                                      \
};

PAR_CLASS(BnBBFMF, BFBnBSolver, 0.5 - abs(0.5 - x(i)))
PAR_CLASS(BnBBFLF, BFBnBSolver, abs(0.5 - x(i)))
PAR_CLASS(BnBBFWF, BFBnBSolver, c(i) * x(i))

PAR_CLASS(BnBDFMF, DFBnBSolver, 0.5 - abs(0.5 - x(i)))
PAR_CLASS(BnBDFLF, DFBnBSolver, abs(0.5 - x(i)))
PAR_CLASS(BnBDFWF, DFBnBSolver, c(i) * x(i))

PAR_CLASS(BnBHMF, HybridBnBSolver, 0.5 - abs(0.5 - x(i)))
PAR_CLASS(BnBHLF, HybridBnBSolver, abs(0.5 - x(i)))
PAR_CLASS(BnBHWF, HybridBnBSolver, c(i) * x(i))

#endif
