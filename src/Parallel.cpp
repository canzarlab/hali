/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

		./hali data0/a1 data0/a1 tmp_outs/test 2 s 1 0 0 20
		
		7BnBBFLF -364
    TOTAL TIME : 		546.771 secs
    453 

    7BnBBFMF -364.002
    TOTAL TIME : 		1080.77 secs 
    453
    
    64 5064
*/

#include "Parallel.h"

GenericBnBSolver* MakeSolver(Graph& t1, Graph& t2, string d, double k, bool dag, int threadid, ParallelSolver& s)
{
	if (threadid < 0 || threadid >= MAX_THREADS)
		return nullptr;
	else if (threadid == 0)
		return new BnBBFMF(t1, t2, d, k, dag, s);
	else if (threadid == 1)
		return new BnBBFLF(t1, t2, d, k, dag, s);
	else if (threadid == 2)
		return new BnBBFWF(t1, t2, d, k, dag, s);
	else if (threadid == 3)
		return new BnBDFMF(t1, t2, d, k, dag, s);
	else if (threadid == 4)
		return new BnBDFLF(t1, t2, d, k, dag, s);
	else if (threadid == 5)
		return new BnBDFWF(t1, t2, d, k, dag, s);
	else if (threadid == 6)
		return new BnBHMF (t1, t2, d, k, dag, s);
	else if (threadid == 7)
		return new BnBHLF (t1, t2, d, k, dag, s);
	else if (threadid == 8)
		return new BnBHWF (t1, t2, d, k, dag, s);
	else if (threadid == 9)
		return new BnBHA  (t1, t2, d, k, dag, s);
	else if (threadid == 10)
		return new BnBBFA (t1, t2, d, k, dag, s);
	else if (threadid == 11)
		return new BnBDFA (t1, t2, d, k, dag, s);
}

ParallelSolver::ParallelSolver(Graph& t1, Graph& t2, string d, double k, bool dag, int nthreads) :
	thr_num(max(1, min(nthreads, MAX_THREADS))), t1(t1), t2(t2), d(d), k(k), dag(dag), sol(0), sys_ub(666)
{
}

double ParallelSolver::GetUBVal()
{
	return sys_ub;
}
 
Vector& ParallelSolver::GetUBVar()
{
	return sys_sol;
}

bool ParallelSolver::Finished()
{
	thr_slock.lock();  
	bool b = sol; 
	thr_slock.unlock();
	return sol; 
}

void ParallelSolver::Callback(string filename, GenericBnBSolver* solver)
{
  #if DEBUG == 1
  solver->debug_file = filename + '_' + typeid(*solver).name();
  #endif
	solver->Solve("");	
	if (!sol && filename != "") 
	{
		thr_slock.lock(); 
		sol = 1; thr_val = this_thread::get_id();
		cout << typeid(*solver).name() << endl;
 		thr_slock.unlock();
		solver->WriteSolution(filename);	
	}
	thr_cond.notify_all();
}

void ParallelSolver::Solve(string filename)
{
	GenericBnBSolver* S[thr_num];
	thread T[thr_num];
	mutex m; 
	
	for (int i = 0; i < thr_num; ++i)
	{
		S[i] = MakeSolver(t1, t2, d, k, dag, i, *this);
		T[i] = thread{&ParallelSolver::Callback, this, filename, S[i]};
	}

	unique_lock<mutex> lock{m};
	thr_cond.wait(lock, [&] { return thr_val != thread::id{}; });

	for (int i = 0; i < thr_num; ++i)
	{
		T[i].join();
		delete S[i];
	}
}

/*void ParallelSolver::Solve(string filename)
{
	GenericBnBSolver* S; thread T; mutex m; 	
	S = MakeSolver(t1, t2, d, k, dag, 1, *this);
	T = thread{&ParallelSolver::Callback, this, filename, S};
	unique_lock<mutex> lock{m};
	thr_cond.wait(lock, [&] { return thr_val != thread::id{}; });
	T.join(); delete S;
}*/

void ParallelSolver::PushUB(Vector& var, double val, GenericBnBSolver& solver)
{
	thr_block.lock();
	if (val < sys_ub)
	{
	  #if DEBUG == 1
		solver.debug_log << "push: " << val << " > " << sys_ub << endl;
		#endif
		sys_ub  = val;
		sys_sol = var;
	}
	thr_block.unlock();
}

void ParallelSolver::PullUB(Vector& var, double& val, GenericBnBSolver& solver)
{
	thr_block.lock();
	if (val > sys_ub)
	{
	  #if DEBUG == 1
		solver.debug_log << "pull: " << val << " < " << sys_ub << endl;
		#endif
		val = sys_ub;
		var = sys_sol;
	}
	thr_block.unlock();
}

#define AGGRESSIVE_VAR(X)                                                                \
double X::VarScore(int i, BnBNode* node)                                                 \
{                                                                                        \
  int m, n, c = 0;                                                                       \
  for (int j = 0; j < t1.GetNumNodes(); ++j)                                             \
    for (int k = 0; k < t2.GetNumNodes(); ++k)                                           \
      if (K[j][k] == i) m = j, n = k;                                                    \
  for (int j = 0; j < t1.GetNumNodes(); ++j)                                             \
    for (int k = 0; k < t2.GetNumNodes(); ++k)                                           \
      if (K[j][k] != -1 && !(node->IsVarFixed(K[j][k])) && !IsNotInConflict(m, j, n, k)) \
        ++c;                                                                             \
  return c;                                                                              \
}

AGGRESSIVE_VAR(BnBHA)
AGGRESSIVE_VAR(BnBBFA)
AGGRESSIVE_VAR(BnBDFA)
