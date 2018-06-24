/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

		./hali data0/a0 data0/a0 align 2 s 1 0 0 10
*/

#include "Parallel.h"

GenericBnBSolver* MakeSolver(Graph& t1, Graph& t2, string d, double k, bool dag, int threadid, ParallelSolver& s)
{
	if (threadid < 0 || threadid > MAX_THREADS - 1)
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
}

ParallelSolver::ParallelSolver(Graph& t1, Graph& t2, string d, double k, bool dag, int nthreads) :
	thr_num(max(1, min(nthreads, MAX_THREADS - 1))), t1(t1), t2(t2), d(d), k(k), dag(dag), sol(0), sys_ub(666)
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
	solver->Solve("");	
	
	if (!sol && filename != "") 
	{
		thr_slock.lock(); 
		sol = 1; thr_val = this_thread::get_id();
		
		cout << typeid(*solver).name() << ' ' << solver->GetObjective() << endl;

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
	S = MakeSolver(t1, t2, d, k, dag, 8, *this);
	T = thread{&ParallelSolver::Callback, this, filename, S};
	unique_lock<mutex> lock{m};
	thr_cond.wait(lock, [&] { return thr_val != thread::id{}; });
	T.join(); delete S;
}*/

void ParallelSolver::PushUB(Vector& var, double val)
{
	thr_block.lock();
	if (val < sys_ub)
	{
		cout << "push: " << val << " > " << sys_ub << endl;
		sys_ub  = val;
		sys_sol = var;
	}
	thr_block.unlock();
}

void ParallelSolver::PullUB(Vector& var, double& val)
{
	thr_block.lock();
	if (val > sys_ub)
	{
		cout << "pull: " << val << " < " << sys_ub << endl;
		val = sys_ub;
		var = sys_sol;
	}
	thr_block.unlock();
}




