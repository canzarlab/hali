/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/

#include "Parallel.h"

#include <sstream>

// atomic<thread::id> val;
// condition_variable cond;

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
	thr_num(max(1, min(nthreads, MAX_THREADS - 1))), t1(t1), t2(t2), d(d), k(k), dag(dag)
{
}

void ParallelSolver::Callback(string filename, GenericBnBSolver* solver)
{
	solver->Solve(filename);
	thr_val = this_thread::get_id();
	thr_cond.notify_all();
}

void ParallelSolver::Solve(string filename)
{
	GenericBnBSolver* S[thr_num];
	mutex m;
	
	for (int i = 0; i < thr_num; ++i)
	{
		S[i] = MakeSolver(t1, t2, d, k, dag, i, *this);
		thread{&ParallelSolver::Callback, this, filename, S[i]}.detach();
	}

	unique_lock<mutex> lock{m};
	thr_cond.wait(lock, [&] { return thr_val != thread::id{}; });

	for (int i = 0; i < thr_num; ++i)
		delete S[i];
}

void ParallelSolver::UpdateUB(Vector& var, double val)
{
	thr_lock.lock();
	if (val < sys_ub)
	{
		sys_ub  = val;
		sys_sol = var;
	}
	thr_lock.unlock();
}

// Solver implementation





