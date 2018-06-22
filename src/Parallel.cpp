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

atomic<thread::id> val;
condition_variable cond;

ParallelSolver::ParallelSolver(Graph& t1, Graph& t2, string d, double k, bool dag, int nthreads) :
	N(nthreads)
{
	cout << "solver created" << endl;
}

void ParallelSolver::Callback(string filename) //, GenericBnBSolver* solver)
{
	val = this_thread::get_id();	
	for (int i = 0; i < (int)1e10; ++i);
	cond.notify_all();
}

void ParallelSolver::Solve(string filename)
{
	cout << "solving" << endl;

	mutex m;
	
	for (int i = 0; i < 10; ++i)
		thread{&ParallelSolver::Callback, this, "test"}.detach();

	unique_lock<mutex> lock{m};
	cond.wait(lock, [&] { return val != thread::id{}; });

	cout << "Thread " << val << " finished first" << endl;
}



