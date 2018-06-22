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
#include <chrono>
#include <mutex>
#include <condition_variable>

#include "BnG.h"

class ParallelSolver
{

	public:
  
	ParallelSolver(Graph& t1, Graph& t2, string d, double k, bool dag, int nthreads);
  
  void Solve(string filename);

	void UpdateUB(Vector& var, double val);
	double  GetUBVal();
	Vector& GetUBVar();

	private:

	void Callback(string filename);
	
	Vector sys_sol;
	double sys_ub;

	int N;
};

class BnBBFMF : public BFBnBSolver // Best first, most fractional
{
	public:

	BnBBFMF(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) : 
		BFBnBSolver(t1, t2, dist, k, dag), par(par)
	{
	}

	void OnUpdateUB(Vector& var, double val)
	{
		par.UpdateUB(var, val);
	}

	protected:
	
	virtual double VarScore(int i) { return 0.5 - abs(0.5 - x(i)); }

	ParallelSolver& par;
};

class BnBBFLF : public BFBnBSolver // Best first, least fractional
{
	public:

	BnBBFLF(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) : 
		BFBnBSolver(t1, t2, dist, k, dag), par(par)
	{
	}

	protected:
	
	virtual double VarScore(int i) { return abs(0.5 - x(i)); }

	ParallelSolver& par;
};

class BnBBFWF : public BFBnBSolver // Best first, weight times fractional
{
	public:

	BnBBFWF(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) : 
		BFBnBSolver(t1, t2, dist, k, dag), par(par)
	{
	}

	protected:
	
	virtual double VarScore(int i) { return x(i) * c(i); }

	ParallelSolver& par;
};

class BnBDFMF : public DFBnBSolver // Depth first, most fractional
{
	public:

	BnBDFMF(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) : 
		DFBnBSolver(t1, t2, dist, k, dag), par(par)
	{
	}

	protected:
	
	virtual double VarScore(int i) { return 0.5 - abs(0.5 - x(i)); }

	ParallelSolver& par;
};

class BnBDFLF : public DFBnBSolver // Depth first, least fractional
{
	public:

	BnBDFLF(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) : 
		DFBnBSolver(t1, t2, dist, k, dag), par(par)
	{
	}

	protected:
	
	virtual double VarScore(int i) { return abs(0.5 - x(i)); }

	ParallelSolver& par;
};

class BnBDFWF : public DFBnBSolver // Depth first, weight times fractional
{
	public:

	BnBDFWF(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) : 
		DFBnBSolver(t1, t2, dist, k, dag), par(par)
	{
	}

	protected:
	
	virtual double VarScore(int i) { return x(i) * c(i); }

	ParallelSolver& par;
};

class BnBHMF : public HybridBnBSolver // Hybrid, most fractional
{
	public:

	BnBHMF(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) : 
		HybridBnBSolver(t1, t2, dist, k, dag), par(par)
	{
	}

	protected:
	
	virtual double VarScore(int i) { return 0.5 - abs(0.5 - x(i)); }

	ParallelSolver& par;
};

class BnBHLF : public HybridBnBSolver // Hybrid, least fractional
{
	public:

	BnBHLF(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) : 
		HybridBnBSolver(t1, t2, dist, k, dag), par(par)
	{
	}

	protected:
	
	virtual double VarScore(int i) { return abs(0.5 - x(i)); }

	ParallelSolver& par;
};

class BnBHWF : public HybridBnBSolver // Hybrid, weight times fractional
{
	public:

	BnBHWF(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) : 
		HybridBnBSolver(t1, t2, dist, k, dag), par(par)
	{
	}

	protected:
	
	virtual double VarScore(int i) { return x(i) * c(i); }

	ParallelSolver& par;
};


#endif
