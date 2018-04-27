/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/

#include "BnG.h"
#include "Timer.h"
#include <iostream>

#define DEBUG 1

BnBNode::BnBNode(vector<ET>& Triplets, size_t rows, size_t cols) : Triplets(Triplets), rows(rows), cols(cols), warm(nullptr)
{
	var_lb.conservativeResizeLike(Vector::Zero(cols));
	var_ub.conservativeResizeLike(Vector::Ones(cols));
}

BnBNode::BnBNode(vector<ET>& Triplets, size_t rows, size_t cols, BnBNode* node) : Triplets(Triplets), rows(rows), cols(cols)
{
	warm = new Vector(node->sol);
	var_lb = node->var_lb;
	var_ub = node->var_ub;
}

BnBNode::~BnBNode()
{
	if (warm) delete warm;
}

bool BnBNode::IsVarFixed(size_t index)
{
	return var_lb(index) == var_ub(index);
}

void BnBNode::FixVar(size_t index, double val)
{
	var_lb(index) = var_ub(index) = val;
}

GenericBnBSolver::GenericBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : LP(t1, t2, dist, k, dag), sys_ub(0)
{	
}

void GenericBnBSolver::Solve(string filename)
{	
	#if DEBUG == 1	
	Timer T; T.start();
	#endif
	
	PushNode(InitNodeFrom(nullptr));

	while (!OpenEmpty())
	{
		BnBNode* open = EvalOpen();
		if (open != nullptr) 
		{
			vector<BnBNode*>* V = EvalBranch(open);
			
			if (V != nullptr)
				PushAll(V);
			else if (sys_ub > open->obj)
					sys_sol = open->sol, sys_ub = open->obj; 			
				
			delete open;
		}
	}
	#if DEBUG == 1	
	T.stop();
	#endif

	x = sys_sol;
	for (size_t i = 0; i < x.size(); ++i)
		x(i) = round(x(i));	

	if (filename != "") 
		WriteSolution(filename); 

	#if DEBUG == 1	
	cout << "OPT: " << -sys_ub << " in " << T.secs() << endl;
	#endif
}


BnBNode* GenericBnBSolver::InitNodeFrom(BnBNode* node)
{
	if (node == nullptr)
		return new BnBNode(Triplets, nr_rows, nr_cols);
	else
		return new BnBNode(Triplets, nr_rows, nr_cols, node);
}

bool GenericBnBSolver::SolveNode(BnBNode* node, double pgtol, double numtol)
{
	Triplets = node->Triplets;
	nr_rows  = node->rows;
	nr_cols  = node->cols;
	x        = (node->warm) ? *(node->warm) : Vector::Zero(nr_cols);

	while(1)
	{
	  SpMat A(nr_rows, nr_cols);
	  A.setFromTriplets(Triplets.begin(), Triplets.end());
	  Vector b = Vector::Ones(nr_rows);   	
		Vector y = Vector::Zero(nr_rows);		
		Vector d = -c;					

	  BranchingJRF simpleJRF(A, b, d, x, y);
		simpleJRF.lo = node->var_lb;
		simpleJRF.hi = node->var_ub;
	  AugmentedLagrangian solver(simpleJRF, 15);
	  solver.setParameter("verbose", false);
	  solver.setParameter("pgtol", pgtol); 
	  solver.setParameter("constraintsTol", 1e-4);

	  SolverStatus status = solver.solve();	
		#if DEBUG == 1		
		if (status == SUBOPTIMAL)
	    cout << "SUBOPTIMAL(" << status << ")" << endl; 
		if (status == UNBOUNDED)
	    cout << "UNBOUNDED(" << status << ")" << endl; 
	  if (status == INFEASIBLE)
	    cout << "INFEASIBLE(" << status << ")" << endl; 
	  if (status == NUM_ERROR)
	    cout << "NUM_ERROR(" << status << ")" << endl; 
		#endif

		if (status == INFEASIBLE) return false;
		if (solver.f() >= sys_ub * (1.0 + numtol)) return false;

		x = Vector::ConstMapType(solver.x(), nr_cols); 
		
		if ((LP::cf == 1 && Add<1>()) || (LP::cf == 2 && (Add<1>() + Add<2>())))
			continue;

		node->obj = solver.f();		
		node->sol = x;

    break;
	}

	#if DEBUG == 1	
	cout << nr_rows << " x " << nr_cols << endl;
	cout << "sol: " << -node->obj << endl;
	cout << endl;
	#endif

	return true;
}

void GenericBnBSolver::PushAll(vector<BnBNode*>* Nodes)
{ 
	for (auto& it : *Nodes) PushNode(it);
	delete Nodes;
}

// ================= TMP\TEST BNB SOLVER ==================== // 

TestBnBSolver::TestBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : GenericBnBSolver(t1, t2, dist, k, dag)
{	
}

void TestBnBSolver::PushNode(BnBNode* node)
{
	Open.push(node);
}

bool TestBnBSolver::OpenEmpty()
{
	return Open.empty();
}

BnBNode* TestBnBSolver::EvalOpen()
{
	BnBNode* node = Open.top();
	Open.pop();
	return SolveNode(node, 0.01, 0.001) ? node : nullptr;
}

BnBNode* TestBnBSolver::MakeNode(BnBNode* init, size_t p, bool b)
{
	BnBNode* node = InitNodeFrom(init);
	node->FixVar(p, b);
	return node;
}

vector<BnBNode*>* TestBnBSolver::EvalBranch(BnBNode* node)
{
	vector<pair<size_t, size_t>> vp;
	size_t p = 0;

	for (size_t i = 0; i < t1.GetNumNodes(); ++i)
	  for (size_t j = 0; j < t2.GetNumNodes(); ++j)
			if (K[i][j] != -1 && IsVarFrac(node->sol(K[i][j])) && !node->IsVarFixed(K[i][j]))				
			  vp.emplace_back(i, j);				

  sort(vp.begin(), vp.end(), [this](pair<size_t, size_t>& p, pair<size_t, size_t>& q)
  {
    return VarScore(K[p.first][p.second]) > VarScore(K[q.first][q.second]);  
  });

  for (int i = 0; !p && i < vp.size(); ++i)
  {
    bool flag = 1;
		for (size_t k = 0; flag && k < t1.GetNumNodes(); ++k)	  	
			for (size_t j = 0; flag && j < t2.GetNumNodes(); ++j)
				if (K[k][j] != -1 && node->IsVarFixed(K[k][j]))						  	
					flag = IsNotInConflict(vp[i].first, k, vp[i].second, j);					
    if (flag) p = 1 + K[vp[i].first][vp[i].second];
  }

	if (vp.size())
	{
		vector<BnBNode*>* V = new vector<BnBNode*>;
		bool b0 = !p; 
		p = (b0) ? K[vp[0].first][vp[0].second] : p - 1; 
		V->push_back(MakeNode(node, p, 0));		
		if (!b0) V->push_back(MakeNode(node, p, 1));
		return V;
	}

	return nullptr;
}
