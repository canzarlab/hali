/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

		./hali inputs/a0 inputs/a999 align 2 s 1 0 0.001 2
*/

#include "BnG.h"
#include "Timer.h"

#define PGTOL 0.1 

BnBNode::BnBNode(vector<ET>& Triplets, size_t rows, size_t cols) : Triplets(Triplets), rows(rows), cols(cols), warm(nullptr), obj(1)
{
	var_lb.conservativeResizeLike(Vector::Zero(cols));
	var_ub.conservativeResizeLike(Vector::Ones(cols));
	#if DEBUG == 1	
	debug_depth = 0;
	#endif
}

BnBNode::BnBNode(BnBNode* node) : Triplets(node->Triplets), rows(node->rows), cols(node->cols), obj(1)
{
	warm = new Vector(node->sol);
	var_lb = node->var_lb;
	var_ub = node->var_ub;
	#if DEBUG == 1	
	debug_depth = 1 + node->debug_depth;
	#endif
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

GenericBnBSolver::GenericBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : LP(t1, t2, dist, k, dag), sys_ub(0), min_c(0), finished(0)
{	
}

void GenericBnBSolver::Solve(string filename)
{	
	#if DEBUG == 1	
	cout << "BnB debug mode on." << endl;
	debug_log = ofstream(filename + ".log");
	debug_nodecnt  = 0;
	debug_genocnt  = 0;
	debug_genotime = 0;
	debug_genomin  = 1e20;
	debug_genomax  = 0;
	#endif

	#if DEBUG == 1	
	Timer T; T.start();
	#endif

	PushNode(InitNodeFrom(nullptr));

	while (!OpenEmpty() && !finished)
	{
		BnBNode* open = EvalOpen();
		if (open != nullptr) 
		{
			vector<BnBNode*>* V = EvalBranch(open);
	
			if (V != nullptr)
				PushAll(V);	
			else if (sys_ub > open->obj)
			{
				sys_sol = open->sol; 
				sys_ub  = open->obj; 
				#if DEBUG == 1	
				debug_log << "BOUND " << sys_ub << " (" << open->debug_nodeid << ")" << endl << endl;
				#endif			
				OnUpdateUB();
			}	

			delete open;
		}
	}
	#if DEBUG == 1	
	T.stop();
	#endif

	if (!finished)
	{
		OnSolverFinish();

		x = sys_sol;
		for (size_t i = 0; i < x.size(); ++i)
			x(i) = round(x(i));	

		if (filename != "") 
			WriteSolution(filename);	 
	}

	#if DEBUG == 1	
	debug_log << "total time: " << T.secs() << endl;
	debug_log << "total nodes: " << debug_nodecnt << endl;
	debug_log << "total geno calls: " << debug_genocnt << endl;
	debug_log << "total geno time: " << debug_genotime << " (" << debug_genomin << ", " << debug_genotime / debug_genocnt << "," << debug_genomax << ")" << endl;
	debug_log << "best: " << -sys_ub << endl;
	#endif
}

BnBNode* GenericBnBSolver::InitNodeFrom(BnBNode* node)
{
	BnBNode* newnode;
	if (node == nullptr)
		newnode = new BnBNode(Triplets, nr_rows, nr_cols);
	else
		newnode = new BnBNode(node);
	#if DEBUG == 1
	newnode->debug_nodeid = debug_nodecnt++;
	newnode->debug_parent = (node == nullptr) ? 0 : node->debug_nodeid;	
	#endif
	return newnode;
}

bool GenericBnBSolver::SolveNode(BnBNode* node, double pgtol, double numtol)
{
	Triplets = node->Triplets;
	nr_rows  = node->rows;
	nr_cols  = node->cols;
	x        = (node->warm) ? *(node->warm) : Vector::Zero(nr_cols);

	#if DEBUG == 1	
	size_t debug_node_genocnt  = 0;
	double debug_node_genotime = 0;
	debug_log << "id: " << node->debug_nodeid << endl;
	debug_log << "parent: " << node->debug_parent << endl;	
	debug_log << "depth: " << node->debug_depth << endl;
	debug_log << "best: " << -sys_ub << endl;
	#endif

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

		#if DEBUG == 1	
		Timer debug_T; debug_T.start();		
		#endif
		OnNodeStart();
		if (finished) return false;
	  SolverStatus status = solver.solve();	
		#if DEBUG == 1	
		debug_T.stop();
		debug_genocnt++;	
		debug_node_genocnt++;
		debug_node_genotime += debug_T.secs();
		debug_genotime += debug_T.secs();
		debug_genomin = min(debug_genomin, (double)debug_T.secs());
		debug_genomax = max(debug_genomax, (double)debug_T.secs());
		if (status == SUBOPTIMAL)
	    debug_log << "SUBOPTIMAL(" << status << ")" << endl; 
		if (status == UNBOUNDED)
	    debug_log << "UNBOUNDED(" << status << ")" << endl; 
	  if (status == INFEASIBLE)
	    debug_log << "INFEASIBLE(" << status << ")" << endl; 
	  if (status == NUM_ERROR)
	    debug_log << "NUM_ERROR(" << status << ")" << endl; 
		#endif
		
		if (status == INFEASIBLE || CheckUB(solver.f(), numtol))
		{
			#if DEBUG == 1	
			debug_log << "geno calls: " << debug_node_genocnt << endl;
			debug_log << "geno time: " << debug_node_genotime << endl;
			debug_log << "CUT" << endl << endl;
			#endif
			cout << "Cut: " << solver.f() << ' ' << sys_ub << endl;
			return false;
		} 

		x = Vector::ConstMapType(solver.x(), nr_cols); 
		
		if ((LP::cf == 1 && Add<1>()) || (LP::cf == 2 && (Add<1>() + Add<2>())))
			continue;

		node->Triplets = Triplets;
		node->rows     = nr_rows;
		node->cols     = nr_cols;
		node->obj      = solver.f();		
		node->sol      = x;

    break;
	}

	#if DEBUG == 1	
	debug_log << "geno calls: " << debug_node_genocnt << endl;
	debug_log << "geno time: " << debug_node_genotime << endl;
	debug_log << node->rows << " x " << node->cols << endl;
	debug_log << "lp: " << -node->obj << endl;
	size_t debug_fracs = 0;
	for (size_t i = 0; i < x.size(); ++i)
		debug_fracs += IsVarFrac(node->sol(i)) && !node->IsVarFixed(i);
	debug_log << "fracs: " << debug_fracs << endl << endl;
	#endif

	return true;
}

vector<BnBNode*>* GenericBnBSolver::EvalBranch(BnBNode* node)
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

BnBNode* GenericBnBSolver::MakeNode(BnBNode* parent, size_t index, double val)
{
	BnBNode* node = InitNodeFrom(parent);
	node->FixVar(index, val);
	return node;
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
	BnBNode* node = Open.top(); Open.pop();
	return SolveNode(node, PGTOL, 0.001) ? node : nullptr;
}

// ================= BF BNB SOLVER ==================== // 

BFBnBSolver::BFBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : GenericBnBSolver(t1, t2, dist, k, dag)
{	
}

void BFBnBSolver::PushNode(BnBNode* node)
{
	Open.push_back(node);
}

bool BFBnBSolver::OpenEmpty()
{
	return !Open.size();
}

BnBNode* BFBnBSolver::EvalOpen()
{
	for (auto& it : Open)
		if (it->obj == 1 && !SolveNode(it, PGTOL, 0.001))
			it->obj = 2;

	sort(Open.begin(), Open.end(), [](BnBNode*& l, BnBNode*& r)
	{ 
		return l->obj > r->obj; 
	});

	BnBNode* node = Open.back(); Open.pop_back();		
	return (node->obj < 1) ? node : nullptr;
}

// ================= DF BNB SOLVER ==================== //

DFBnBSolver::DFBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : GenericBnBSolver(t1, t2, dist, k, dag)
{	
}

void DFBnBSolver::PushNode(BnBNode* node)
{
	Open.push_back(node);
}

bool DFBnBSolver::OpenEmpty()
{
	return !Open.size();
}

BnBNode* DFBnBSolver::EvalOpen()
{
	BnBNode* lt = Open.back(); Open.pop_back();
	if (Open.size() > 1)
	{
		BnBNode* rt = Open.back(); 
		if (lt->obj == 1 && !SolveNode(lt, PGTOL, 0.001)) lt->obj = 2;
		if (lt->warm == rt->warm) // TODO this is horrible.
		{
			Open.pop_back();
			if (rt->obj == 1 && !SolveNode(rt, PGTOL, 0.001)) rt->obj = 2;
			if (lt->obj > rt->obj) swap(lt, rt);
			if (rt->obj < 1) Open.push_back(rt); 		
		}		
	}
	return (lt->obj < 1 || SolveNode(lt, PGTOL, 0.001)) ? lt : nullptr;
}

// ================= Hybrid BNB SOLVER ==================== //

HybridBnBSolver::HybridBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : GenericBnBSolver(t1, t2, dist, k, dag)
{	
}

void HybridBnBSolver::PushNode(BnBNode* node)
{
	Open.push_back(node);
}

bool HybridBnBSolver::OpenEmpty()
{
	return !Open.size();
}

BnBNode* HybridBnBSolver::EvalOpen() // TODO this is horrible.
{
	if (!GetObjective())
	{
		BnBNode* lt = Open.back(); Open.pop_back();
		if (Open.size() > 1)
		{
			BnBNode* rt = Open.back(); 
			if (lt->obj == 1 && !SolveNode(lt, PGTOL, 0.001)) lt->obj = 2;
			if (lt->warm == rt->warm)
			{
				Open.pop_back();
				if (rt->obj == 1 && !SolveNode(rt, PGTOL, 0.001)) rt->obj = 2;
				if (lt->obj > rt->obj) swap(lt, rt);
				if (rt->obj < 1) Open.push_back(rt); 		
			}		
		}
		return (lt->obj < 1 || SolveNode(lt, PGTOL, 0.001)) ? lt : nullptr;
	}
	else
	{
		for (auto& it : Open)
			if (it->obj == 1 && !SolveNode(it, PGTOL, 0.001))
				it->obj = 2;

		sort(Open.begin(), Open.end(), [](BnBNode*& l, BnBNode*& r)
		{ 
			return l->obj > r->obj; 
		});

		BnBNode* node = Open.back(); Open.pop_back();		
		return (node->obj < 1) ? node : nullptr;
	}
}

