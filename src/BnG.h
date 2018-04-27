/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/

#ifndef BNG_H
#define BNG_H

#include "LP.h"
#include <vector>
#include <stack>

/*
	Class BnBNode

	It holds all the information required to evaluate LP at a single BnB node.
*/
class BnBNode
{
	public:

	// This constructor creates a node without a parent (eg. for creating BnB root node). 
	BnBNode(vector<ET>& Triplets, size_t rows, size_t cols); 

	// This constructor creates a node with a parent (last argument) and warm starts from it.
	BnBNode(vector<ET>& Triplets, size_t rows, size_t cols, BnBNode* node);

	~BnBNode();

	// Returns whether has BnB fixed a variable at index.
	bool       IsVarFixed(size_t index);

	// Fixes a variable at index to value val.
	void       FixVar(size_t index, double val);

	vector<ET> Triplets; // Stores the system matrix.
	size_t     rows;     // Stores the number of rows of system matrix.
	size_t     cols;		 // Stores the number of columns of system matrix.

	Vector     var_lb;   // Stores lower bounds for the variables.
	Vector     var_ub;	 // Stores upper bounds for the variables.
	
	Vector*    warm;     // Stores the warm start vector x.
	Vector     sol;      // Stores the solution vector x after the node has been evaluated.
	double     obj;      // Stores the objective value after the node has been evaluated.
};

/*
	Class GenericBnBSolver

	This is the template for the BnB solver which takes care of the basic functionallity.

	The idea.
	1. Create Open (set of BnBNode pointers)
	2. Open.Push(BnB_Root_Node)
	3. WHILE (NOT Open.Empty())
	4. 	Let (Node \in Open) be a solved node
	5.  IF there are fractionals in Node, add its children to Open

	Open is an arbitrary data structure for storing nodes waiting to be evaluated.
	In steps 4 and 5, different branching strateges may be applied. 	
*/
class GenericBnBSolver : public LP
{
	public:

	// Same constructor as for the LP solver.
	GenericBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag);

	// Solves the problem and writes the solution to the file. If no file is required pass an empty string.
	void Solve(string filename) override; 

	Vector GetSolution()  { return sys_sol; } // Returns the current solution (copy).
	double GetObjective() { return sys_ub; }  // Returns the current objective function value.

	protected:

	// Creates a BnBNode using current solver data and a designated parent node. 
	// If no parent node is needed, pass a nullptr.
	BnBNode*                  InitNodeFrom(BnBNode* node);

	// Solves the LP in a node using specified pgtol and numtol.
	//	pgtol  - geno solver parameter
	//	numtol - numerical offset used when cutting branches (cut if obj >= best_obj * (1 + numtol)) 
	virtual bool              SolveNode  (BnBNode* node, double pgtol, double numtol);	

	// Push a node into the Open set.
	virtual void              PushNode   (BnBNode* node) { }

	// Is the Open set empty? If it is, Solver will terminate.	
	virtual bool              OpenEmpty  ()              { return true; }	

	// Evaluate the Open set and decide which node to branch.	
	virtual BnBNode*          EvalOpen   ()              { return nullptr; }

	// Branch the node and return a vector of children nodes.
	// Return nullptr when there is nothing to branch into.
	virtual vector<BnBNode*>* EvalBranch (BnBNode* node) { }

	// Auxilliary function which decides whether a variable is to be considered fractional.
	virtual bool              IsVarFrac  (double val)    { return val > 0.001 && val < 0.999; }         

	private:

	// Pushes all nodes to the Open set.
	void PushAll(vector<BnBNode*>* Nodes);

	double sys_ub;  // Stores the best current upper bound.
	Vector sys_sol; // Stores the best current solution.
};

/*
	class TestBnBSolver
	
	Example on how to create BnB solvers using the generic solver. It is equivallent to hali option 2. 
*/
class TestBnBSolver : public GenericBnBSolver
{
	public:

	TestBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag);

	private:

	// Override when neccessary.
	void              PushNode   (BnBNode* node) override;
	bool              OpenEmpty  ()              override;
	BnBNode*          EvalOpen   ()              override;
	vector<BnBNode*>* EvalBranch (BnBNode* node) override;
	
	BnBNode* MakeNode(BnBNode* init, size_t p, bool b);

	// Fractionallity of variable i.
	double VarScore(int i) { return 0.5 - abs(0.5 - x(i)); } //{ return c(i); } // 

	stack<BnBNode*> Open; // Open set defined as a stack.
};

#endif
