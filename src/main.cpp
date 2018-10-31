/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "Timer.h"
#include "Greedy.h"
#include "LP.h"
#include "BnB.h"
#include "BnG.h"
#include "LPInt.h"
#include "LPCP.h"
#include "LPFInt.h"
#include "Similarity.h"
#include "Parallel.h"

#include <iostream>
#include <numeric>
int Solver::cf;
bool Solver::tt;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<vvvi> vvvvi;
typedef vector<string> vs;
typedef vector<vs> vvs;

int GetNumTrees(newick_node* root, vn& N)
{
    int k = 0, s = 1, f = 1;
    N.push_back(root);
    for (newick_child* child = root->child; child; child = child->next)
    {
        f *= ++k;
        s *= GetNumTrees(child->node, N);
    }
    return s * f;
}

void GenPerms(vn& v, int k, vvi S, vvvi& P)
{
    if (k == v.size())
    {
        P.push_back(S);
        return;
    }
    int c = 0;
    for (newick_child* child = v[k]->child; child; child = child->next)
        c++;

    vector<int> I(c);
    iota(I.begin(), I.end(), 0);
    do
    {
        S.push_back(I);
        GenPerms(v, k + 1, S, P);
        S.pop_back();
    } while (next_permutation(I.begin(), I.end()));
}

int MakeTree(newick_node* root, newick_node* nroot, int k, vvi& I)
{
    vector<newick_node*> c;
    for (newick_child* child = root->child; child; child = child->next)
        c.push_back(child->node);

    vector<newick_node*> nc;
    newick_child** childptr = &nroot->child;
    vi J(I[k].size());
    for (int i = 0; i < I[k].size(); ++i)
    {
        newick_node* croot = new newick_node(c[I[k][i]]->taxon);
        *childptr = new newick_child(croot);
        childptr = &(*childptr)->next;
        nc.push_back(croot);
        J[I[k][i]] = i;
    }
    int s = 1;
    for (int i = 0; i < c.size(); ++i)
        s += MakeTree(c[i], nc[J[i]], k + s, I);
    return s;
}

int CalcSubtreeSizes(newick_node* root, map<newick_node*, int>& S)
{
    int s = 1;
    for (newick_child* child = root->child; child; child = child->next)
        s += CalcSubtreeSizes(child->node, S);
    return S[root] = s;
}

int hack5(vvvvi& DP, vvi& T1, vvi& T2, vvs& L1, vvs& L2, int n, int m, int i, int j, int k, int l)
{
    if (i + j == n)
        return m - k - l;
    if (k + l == m)
        return n - i - j;

    int p = n - T1[i][j + 1] - j, q = m - T2[k][l + 1] - l;
    int x = i ? p + 1 : p + 2, y = i ? j : j - 1;
    int xx = k ? q + 1 : q + 2, yy = k ? l : l - 1;
    if (i == 0 && j == 0)
        x = 1, y = 0;
    if (k == 0 && l == 0)
        xx = 1, yy = 0;

    assert(DP[i][j + 1][k][l] != -1);
    assert(DP[i][j][k][l + 1] != -1);
    assert(DP[i][j + T1[i][j + 1]][k][l + T2[k][l + 1]] != -1);
    assert(DP[x][y][xx][yy] != -1);
    assert(!L1[i][j].empty() && !L2[k][l].empty());

    clog << i << ' ' << j << ' ' << k << ' ' << l << endl << string(7, '=') << endl;
    clog << i << ' ' << j + 1 << ' ' << k << ' ' << l << endl;
    clog << i << ' ' << j << ' ' << k << ' ' << l + 1 << endl;
    clog << i << ' ' << j + T1[i][j + 1] << ' ' << k << ' ' << l + T2[k][l + 1] << endl;
    clog << x << ' ' << y << ' ' << xx << ' ' << yy << endl;
    clog << endl;

    int o = DP[i][j + 1][k][l] + 1;
    int a = DP[i][j][k][l + 1] + 1;
    int b = DP[i][j + T1[i][j + 1]][k][l + T2[k][l + 1]];
    int c = DP[x][y][xx][yy];
    int d = L1[i][j + 1] != L2[k][l + 1];
    int mm = min(min(o, a), b + c + d);
    return mm;
}

newick_node* GetLastRemoved(newick_node* root, int i, int j)
{
    list<newick_node*> D(1, root);
    while (i--)
    {
        root = D.front();
        D.pop_front();
        auto it = D.begin();
        for (newick_child* child = root->child; child; child = child->next)
            D.insert(it, child->node);
    }
    while (j--)
    {
        root = D.back();
        D.pop_back();
        for (newick_child* child = root->child; child; child = child->next)
            D.push_back(child->node);
    }
    return root;
}

void GenTables(newick_node* root, int n, vvs& L, vvi& T)
{
    map<newick_node*, int> S;
    CalcSubtreeSizes(root, S);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (newick_node* node = (i + j <= n ? GetLastRemoved(root, i, j) : nullptr))
                L[i][j] = node->taxon, T[i][j] = S[node];
}

int OrderedEditDist(Tree& t1, newick_node* t2, int m)
{
    int n = t1.GetNumNodes();
    vvs L1(n, vs(n)), L2(m, vs(m));
    vvi T1(n + 1, vi(n + 1, 1)), T2(m + 1, vi(m + 1, 1));
    GenTables(t1.GetRoot(), n, L1, T1);
    GenTables(t2, m, L2, T2);
    vvvvi DP(n + 1, vvvi(n + 1, vvi(m + 1, vi(m + 1, -1))));
    for (int i = n; i >= 0; --i)
        for (int j = n; j >= 0; --j)
            for (int k = m; k >= 0; --k)
                for (int l = m; l >= 0; --l)
                    if (i + j <= n && k + l <= m)
                        DP[i][j][k][l] = hack5(DP, T1, T2, L1, L2, n, m, i, j, k, l);

    return DP[0][0][0][0];
}

void DeallocTree(newick_node* root)
{
    for (newick_child* child = root->child; child; child = child->next)
        DeallocTree(child->node);
    delete root;
}

void EditDist(Graph& g1, Graph& g2)
{
    vn n1, n2;
    int k1 = GetNumTrees(g1.GetRoot(), n1), k2 = GetNumTrees(g2.GetRoot(), n2);
    Tree& t1 = (Tree&)(k1 >= k2 ? g1 : g2), &t2 = (Tree&)(k1 >= k2 ? g2 : g1);
    if (k2 > k1)
        swap(n1, n2);

    vvvi P;
    GenPerms(n2, 0, {}, P);
    vector<newick_node*> trees;
    for (int i = 0; i < P.size(); ++i)
    {
        newick_node* nroot = new newick_node(n2[0]->taxon);
        MakeTree(n2[0], nroot, 0, P[i]);
        trees.push_back(nroot);
    }

    int m = numeric_limits<int>::max();
    for (int i = 0; i < P.size(); ++i)
    {
        print_tree(t1.GetRoot(), cout);
        cout << ' ';
        print_tree(trees[i], cout);
        cout << ' ';
        int nm = OrderedEditDist(t1, trees[i], t2.GetNumNodes());
        m = min(m, nm);
        cout << nm << endl;
        DeallocTree(trees[i]);
    }
    cout << "MIN: " << m << endl;
}

Solver* MakeSolver(Graph& t1, Graph& t2, int argc, char** argv)
{
    int s = stoi(argv[argc - 1]);
    Solver::cf = stoi(argv[4 + (argc == 9) + 2 * (argc == 12)]);
    Solver::tt = argc == 12;
    string d = argv[5 + (argc == 9) + 2 * (argc == 12)];
    double k = stod(argv[6 + (argc == 9) + 2 * (argc == 12)]);
    double c = (s != 2) ? 0 : stod(argv[8 + 2 * (argc == 12)]);
    var_eps = (argc == 9) ? 0 : stod(argv[7 + 2 * (argc == 12)]);
    bool dag = argc == 9;

    assert(LP::cf >= 0 && LP::cf <= 2);
    assert(d == "j" || d == "s");
    assert(s >= 0 && s <= 9);

    if (s == 0)
        return new Greedy(t1, t2, d, k, dag);
    else if (s == 1)
        return new LP(t1, t2, d, k, dag);
    else if (s == 2)
        EditDist(t1, t2);
    else if (s == 3)
        return new LPCP(t1, t2, d, k, dag);
    else if (s == 4)
        return new DFBnBSolver(t1, t2, d, k, dag); // was test, to delete
    else if (s == 5)
        return new LPInt(t1, t2, d, k, dag);
    else if (s == 6)
        return new LPFInt(t1, t2, d, k, dag);
    else if (s == 7)
        return new BFBnBSolver(t1, t2, d, k, dag);
    else if (s == 8)
        return new DFBnBSolver(t1, t2, d, k, dag);
    else if (s == 9)
        return new HybridBnBSolver(t1, t2, d, k, dag);
    return nullptr;
}

Graph* MakeDAG(const char* f1, const char* f2, int s)
{
    return ((s == 0) ? new DAG(f1, f2) : new LDAG(f1, f2))->Init();
}

pair<Graph*, Graph*> MakeGraphs(int argc, char** argv)
{
    if (argc == 10)
        return {new Tree(argv[1]), new Tree(argv[2])};
    else if (argc == 12)
        return {new Tree(argv[1], argv[2]), new Tree(argv[3], argv[4])};

    int s = stoi(argv[argc - 1]);
    return {MakeDAG(argv[1], argv[2], s), MakeDAG(argv[3], nullptr, s)};
}

int main(int argc, char** argv)
{
    if (argc < 9 || argc > 12)
    {
        cout << "tree usage: " << argv[0] << " <filename.newick> <filename.newick> <align> <constraints> <weightfunc> <k> <vareps> <coneps> <solver>" << endl;
        cout << "dag usage: " << argv[0] << " <yeastnet> <mapping> <go> <align> <constraints> <weightfunc> <k> <solver>" << endl;
        cout << "tree usage (2): " << argv[0] << " <tree> <map> <tree> <map> <align> <constraints> <weightfunc> <k> <vareps> <coneps> <solver>" << endl;
        return EXIT_FAILURE;
    }

    Timer T;
    T.start();
    Graph *t1, *t2;
    tie(t1, t2) = MakeGraphs(argc, argv);
    if (stoi(argv[argc - 1]) > 9)
    {
        int s = stoi(argv[argc - 1]);
        Solver::cf = stoi(argv[4 + (argc == 9) + 2 * (argc == 12)]);
        Solver::tt = argc == 12;
        string d = argv[5 + (argc == 9) + 2 * (argc == 12)];
        double k = stod(argv[6 + (argc == 9) + 2 * (argc == 12)]);
        double c = (s != 2) ? 0 : stod(argv[8 + 2 * (argc == 12)]);
        var_eps = (argc == 9) ? 0 : stod(argv[7 + 2 * (argc == 12)]);

        assert(LP::cf >= 0 && LP::cf <= 2);
        assert(d == "j" || d == "s");

        ParallelSolver(*t1, *t2, d, k, argc == 9, s - 9).Solve(argv[3 + (argc == 9) + 2 * (argc == 12)]);
    }
    else
    {
        Solver* solver = MakeSolver(*t1, *t2, argc, argv);
        if (solver) solver->Solve(argv[3 + (argc == 9) + 2 * (argc == 12)]);
        delete solver;
        delete t1;
        delete t2;
    }
    T.stop();
    cout << "TIME: " << T.secs() << " secs" << endl;
}
