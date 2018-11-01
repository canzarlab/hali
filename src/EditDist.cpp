/*
    Copyright (C) 2018 Mislav Blažević, Luka Borozan

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include <iostream>
#include <numeric>
#include <limits>
#include <cassert>
#include <algorithm>
#include "Graph.h"

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<vvvi> vvvvi;
typedef vector<string> vs;
typedef vector<vs> vvs;
typedef tuple<int, int, int, int> i4;
typedef map<i4, pair<vector<i4>, int>> bt_map;

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

void GenPerms(vn& v, int k, vvi& S, vvvi& P)
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
    vn c;
    for (newick_child* child = root->child; child; child = child->next)
        c.push_back(child->node);

    vn nc;
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

int Distance(vvvvi& DP, bt_map& P, vvi& T1, vvi& T2, vvs& L1, vvs& L2, int n, int m, int i, int j, int k, int l)
{
    if (i + j == n)
    {
        P[{i, j, k, l}] = {{}, -2};
        return m - k - l;
    }
    if (k + l == m)
    {
        P[{i, j, k, l}] = {{}, -1};
        return n - i - j;
    }

    int p = n - T1[i][j + 1] - j, q = m - T2[k][l + 1] - l;
    int x = i ? p + 1 : p + 2, y = i ? j : j - 1;
    int xx = k ? q + 1 : q + 2, yy = k ? l : l - 1;
    if (i == 0 && j == 0)
        x = 1, y = 0;
    if (k == 0 && l == 0)
        xx = 1, yy = 0;

    int o = DP[i][j + 1][k][l] + 1;
    int a = DP[i][j][k][l + 1] + 1;
    int b = DP[i][j + T1[i][j + 1]][k][l + T2[k][l + 1]];
    int c = DP[x][y][xx][yy];
    int d = L1[i][j + 1] != L2[k][l + 1];
    int mm = min(min(o, a), b + c + d);
    if (mm == b + c + d)
        P[{i, j, k, l}] = {{{x, y, xx, yy}, {i, j + T1[i][j + 1], k, l + T2[k][l + 1]}}, 0};
    else if (mm == o)
        P[{i, j, k, l}] = {{{i, j + 1, k, l}}, 1};
    else
        P[{i, j, k, l}] = {{{i, j, k, l + 1}}, 2};
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
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            if (newick_node* node = (i + j <= n ? GetLastRemoved(root, i, j) : nullptr))
                L[i][j] = node->taxon, T[i][j] = S[node];
}

int OrderedEditDist(bt_map& P, vvs& L1, vvs& L2, vvi& T1, vvi& T2, int n, int m)
{
    vvvvi DP(n + 1, vvvi(n + 1, vvi(m + 1, vi(m + 1, -1))));
    for (int i = n; i >= 0; --i)
        for (int j = n; j >= 0; --j)
            for (int k = m; k >= 0; --k)
                for (int l = m; l >= 0; --l)
                    if (i + j <= n && k + l <= m)
                        DP[i][j][k][l] = Distance(DP, P, T1, T2, L1, L2, n, m, i, j, k, l);

    return DP[0][0][0][0];
}

void DeallocTree(newick_node* root)
{
    for (newick_child* child = root->child; child; child = child->next)
        DeallocTree(child->node);
    delete root;
}

void PrintMatching(bt_map& BP, vvs& L1, vvs& L2, i4 r, bool swp)
{
    int i, j, k, l, op = BP[r].second;
    tie(i, j, k, l) = r;
    int n = L1.size() - 1, m = L2.size() - 1;
    switch (op)
    {
        case -2:
            for (int nl = l + 1; nl <= m - k; ++nl)
                cout << (!swp ? "INS " : "DEL ") << L2[k][nl] << '\n';
            break;
        case -1:
            for (int nj = j + 1; nj <= n - i; ++nj)
                cout << (!swp ? "DEL " : "INS ") << L1[i][nj] << '\n';
            break;
        case 0:
        {
            string x = L1[i][j + 1], y = L2[k][l + 1];
            cout << "MATCH " << (swp ? y : x) << ' ' << (swp ? x : y) << '\n';
            break;
        }
        case 1:
            cout << (!swp ? "DEL " : "INS ") << L1[i][j + 1] << '\n';
            break;
        case 2:
            cout << (!swp ? "INS " : "DEL ") << L2[k][l + 1] << '\n';
            break;
    }
    for (i4 p : BP[r].first)
        PrintMatching(BP, L1, L2, p, swp);
}

void EditDist(Graph& g1, Graph& g2)
{
    vn n1, n2;
    int k1 = GetNumTrees(g1.GetRoot(), n1), k2 = GetNumTrees(g2.GetRoot(), n2);
    Tree& t1 = (Tree&)(k1 >= k2 ? g1 : g2), &t2 = (Tree&)(k1 >= k2 ? g2 : g1);
    bool swp = false;
    if (k2 > k1)
        swap(n1, n2), swp = true;

    vvvi P;
    vvi S;
    GenPerms(n2, 0, S, P);
    bt_map BP;
    vvs BL1, BL2;
    int mm = numeric_limits<int>::max();
    for (int i = 0; i < P.size(); ++i)
    {
        newick_node* nroot = new newick_node(n2[0]->taxon);
        MakeTree(n2[0], nroot, 0, P[i]);
        int n = n1.size(), m = n2.size();
        vvs L1(n + 1, vs(n + 1)), L2(m + 1, vs(m + 1));
        vvi T1(n + 1, vi(n + 1, 1)), T2(m + 1, vi(m + 1, 1));
        GenTables(t1.GetRoot(), n, L1, T1);
        GenTables(t2.GetRoot(), m, L2, T2);
        bt_map NP;
        int nm = OrderedEditDist(NP, L1, L2, T1, T2, n, m);
        if (nm < mm)
            mm = nm, BP = NP, BL1 = L1, BL2 = L2;
        DeallocTree(nroot);
    }
    PrintMatching(BP, BL1, BL2, {0, 0, 0, 0}, swp);
    clog << "DIST: " << mm << endl;
}
