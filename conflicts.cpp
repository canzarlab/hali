#include <iostream>
#include <vector>
#include <fstream>
#include <utility>
#include <tuple>
#include "newick.h"
using namespace std;

typedef vector<int> vi;
typedef pair<int, int> ii;

class segtree
{
    const vi& A;
    vi st;
    int n;
public:
    segtree(const vi& A) : A(A), n(A.size())
    {
        st.assign(4 * n, 0);
        build(1, 0, n - 1);
    }

    int rmq(int i, int j)
    {
        return rmq(1, 0, n - 1, i, j);
    }

private:
    int left(int p)
    {
        return p << 1;
    }

    int right(int p)
    {
        return (p << 1) + 1;
    }

    void build(int p, int L, int R)
    {
        if (L == R)
            st[p] = L;
        else
        {
            build(left(p), L, (L + R) / 2);
            build(right(p), (L + R) / 2 + 1, R);
            int p1 = st[left(p)], p2 = st[right(p)];
            st[p] = (A[p1] <= A[p2]) ? p1 : p2;
        }
    }

    int rmq(int p, int L, int R, int i, int j)
    {
        if (i > R || j < L) return -1;
        if (L >= i && R <= j) return st[p];

        int p1 = rmq(left(p), L, (L + R) / 2, i, j);
        int p2 = rmq(right(p), (L + R) / 2 + 1, R, i, j);

        if (p1 == -1) return p2;
        if (p2 == -1) return p1;
        return (A[p1] <= A[p2]) ? p1 : p2;
    }
};

class tree
{
    vector<newick_node*> nodes;
    newick_node* root;
    segtree* seg;
    vi L, E, H;
    int size, ntax, idx;
public:
    vector<newick_node*> leaves;
    vi matches;

    tree(newick_node* root) :
        root(root),
        size(0),
        ntax(1),
        idx(0)
    {
        init1(root);
        nodes.resize(size);
        matches.resize(size, -1);
        L.resize(2 * size);
        E.resize(2 * size);
        H.resize(size, -1);
        leaves.reserve(ntax - 1);
        root->taxoni = 0;
        root->taxon = "0";
        init2(root, 0);
        seg = new segtree(L);
    }

    ~tree()
    {
        delete seg;
        delete root;
    }

    int lca(int p, int q)
    {
        if (H[q] < H[p]) swap(p, q);
        return E[seg->rmq(H[p], H[q])];
    }

    int lca(const vi& v, int i = 0)
    {
        if (v.size() - i == 1)
            return v[i];
        return lca(v[i], lca(v, i + 1));
    }

private:
    void init1(newick_node* node)
    {
        size++;
        if (!node->child)
            ntax++;

        for (newick_child* child = node->child; child; child = child->next)
        {
            child->node->parent = node;
            init1(child->node);
        }
    }

    void init2(newick_node* node, int depth)
    {
        if (node->taxoni == -1)
        {
            nodes[node->taxoni = ntax++] = node;
            node->taxon = to_string(node->taxoni);
        }
        else
            nodes[node->taxoni] = node;

        if (!node->child)
            leaves.push_back(node);

        H[node->taxoni] = idx;
        E[idx] = node->taxoni;
        L[idx++] = depth;
        for (newick_child* child = node->child; child; child = child->next)
        {
            init2(child->node, depth + 1);
            E[idx] = node->taxoni;
            L[idx++] = depth;
        }
    }
};

class cc
{
    tree *t1, *t2;
    int c1, c2;
public:    
    // TODO: n should probaby be number of aligns in file so we do a batch comparison instead of this
    cc(const char* fn1, const char* fn2, const char* fn3, int n) :
        t1(new tree(load_tree(fn1))),
        t2(new tree(load_tree(fn2)))
    {
        string in;
        ifstream f(fn3);
        for (int i = 0; i < n; ++i)
            while (getline(f, in) && in != "-----------------");

        for (int i = 0; i < 5; ++i)
            getline(f, in);

        while (f && f.peek() != '-')
        {
            int l = t1->lca(get_cluster(f));
            int r = t2->lca(get_cluster(f));
            t1->matches[l] = r;
            t2->matches[r] = l;
            getline(f, in);
            getline(f, in);
        }
    }

    ~cc()
    {
        delete t1;
        delete t2;
    }

    ii get_c()
    {
        int c1, c2, d1, d2;
        tie(c1, c2) = get_c_i();
        tie(d1, d2) = get_c_i();
        return make_pair(c1 + d1, c2 + d2);
    }

private:
    ii get_c_i()
    {
        c1 = c2 = 0;
        for (newick_node* leaf : t1->leaves)
            for (newick_node* node = leaf; node->parent; node = node->parent)
                for (newick_node* nodeup = node->parent; nodeup; nodeup = nodeup->parent)
                    check_c(node->taxoni, nodeup->taxoni);
        swap(t1, t2);
        return make_pair(c1, c2);
    }

    void check_c(int p, int q)
    {
        int mp = t1->matches[p];
        int mq = t1->matches[q];
        if (mp == -1 || mq == -1)
            return;

        int lp = t1->lca(p, q);
        int rp = t2->lca(mp, mq);
        if (rp != mp && rp != mq) c2++;
        else if (t1->matches[lp] != rp) c1++;
    }

    vi get_cluster(ifstream& file)
    {
        vi v;
        while (file.get() != '(');
        while (file.peek() != ')')
        {
            int l;
            file >> l;
            v.push_back(l);
            if (file.peek() == ',')
                file.ignore();
        }
        return v;
    }
};

int main(int argc, char** argv)
{
    if (argc != 5)
    {
        cout << "usage " << argv[0] << " <t1.newick> <t2.newick> <align> <n>\n";
        return EXIT_FAILURE;
    }
    int c1, c2;
    tie(c1, c2) = cc(argv[1], argv[2], argv[3], stoi(argv[4])).get_c();
    cout << c1 << " " << c2 << '\n';
}
