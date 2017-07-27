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
public:
    class lca
    {
        vi L, E, H;
        int idx;
        tree* t;
        segtree* seg;
    public:
        lca(tree* t) :
            L(2 * t->size()),
            E(2 * t->size()),
            H(t->size(), -1),
            idx(0),
            t(t)
        {
            init(0, 0);
            seg = new segtree(L);
        }

        ~lca()
        {
            delete seg;
        }

        int query(int p, int q)
        {
            if (H[q] < H[p]) swap(p, q);
            return E[seg->rmq(H[p], H[q])];
        }

        int query(const vi& v, int i = 0)
        {
            if (v.size() - i == 1)
                return v[i];
            return query(v[i], query(v, i + 1));
        }

    private:
        void init(int cur, int depth)
        {
            H[cur] = idx;
            E[idx] = cur;
            L[idx++] = depth;
            for (newick_child* child = t->get(cur)->child; child; child = child->next)
            {
                init(child->node->taxoni, depth + 1);
                E[idx] = cur;
                L[idx++] = depth;
            }
        }
    };

    tree(newick_node* root) : root(root), sz(0), ntax(1)
    {
        init1(root);
        nodes.resize(size());
        matches.resize(size(), -1);
        leaves.reserve(ntax - 1);
        root->taxoni = 0;
        init2(root);
        tlca = new lca(this);
    }

    ~tree()
    {
        delete tlca;
        delete root;
    }

    lca* get_lca()
    {
        return tlca;
    }

    void match(int l, int r)
    {
        matches[l] = r;
    }

    int get_match(int p)
    {
        return matches[p];
    }

    int size()
    {
        return sz;
    }

    int n_leaves()
    {
        return leaves.size();
    }

    newick_node* get(int i)
    {
        return nodes[i];
    }

private:
    void init1(newick_node* node)
    {
        sz++;
        if (!node->child)
            ntax++;

        for (newick_child* child = node->child; child; child = child->next)
        {
            child->node->parent = node;
            init1(child->node);
        }
    }

    void init2(newick_node* node)
    {
        if (node->taxoni != -1)
            nodes[node->taxoni] = node;
        else
            nodes[node->taxoni = ntax++] = node;

        node->taxon = to_string(node->taxoni);
        if (!node->child)
            leaves.push_back(node);

        for (newick_child* child = node->child; child; child = child->next)
            init2(child->node);
    }

    newick_node* root;
    vector<newick_node*> nodes;
    vector<newick_node*> leaves;
    vector<int> matches;
    int sz, ntax;
    lca* tlca;
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
            int l = t1->get_lca()->query(get_cluster(f));
            int r = t2->get_lca()->query(get_cluster(f));
            t1->match(l, r);
            t2->match(r, l);
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
        swap(t1, t2);
        tie(d1, d2) = get_c_i();
        swap(t1, t2);
        return make_pair(c1 + d1, c2 + d2);
    }

private:
    ii get_c_i()
    {
        c1 = c2 = 0;
        for (int i = 1; i < t1->n_leaves(); ++i)
            for (newick_node* node = t1->get(i); node->parent; node = node->parent)
                for (newick_node* nodeup = node->parent; nodeup; nodeup = nodeup->parent)
                    check_c(node->taxoni, nodeup->taxoni);
        return make_pair(c1, c2);
    }

    void check_c(int p, int q)
    {
        int mp = t1->get_match(p);
        int mq = t1->get_match(q);
        if (mp == -1 || mq == -1)
            return;

        int lp = t1->get_lca()->query(p, q);
        int rp = t2->get_lca()->query(mp, mq);
        if (rp != mp && rp != mq) c2++;
        else if (t1->get_match(lp) != rp) c1++;
    }

    vi get_cluster(ifstream& file)
    {
        while (file.peek() != '(')
            file.ignore();
        file.ignore();

        vi v;
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
