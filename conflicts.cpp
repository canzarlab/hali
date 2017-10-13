#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include <utility>
#include <tuple>
#include "newick.h"
#include "Similarity.h"
using namespace std;

typedef vector<int> vi;
typedef pair<int, int> ii;
typedef vector<vb> vbb;

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

class graph
{
protected:
    vector<newick_node*> nodes;
    newick_node* root;
    int size, ntax;
public:
    vector<newick_node*> leaves;
    vi matches;

    graph(newick_node* root) :
        root(root),
        size(0),
        ntax(1)
    {
        init1(root);
        nodes.resize(size);
        matches.resize(size, -1);
    }

    virtual void init()
    {
        leaves.reserve(ntax - 1);
        root->taxoni = 0;
        //root->taxon = "0";
        init2(root, 0);
    }

    virtual ~graph()
    {
        dealloc_dag(root, size);
    }

    int get_size() const
    {
        return size;
    }

protected:
    virtual void pre(newick_node* node, int depth) = 0;
    virtual void post(newick_node* node, newick_node* child, int depth) = 0;

private:
    void init1(newick_node* node)
    {
        size++;
        if (!node->child)
            ntax++;

        for (newick_child* child = node->child; child; child = child->next)
        {
            newick_node* cnode = child->node;
            newick_parent** parentptr = &cnode->parent;
            if (!cnode->parent)
                init1(cnode);

            while (*parentptr) parentptr = &(*parentptr)->next;
            *parentptr = new newick_parent(node);
        }
    }

    void init2(newick_node* node, int depth)
    {
        if (node->taxoni >= 0 && nodes[node->taxoni])
            return;

        if (node->taxoni == -1)
        {
            nodes[node->taxoni = ntax++] = node;
            //node->taxon = to_string(node->taxoni);
        }
        else
            nodes[node->taxoni] = node;

        if (!node->child)
            leaves.push_back(node);

        pre(node, depth);
        for (newick_child* child = node->child; child; child = child->next)
        {
            init2(child->node, depth + 1);
            post(node, child->node, depth);
        }
    }
};

class tree : public graph
{
    segtree* seg;
    vi L, E, H;
    int idx;
public:
    tree(newick_node* root) : graph(root), idx(0)
    {
    }

    ~tree()
    {
        delete seg;
    }

    virtual void init()
    {
        L.resize(2 * size);
        E.resize(2 * size);
        H.resize(size, -1);
        graph::init();
        seg = new segtree(L);
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

protected:
    virtual void pre(newick_node* node, int depth)
    {
        H[node->taxoni] = idx;
        E[idx] = node->taxoni;
        L[idx++] = depth;
    }

    virtual void post(newick_node* node, newick_node* child, int depth)
    {
        E[idx] = node->taxoni;
        L[idx++] = depth;
    }
};

class dag : public graph
{
    vector<set<int> > D;
public:
    dag(newick_node* root) : graph(root), D(size)
    {
    }

    virtual void init()
    {
        ntax = 1;
        graph::init();
    }

    bool dsc(int p, int q)
    {
        return D[q].find(p) != D[q].end();
    }

protected:
    virtual void pre(newick_node* node, int depth)
    {
    }

    virtual void post(newick_node* node, newick_node* child, int depth)
    {
        set<int>& cr = D[child->taxoni];
        D[node->taxoni].insert(cr.begin(), cr.end());
        D[node->taxoni].insert(child->taxoni);
    }
};

class cc
{
protected:
    graph *t1, *t2;
    int c1, c2;
public:
    cc(graph* t1, graph* t2) : t1(t1), t2(t2)
    {
    }

    virtual ~cc()
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
        vbb C(t1->get_size(), vb(t1->get_size()));
        for (newick_node* leaf : t1->leaves)
            dfs_c_i(leaf, leaf, C);
        swap(t1, t2);
        return make_pair(c1, c2);
    }

    void dfs_c_i(newick_node* node, newick_node* rnode, vbb& C)
    {
        if (node != rnode)
            check_c(node->taxoni, rnode->taxoni);

        C[node->taxoni][rnode->taxoni] = true;
        for (newick_parent* parent = node->parent; parent; parent = parent->next)
        {
            newick_node* pn = parent->node;
            if (!C[pn->taxoni][rnode->taxoni])
                dfs_c_i(pn, rnode, C);
            if (!C[pn->taxoni][pn->taxoni])
                dfs_c_i(pn, pn, C);
        }
    }

    virtual void check_c(int p, int q) = 0;
};

class cct : public cc
{
public:
    // TODO: n should probaby be number of aligns in file so we do a batch comparison instead of this
    cct(const char* fn1, const char* fn2, const char* fn3, int n) :
        cc(new tree(load_tree(fn1)), new tree(load_tree(fn2)))
    {
        t1->init();
        t2->init();
        string in;
        ifstream f(fn3);
        for (int i = 0; i < n; ++i)
            while (getline(f, in) && in != "-----------------");

        for (int i = 0; i < 5; ++i)
            getline(f, in);

        tree* t1 = (tree*)this->t1;
        tree* t2 = (tree*)this->t2;
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

protected:
    virtual void check_c(int p, int q)
    {
        int mp = t1->matches[p];
        int mq = t1->matches[q];
        if (mp == -1 || mq == -1)
            return;

        tree* t1 = (tree*)this->t1;
        tree* t2 = (tree*)this->t2;
        int lp = t1->lca(p, q);
        int rp = t2->lca(mp, mq);
        if (rp != mp && rp != mq) c2++;
        else if (t1->matches[lp] != rp) c1++;
    }

private:
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

class ccd : public cc
{
    mnls cl1, cl2;
    msn A, B;
public:
    ccd(const char* fn1, const char* fn2, const char* fn3, const char* fn4) :
        cc(nullptr, nullptr)
    {
        (t1 = new dag(load_dag(fn1, fn2, cl1, A)))->init();
        (t2 = new dag(load_dag(fn3, nullptr, cl2, B)))->init();
        ifstream f(fn4);
        string a, b;
        double p, w = 0;
        while (getline(f, a), a != "Matched");
        while (f >> a >> b >> p)
        {
            newick_node* ln = A[a];
            newick_node* rn = B[b];
            int l = ln->taxoni;
            int r = rn->taxoni;
            t1->matches[l] = r;
            t2->matches[r] = l;
            w += JaccardSim(cl1[ln], cl2[rn], 1);
        }
        cout << "Weight of matching: " << w << endl;
    }

    virtual void check_c(int p, int q)
    {
        int mp = t1->matches[p];
        int mq = t1->matches[q];
        if (mp == -1 || mq == -1)
            return;

        dag* t1 = (dag*)this->t1;
        dag* t2 = (dag*)this->t2;
        if (!t2->dsc(mp, mq) && !t2->dsc(mq, mp))
            c2++;
        else if (t1->dsc(q, p) == t2->dsc(mp, mq))
            c1++;
    }
};

int main(int argc, char** argv)
{
    int c1, c2, opt = (argc == 6 ? stoi(argv[5]) : -1);
    if (opt > 1 || opt < 0)
    {
        cout << "tree usage " << argv[0] << " <t1.newick> <t2.newick> <align> <n> 0" << endl;
        cout << "dag usage: " << argv[0] << " <yeastnet> <mapping> <go> <align> 1" << endl;
        return EXIT_FAILURE;
    }
    else if (opt == 0)
        tie(c1, c2) = cct(argv[1], argv[2], argv[3], stoi(argv[4])).get_c();
    else
        tie(c1, c2) = ccd(argv[1], argv[2], argv[3], argv[4]).get_c();
    cout << "C1 violations: " << c1 << endl;
    cout << "C2 violations: " << c2 << endl;
}
