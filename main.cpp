#include "PhylogeneticTree.h"
#include "young/solver.h"
#include <iostream>
#include <fstream>
#include <sstream>

class LP
{
    typedef vector<vector<double> > DPT;
public:
    LP(const char* p1, const char* p2, double e) : t1(p1), t2(p2), solver(Solver::create(e))
    {
        DP.resize(t1.GetNumNodes());
        for (auto& v : DP)
            v.resize(t2.GetNumNodes());
    }

    ~LP()
    {
        delete solver;
    }

    void GenMatchingConstraints(const char* path)
    {
        ifstream SimFile(path);
        if (!SimFile)
        {
            cout << "Failed to open " << path << endl;
            exit(EXIT_FAILURE);
        }

        vector<vector<double> > C;
        string line;
        for (int i = 0, k = 0; getline(SimFile, line) && !line.empty(); ++i)
        {
            istringstream ss(line);
            if (!t1.NodeExists(i + 1))
                continue;

            t1.N[i] = k++;        
            C.push_back(vector<double>());
            double w;
            for (int j = 0, l = 0; ss >> w; ++j)
            {
                if (!t2.NodeExists(j + 1))
                    continue;

                t2.N[j] = l++;
                C.back().push_back(w);
            }
        }
        
        int n = t1.GetNumNodes(), m = t2.GetNumNodes(), k = 0;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                if (C[i][j] != 0)
                {
                    solver->add_entry(i, m * i + j - k, 1. / C[i][j]);
                    solver->add_entry(n + j, i * m + j - k, 1. / C[i][j]);
                }
                else
                    ++k;
            }
        }
    }

    void GenCrossingConstraints()
    {
        dfs1(t1.GetRoot());
    }
    
    void Solve()
    {    
        solver->done_adding_entries();
        if (!solver->solve())
            cout << "Not solved" << endl;
    }
    
private:
    double GetWeight(newick_node* nodel, newick_node* noder)
    {
        return 0;
    }

    void dfs1(newick_node* node)
    {
        dfs2(t2.GetRoot(), node);
        for (newick_child* child = node->child; child; child = child->next)
            dfs1(child->node);
    }
    
    double dfs2(newick_node* node, newick_node* nodel)
    {
        double mx = 0;
        for (newick_child* child = node->child; child; child = child->next)
            mx = max(mx, dfs2(child->node, nodel));
        return DP[t1.GetIndex(nodel)][t2.GetIndex(node)] = mx + GetWeight(nodel, node);
    }

    DPT DP;
    PhylogeneticTree t1, t2;
    Solver* solver;
};

int main(int argc, char** argv)
{
    if (argc != 5)
    {
        cout << "usage: " << argv[0] << " <filename.newick> <filename.newick> <filename.sim> <epsilon>" << endl;
        return EXIT_FAILURE;
    }

    LP lp(argv[1], argv[2], stod(argv[4]));
    lp.GenMatchingConstraints(argv[3]);
    lp.Solve();
    lp.GenCrossingConstraints();    
}
