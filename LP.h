#ifndef LP_H
#define LP_H

#include "Solver.h"
#include "Geno.h"

extern int c;

class LP : public Solver
{
public:
    LP(Graph& t1, Graph& t2, string d, double k, bool dag);
    ~LP();

    void Solve();
    void WriteSolution(string fileName);

    static int cf;
private:    
    template<class T>
    int Add()
    {
        int row_old = nr_rows;
        T c12(t1, t2, K, x, false);
        nr_rows += c12.AddTriplets(Triplets, nr_rows);
        T c21(t2, t1, K, x, true);
        nr_rows += c21.AddTriplets(Triplets, nr_rows);
        return nr_rows - row_old;
    }

    void MatchingConstraints();
    void SolveLP();

    int GetMax(newick_node* node, int& hmax);
    float SymdifDist(float weight);
    float JaccardDist(float weight);

    vector<ET> Triplets;
    Vector x, y;
    vvi K;
    Vector c;
    int nr_rows, nr_cols;
    int cnt;
};

#endif
