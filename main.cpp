#include "PhylogeneticTree.h"
#include "young/solver.h"
#include <fstream>

int main(int argc, char** argv)
{
    if (argc != 5)
    {
        cout << "usage: " << argv[0] << " <filename.newick> <filename.newick> <filename.sim> <epsilon>" << endl;
        return EXIT_FAILURE;
    }

    PhylogeneticTree t1(argv[1]);
    PhylogeneticTree t2(argv[2]);
    int n = t1.getNumNodes(), m = t2.getNumNodes();

    ifstream SimFile(argv[3]);
    double* C = new double[n * m];
    for (int i = 0; i < n * m; ++i)
        SimFile >> C[i];

    double epsilon = stod(argv[4]);
    Solver* solver = Solver::create(epsilon);

    int k = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            if (C[i * n + j] != 0)
                solver->add_entry(i, m * i + j - k, 1. / C[i * n + j]);
            else
                ++k;

    k = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            if (C[i * n + j] != 0)
                solver->add_entry(n + j, i * m + j - k, 1. / C[i * n + j]);
            else
                ++k;
    delete[] C;

    solver->done_adding_entries();

    if (solver->solve())
    {
        double col_value = 0, row_value = 0;
        for (int row = 0; row < solver->n_rows();  ++row) row_value += solver->value_of_row_variable(row);
        for (int col = 0; col < solver->n_cols();  ++col) col_value += solver->value_of_col_variable(col);
        cout << row_value << " " << col_value << endl;
    }
    else
        cout << "Not solved" << endl;

    delete solver;
}
