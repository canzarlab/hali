#include "LPFInt.h"

LPFInt::LPFInt(Graph& t1, Graph& t2, string d, double k, bool dag) : LPInt(t1, t2, d, k, dag)
{
}

void LPFInt::Solve(string filename)
{
    LP::Solve(filename);
    LPInt::Solve(filename);
}
