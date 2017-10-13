#ifndef SIMILARITY_H
#define SIMILARITY_H

#include <list>
#include <string>
using namespace std;

typedef list<string> ls;

double JaccardSim(const ls& L1, const ls& L2, double k);
double SymdifSim(const ls& L1, const ls& L2);

#endif
