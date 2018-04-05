/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef BNG_H
#define BNG_H

#include "LP.h"

class BnG : public LP
{
public:
    BnG(Graph& t1, Graph& t2, string d, double k, bool dag);
    virtual void Solve(string filename) override;

private:
    void  Cleanup(size_t nr_t, size_t nr_r);
    bool  SolveLP() override;
    float Geno();

    vector<bool> sys_x;
	Vector		 sys_lo;
	Vector		 sys_hi;
};

#endif
