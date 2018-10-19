/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "Timer.h"
#include "Greedy.h"
#include "LP.h"
#include "BnB.h"
#include "BnG.h"
#include "LPInt.h"
#include "LPCP.h"
#include "LPFInt.h"
#include "Similarity.h"
#include "Parallel.h"

#include <iostream>
using namespace std;
// global varialbe for cost matrix name 
std::string costMatrixFileName;
int Solver::cf;
bool Solver::tt;

Solver* MakeSolver(Graph& t1, Graph& t2, int argc, char** argv)
{
    int s = stoi(argv[argc - 1]);
    Solver::cf = stoi(argv[4 + (argc == 9) + 2 * (argc == 12)]);
    Solver::tt = argc == 12;
    string d = argv[5 + (argc == 9) + 2 * (argc == 12)];
    double k = stod(argv[6 + (argc == 9) + 2 * (argc == 12)]);
    double c = (s != 2) ? 0 : stod(argv[8 + 2 * (argc == 12)]);
    var_eps = (argc == 9) ? 0 : stod(argv[7 + 2 * (argc == 12)]);

    assert(LP::cf >= 0 && LP::cf <= 2);
    assert(d == "j" || d == "s" || d == "e");
    assert(s >= 0 && s <= 9);

    if (s == 0)
        return new Greedy(t1, t2, d, k, argc == 9);
    else if (s == 1)
        return new LP(t1, t2, d, k, argc == 9);
    else if (s == 2)
        return new BnB(t1, t2, d, k, argc == 9, c);
    else if (s == 3)
        return new LPCP(t1, t2, d, k, argc == 9);
    else if (s == 4)
        return new TestBnBSolver(t1, t2, d, k, argc == 9);
    else if (s == 5)
        return new LPInt(t1, t2, d, k, argc == 9);
		else if (s == 6) 
    	return new LPFInt(t1, t2, d, k, argc == 9);
		else if (s == 7)
			return new BFBnBSolver(t1, t2, d, k, argc == 9);
		else if (s == 8)
			return new DFBnBSolver(t1, t2, d, k, argc == 9);
		else // if (s == 9)
			return new HybridBnBSolver(t1, t2, d, k, argc == 9);
}

Graph* MakeDAG(const char* f1, const char* f2, int s)
{
    return ((s == 0) ? new DAG(f1, f2) : new LDAG(f1, f2))->Init();
}

pair<Graph*, Graph*> MakeGraphs(int argc, char** argv)
{
    if (argc == 10)
        return {new Tree(argv[1]), new Tree(argv[2])};
    else if (argc == 12)
        return {new Tree(argv[1], argv[2]), new Tree(argv[3], argv[4])};

    int s = stoi(argv[argc - 1]);
    return {MakeDAG(argv[1], argv[2], s), MakeDAG(argv[3], nullptr, s)};
}

void runHali(int argc, char** argv){
    
    for (int i=0; i < argc; ++i){
        cout << argv[i] << " ";
        if (i == 4)
            cout << costMatrixFileName << " ";
    }
    cout << "\n";
    
    Timer T;
    T.start();
    Graph *t1, *t2;
    tie(t1, t2) = MakeGraphs(argc, argv);
    if (stoi(argv[argc - 1]) > 9)
    {
        int s = stoi(argv[argc - 1]);
        Solver::cf = stoi(argv[4 + (argc == 9) + 2 * (argc == 12)]);
        Solver::tt = argc == 12;
        string d = argv[5 + (argc == 9) + 2 * (argc == 12)];
        double k = stod(argv[6 + (argc == 9) + 2 * (argc == 12)]);
        double c = (s != 2) ? 0 : stod(argv[8 + 2 * (argc == 12)]);
        var_eps = (argc == 9) ? 0 : stod(argv[7 + 2 * (argc == 12)]);
//         cout << "First: " << var_eps << endl;
        CSVReader reader(costMatrixFileName);
        std::vector<std::vector<double>> cost_matrix = reader.getDoubleData();
        int n = (int) cost_matrix.size();
        int m = (int) cost_matrix[0].size();
        std::vector<double> weight_list;
        for (int i=0; i < n-1; ++i)
        {
            for (int j=0; j < m-1; ++j){
                if (cost_matrix[i][m-1] + cost_matrix[n-1][j] - cost_matrix[i][j] > 0.0)
                    weight_list.push_back(cost_matrix[i][m-1] + cost_matrix[n-1][j] - cost_matrix[i][j]);
            }
        }
        std::sort(weight_list.begin(), weight_list.end());
        int id = (int) (var_eps*((double) weight_list.size()));
//         cout << "id:"<< id << endl;
        if (id > 0 && id < (int) weight_list.size())
            var_eps = std::max(0.0 ,weight_list.at(id));
        else 
            var_eps = std::max(0.0, weight_list.at(0));
        
        var_eps = 0.0;
        cout << "var_eps: " << var_eps << endl;
        assert(LP::cf >= 0 && LP::cf <= 2);
        assert(d == "j" || d == "s" || d == "e");
        double optVal = 0.0;
        std::vector<std::pair<int, int>> hali_sol = ParallelSolver(*t1, *t2, d, k, argc == 9, s - 9).Solve(argv[3 + (argc == 9) + 2 * (argc == 12)]);
        int false_align = 0;
        for (auto &u: hali_sol){
            if (u.first != u.second){
                false_align++;
            }
        }
        cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> OPT: "<< optVal <<" Total number of false alignments: " << false_align << endl;
        cout << "=============================================================================\n";
    }		
    else
    {
        Solver* solver = MakeSolver(*t1, *t2, argc, argv);
        solver->Solve(argv[3 + (argc == 9) + 2 * (argc == 12)]);
        delete solver;
        delete t1;
        delete t2;
    } 
    T.stop();
    cout << "TIME: " << T.secs() << " secs" << endl;
}
int main(int uargc, char** uargv)
{
    if (uargc != 2)
    {
        cout << "usage: " << uargv[0] << " <up_to_instance>" << endl;
        return EXIT_FAILURE;
    }
    int max_instances = stoi(uargv[1]);
    
    int argc = 12;
    char* argv[12] = {"./hali", "t1.tree", "t1.map", "t2.tree", "t2.map", "t_output", "2", "e", "0", "0.2", "0.01", "35"};
    Timer T;
    T.start();
    
    // change the parser for cell hali inputs. 
    for (int i = 0; i < max_instances; ++i)
    {
        for (int j=i+1; j < max_instances; ++ j)
        {
            cout << "Compute Hali for tree:" << i << " and " << "tree:"<<j<< endl; 
            //inputs
            std::string t1Name = "t_"+to_string(i)+".tree";
            std::string t2Name = "t_"+to_string(j)+".tree";
            std::string t1Map  = "t_"+to_string(i)+".map";
            std::string t2Map  = "t_"+to_string(j)+".map";
            std::string dtwCostName = "dtw_cost_matrix_"+to_string(i)+"_"+to_string(j)+".csv";
            std::string maxCostName = "max_cost_matrix_"+to_string(i)+"_"+to_string(j)+".csv";
            std::string averCostName = "aver_cost_matrix_"+to_string(i)+"_"+to_string(j)+".csv";
            // outputs
            std::string dtwOut = "dtw_output_"+to_string(i)+"_"+to_string(j)+".csv";
            std::string maxOut = "max_output_"+to_string(i)+"_"+to_string(j)+".csv";
            std::string averOut = "aver_output_"+to_string(i)+"_"+to_string(j)+".csv";
            
            argv[1] = const_cast<char*>(t1Name.c_str());
            argv[2] = const_cast<char*>(t1Map.c_str());
            argv[3] = const_cast<char*>(t2Name.c_str());
            argv[4] = const_cast<char*>(t2Map.c_str());
            cout << "DTW:\n";
            argv[5] = const_cast<char*>(dtwOut.c_str());
            costMatrixFileName = dtwCostName;
            runHali(argc, argv);
            
            cout << "Average:\n";
            argv[5] = const_cast<char*>(averOut.c_str());
            costMatrixFileName = averCostName;
            runHali(argc, argv);
            
            cout << "Maximum:\n";
            argv[5] = const_cast<char*>(maxOut.c_str());
            costMatrixFileName = maxCostName;
            runHali(argc, argv);
            cout << "===========================================***************************************==================================\n\n";
        }
    }
    T.stop();
    cout << "TIME: " << T.secs() << " secs" << endl;
}
