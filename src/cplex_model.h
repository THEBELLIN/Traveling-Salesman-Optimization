#ifndef H_CPLEX_MODEL
#define H_CPLEX_MODEL

#include "cplex.h"
#include "TSP.h"
#include "utility.h"
#include "heuristics.h"

//================functions=================
int TSPopt(Instance*);
int xpos(int, int, Instance*);
void build_model(Instance*, CPXENVptr, CPXLPptr);
void build_sol(const double*, Instance*, int*, int*, int*);
double mip_value(CPXENVptr, CPXLPptr);
void mip_delete_all_mipstarts(CPXENVptr, CPXLPptr);
int mip_solution_available(CPXENVptr, CPXLPptr);
int mip_solved_to_optimality(CPXENVptr, CPXLPptr);
int mip_infeasible(CPXENVptr, CPXLPptr);
void benders_loop(Instance*, CPXENVptr, CPXLPptr);
void add_SEC(Instance*, const int, int*, CPXENVptr, CPXLPptr);
#endif
