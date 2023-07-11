#ifndef H_HARD_FIXING
#define H_HARD_FIXING

#include <ilcplex/cplex.h>
#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include "cplex_model.h"
#include <time.h>


//================functions=================
void hard_fixing(Instance*, double, double);
void add_mip_start(Instance*, CPXENVptr, CPXLPptr);
void benders_loop2(Instance*, CPXENVptr, CPXLPptr);
#endif