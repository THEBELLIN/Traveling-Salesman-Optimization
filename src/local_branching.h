#ifndef H_LOCAL_BRANCHING
#define H_LOCAL_BRANCHING

#include <cplex.h>
#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include "cplex_model.h"
#include "hard_fixing.h"
#include <time.h>


//================functions=================
void local_branching(Instance*);
#endif