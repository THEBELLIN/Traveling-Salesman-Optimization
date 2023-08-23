#ifndef ANNEALING_H
#define ANNEALING_H

#include "TSP.h"
#include "utility.h"
#include "heuristics.h"


//============ CONSTANTS ============
#define COOLING_FACTOR 0.99
#define INITIAL_TEMPERATURE 100.0

//============ FUNCTIONS ============
void simulated_annealing2(Instance*);
void random_two_opt_move2(int*,Instance*);

#endif
