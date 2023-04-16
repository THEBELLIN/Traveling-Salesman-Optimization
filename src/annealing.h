#ifndef ANNEALING_H
#define ANNEALING_H

#include "TSP.h"
#include "utility.h"
#include "heuristics.h"

void simulated_annealing(Instance*, int);
random_two_opt_move(int*,Instance*);
