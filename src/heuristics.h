#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "TSP.h"

Point* convex_hull(Instance*);
void nearest_neighbor(Instance*, int);
int nearest_neighbor_allstart(Instance*);
void nearest_neighbor_grasp(Instance*, int, double, double);

#endif