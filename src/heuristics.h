#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "TSP.h"
#include "utility.h"
#include "genetics.h"

//==============SHORTCUTS=================
#define NORMAL {.opt=NORM, .p1=1, .p2=0};
#define GRASP2(p) {.opt=GRASP_2, .p1=p, .p2=0};
#define GRASP3(p, pp) {.opt=GRASP_3, .p1=p, .p2=pp};

//================= CONSTANTS =======================
#define COOLING_FACTOR 0.6
#define INITIAL_TEMPERATURE 70.0


//====================FUNCTIONS=======================
int* points_to_indexes(Instance*, Point*, int, int*);
void add_in_position(int, int, int*, int);
void extra_mileage(Instance*);
void extra_mileage_det(Instance*);
void extra_mileage_grasp2(Instance*);
void extra_mileage_grasp3(Instance*);
void nearest_neighbor(Instance*);
void nearest_neighbor_det(Instance*, int);
void nearest_neighbor_grasp2(Instance*, int);
void nearest_neighbor_grasp3(Instance*, int);
void grasp_iter(Instance*);
void grasp_extra_iter(Instance*);
int nearest_neighbor_allstart(Instance*);
Point* convex_hull(Point*, int, int*);
void nearest_neighbor_grasp_random(Instance*, int);
void next_bestsol(Instance*, int);
void tabu_search(Instance*);
void two_opt(Instance*);
void vns(Instance*);
void kick(Instance*, int*);
int checkNonEqual(int*, int, int);
double random_two_opt_move(Instance*, int*);
void simulated_annealing(Instance*);


//TODO SPLIT IN DIFFERENT FILES

#endif