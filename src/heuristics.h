#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "TSP.h"
#include "utility.h"
#include "genetics.h"

//==============SHORTCUTS=================
#define NORMAL {.opt=NORM, .p1=1, .p2=0};
#define GRASP2(p) {.opt=GRASP_2, .p1=p, .p2=0};
#define GRASP3(p, pp) {.opt=GRASP_3, .p1=p, .p2=pp};

//=================ENUMS===========================


//=================STRUCTS============================

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
int nearest_neighbor_allstart(Instance*);
Point* convex_hull(Point*, int, int*);
void nearest_neighbor_grasp_random(Instance*, int);
void next_bestsol(Instance*, int);
void tabu_search(Instance*); 
void vns(Instance*);
void two_opt(Instance*);
void kick(Instance*, int*);

//TODO SPLIT IN DIFFERENT FILES

#endif