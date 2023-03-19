#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "TSP.h"
#include "utility.h"

#define NORMAL {.opt=NORM, .p1=1, .p2=0};
#define GRASP2(p) {.opt=GRASP_2, .p1=p, .p2=0};
#define GRASP3(p, pp) {.opt=GRASP_3, .p1=p, .p2=pp};

//option to select the method for extra mileage algorithm  start
typedef enum 
{
	RAND,		//select 2 nodes at random for the strating set
	MAX_DIST,	//select the 2 most distant nodes for the starting set
	CONV_HULL	//start with the convex hull
}em_start;

//option to select the method for extra mileage algorithm 
typedef enum
{
	NORM,		//always select the best option at each step -> deterministic
	GRASP_2,	//use the second best option with probability 1-p
	GRASP_3		//use the 1st option with probability p1, the 2nd with p2, 3rd with 1-p1-p2
}em_opt;

typedef struct
{
	em_opt opt;
	double p1;
	double p2;
}em_options;

int* points_to_indexes(Instance*, Point*, int);
void add_in_position(int, int, int*, int);
void extra_mileage(Instance*, int, em_options*);
void extra_mileage_det(Instance*, int);
void extra_mileage_grasp2(Instance*, em_start, double);
void extra_mileage_grasp3(Instance*, em_start, double, double);
void nearest_neighbor_grasp2(Instance*, int, double, double);
void nearest_neighbor_grasp(Instance*, int, double);
int nearest_neighbor_allstart(Instance*);
void nearest_neighbor(Instance*, int);
Point* convex_hull(Point*, int, int*);
void nearest_neighbor_grasp_random(Instance*, int, double);
void two_opt_move(Instance*);
#endif