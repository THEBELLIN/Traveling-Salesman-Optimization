#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "TSP.h"

#define NORMAL {.opt=NORM, .p1=1, .p2=0};

//option to select the method for extra mileage algorithm  start
typedef enum
{
	RAND,		//select 2 nodes at random for the strating set
	MAX_DIST	//select the 2 most distant nodes for the starting set
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
	int p1;
	int p2;
}em_options;

Point* convex_hull(Instance*);
void add_in_position(int, int, int*, int);
void extra_mileage(Instance*, int, em_options*);
void extra_mileage_det(Instance*, int);

#endif