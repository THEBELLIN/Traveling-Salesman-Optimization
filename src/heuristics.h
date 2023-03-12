#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "TSP.h"

Point* convex_hull(Instance*);
void extra_mileage(Instance*, int);

//option to select the method for extra mileage algorithm  start
typedef enum
{
	RAND,		//select 2 nodes at random for the strating set
	MAX_DIST	//select the 2 most distant nodes for the starting set
}em_start;

//option to select the method for extra mileage algorithm 
typedef enum
{
	RAND,		//select 2 nodes at random for the strating set
	MAX_DIST	//select the 2 most distant nodes for the starting set
}em_opt;

#endif