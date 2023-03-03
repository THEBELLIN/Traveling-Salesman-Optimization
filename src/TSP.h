#ifndef VRP_H_ //avoid defining twice if included twice

#define VRP_H_

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
	double x, y;
}Point;

typedef struct
{
	Point* points;
	int nnodes;
}Instance;

#endif
