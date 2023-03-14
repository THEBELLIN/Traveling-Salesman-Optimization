#ifndef VRP_H_ //avoid defining twice if included twice

#define VRP_H_

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

//===========constants=================
#define XSMALL 1e-5
#define EPSILON 1e-9
#define INF_DOUBLE 1e50
//#define VERBOSE 1			//  0, 10, 20, 30, 40, 50

//===============structs===============
typedef struct
{
	double x;
	double y;
}Point;

typedef struct
{
	int from;
	int to;
}edge;

typedef struct
{
	Point* points;
	int nnodes;
	double* cost;

	char inputfile[1000];
	char outputfile[1000];
	int randomseed;

	double tstart;
	double bestcost, tbest, bestlb;
	int* bestsol;
	int verbose;			//  0, 10, 20, 30, 40, 50

}Instance;

//=====================functions=======================================
void parse_TSPLIB(Instance*);
void print_error(const char*);
void initialize_instance(Instance*);
void parse_args(Instance*, int, char**);
void print_points(Instance*);
Point* generate_random_points(int);
Point* generate_random_points_range(int, double, double);
void free_instance(Instance*);
double distance(Point*, Point*);
void plot_generator(Instance*);

#endif
