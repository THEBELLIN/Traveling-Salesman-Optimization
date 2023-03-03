#ifndef VRP_H_ //avoid defining twice if included twice

#define VRP_H_

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>


#define XSMALL 1e-5
#define EPSILON 1e-9
#define VERBOSE 1			//  0, 10, 20, 30, 40, 50

//TODO L implement basic stuff, basic system, data structure, reading instance TSPLIB,  parse command line, 
//TODO M plot instance, generate random data
//DONE
typedef struct
{
	double x, y;
}Point;

typedef struct
{
	Point* points;
	int nnodes;

	char inputfile[1000];
	char outputfile[1000];
	int randomseed;

	double tstart;
	double zbest, tbest, bestlb;
	int* bestsol;

}Instance;

inline int imax(int i1, int i2) { return (i1 > i2) ? i1 : i2; }
inline double dmin(double d1, double d2) { return (d1 < d2) ? d1 : d2; }
inline double dmax(double d1, double d2) { return (d1 > d2) ? d1 : d2; }

#endif
