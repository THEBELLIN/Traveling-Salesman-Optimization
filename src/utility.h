#ifndef UTILITY_H
#define UTILITY_H

#include "TSP.h"

#define MALLOC(n, type) ((type *) malloc (n * sizeof(type)))
#define CALLOC(n, type) ((type *) calloc (n, sizeof(type)))
#define REALLOC(ptr, n, type) (realloc(ptr, n * sizeof(type)))
#define COST(i1, i2) (inst->cost[i1 * inst->nnodes + i2])

//==================inline functions===================================
inline int imax(int i1, int i2) { return (i1 > i2) ? i1 : i2; }
inline double dmin(double d1, double d2) { return (d1 < d2) ? d1 : d2; }
inline double dmax(double d1, double d2) { return (d1 > d2) ? d1 : d2; }

//====================functions============================
int comparePoints(const Point*, const Point*);
bool ccw(const Point*, const Point*, const Point*);
void swap(int*, int, int);
void initialize_cost(Instance*);
bool is_equal_double(double, double);
bool is_equal_points(Point*, Point*);
double rand01();
void invert_nodes(int*, int, int);

#endif
