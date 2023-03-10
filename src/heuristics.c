#include "TSP.h";
#include "heuristics.h";

Point* convex_hull(Instance* inst)
{
	if (inst->nnodes < 1)
		return NULL;

	Point* hull;
	hull = malloc(inst->nnodes * sizeof(Point));
}