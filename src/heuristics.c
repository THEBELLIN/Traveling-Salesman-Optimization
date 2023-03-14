#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include <assert.h>

Point* convex_hull(Instance* inst, int* hsize) //mostly found online
{
	if (inst->nnodes < 1)
		return NULL;

	Point* hull;
	int capacity = 4;
    int size = 0;
	hull = MALLOC(capacity, Point);

	qsort(inst->points, inst->nnodes, sizeof(Point), comparePoints);

    /* lower hull */
    for (int i = 0; i < inst->nnodes; ++i) 
    {
        while (size >= 2 && !ccw(&hull[size - 2], &hull[size - 1], &inst->points[i]))
            --size;
        if (size == capacity) 
        {
            capacity *= 2;
            hull = REALLOC(hull, capacity, Point);
        }
        assert(size >= 0 && size < capacity);
        hull[size++] = inst->points[i];
    }

    /* upper hull */
    int t = size + 1;
    for (int i = inst->nnodes - 1; i >= 0; i--) {
        while (size >= t && !ccw(&hull[size - 2], &hull[size - 1], &inst->points[i]))
            --size;
        if (size == capacity) {
            capacity *= 2;
            hull = REALLOC(hull, capacity, Point);
        }
        assert(size >= 0 && size < capacity);
        hull[size++] = inst->points[i];
    }
    --size;
    assert(size >= 0);
    hull = xrealloc(hull, size * sizeof(Point));
    *hsize = size;
    return hull;
}





