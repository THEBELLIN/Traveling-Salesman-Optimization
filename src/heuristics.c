#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include <assert.h>

//returns the convex hull of points in instance
//NOTE: THIS REARRANGES POINTS TO BE SORTED -> INDEXES CHANGE
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
    hull = REALLOC(hull, size, Point);
    *hsize = size;
    return hull;
}

//doesn't check for memory access overflow
void add_in_position(int h, int pos, int* arr, int len)
{
    if (h != len)
        swap(arr, h, len);
    //make space
    for(int i = len; i > pos; i--)
        swap(arr, i, i - 1);
}

void extra_mileage(Instance* inst, em_start start, em_options* options)
{
    if (options->opt == NORM)
        extra_mileage_det(inst, start);
    else
        print_error("Wrong parameter for extra mileage algorithm");
}

void extra_mileage_det(Instance* inst, em_start start)
{
    //construct cost matrix if not already done
    if (inst->cost == NULL)
        print_error("Cost matrix not initialized");

    //select startind indexes
    int i1, i2;
    int n = inst->nnodes;

    if (start == RAND)
    {
        i1 = rand() % n;
        i2 = i1;
        while(i2==i1)
            i2 = rand() % n;
    }
    else if (start == MAX_DIST)
    {
        int maxcost = 0;
        i1 = 0;
        i2 = 0;
        for (int i = 1; i < (n * n); i++)
        {
            if (inst->cost[i] > maxcost)
            {
                maxcost = inst->cost[i];
                i1 = i / n;
                i2 = i % n;
            }
        }
    }
    else
        print_error("Wrong starting flag in extra mileage");

    //extra mileage loop
    double curr_cost = 0;
    int current_nodes = 2;
    for (int i = 0; i < n; i++)
        inst->bestsol[i] = i;
    inst->bestsol[n] = i1;
    //set current tour
    swap(inst->bestsol, 0, i1);
    swap(inst->bestsol, 1, i2);
    swap(inst->bestsol, 2, n);
    curr_cost += 2 * COST(i1, i2);
    //untill all nodes are added
    while (current_nodes < n)
    {
        double new_cost = curr_cost;
        double min_new_cost = INF_DOUBLE;
        int place;
        int new_point_index;
        //for every edge in current tour
        for (int i = 0; i < current_nodes; i++) //considering edge inst->bestsol[i]->inst->bestsol[i+1]
        {
            edge e;
            e.from = inst->bestsol[i];
            e.to = inst->bestsol[i + 1];
            //for all points not already considered
            for (int j = current_nodes + 1; j < n + 1; j++) //considering point inst->bestsol[j]
            {
                new_cost = curr_cost - COST(e.from, e.to) + COST(e.from, inst->bestsol[j]) + COST(inst->bestsol[j], e.to);
                if (new_cost < min_new_cost)
                {
                    min_new_cost = new_cost;
                    new_point_index = j;
                    place = i + 1;
                }
            }
        }
        curr_cost = min_new_cost;
        current_nodes++;
        add_in_position(new_point_index, place, inst->bestsol, current_nodes+1);
        //DEBUG
        /*for (int i = 0; i < current_nodes + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf(" || ");
        for (int i = current_nodes+1; i < n + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf("\n");*/
    }

    //save solution
    inst->bestcost = curr_cost;
}