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


//greedy NN given a starting point O(n^2)
void nearest_neighbor(Instance* inst, int start) {
    if (start < 0)
        print_error("Invalid choice of the start node");
    // allocates the array with the current tour and the unvisited nodes
    int* current = (int*)calloc(inst->nnodes, sizeof(int));
    int n = inst->nnodes;
    int len = 0;
    double tourCost = 0;
    int last = start;
    double min = MAX;
    int best_pos;
    // initializes the array of costs 
    initialize_cost(inst);
    //initialize current with the trivial permutation 0-1-2....-n
    for (int i = 0; i < inst->nnodes; i++) {
        current[i] = i;
    }
    // moves at the beginnig the first cover node which is the starting node
    swap(current, 0, start);
    len++;
    for (int i = 1; i < inst->nnodes; i++) {
        // finds minimum cost node from the last visited node of the tour
        for (int j = len; j < inst->nnodes; j++) {
            if (inst->cost[last * n + current[j]] < min) {
                min = inst->cost[last * n + current[j]];
                best_pos = j;
            }
        }
        // update the cost and rearranges the current array
        tourCost = tourCost + min;
        swap(current, len, best_pos);
        min = MAX;
        last = current[len];
        len++;
    }
    tourCost += inst->cost[start * n + last];
    inst->bestsol = current;
    inst->zbest = tourCost;
}

// greedy NN with O(n^3) finds the best starting point and its solution
int nearest_neighbor_allstart(Instance* inst) {
    double best_current_cost = MAX;
    int best_start;
    for (int i = 0; i < inst->nnodes; i++) {
        nearest_neighbor(inst, i);
        if (inst->zbest < best_current_cost) {
            best_current_cost = inst->zbest;
            best_start = i;
            printf("The solution has been updated\n new solution: ");
            for (int j = 0; j < inst->nnodes; j++) {
                printf("%d ", inst->bestsol[j]);
            }
            printf("\n The total cost is now: %f\n", inst->zbest);
        }
    }
    nearest_neighbor(inst, best_start);
    return best_start;
}


// grasp NN given a starting point and the 2 probabities
void nearest_neighbor_grasp(Instance* inst, int start, double p2, double p3) {
    if (start < 0)
        print_error("Invalid choice of the start node");
    int* current = (int*)calloc(inst->nnodes, sizeof(int));
    int n = inst->nnodes;
    int len = 0;
    double tourCost = 0;
    int last = start;
    double min = MAX;
    double min2 = MAX;
    double min3 = MAX;
    int best_pos = -1;
    int best_pos2 = -1;
    int best_pos3 = -1;
    initialize_cost(inst);
    //initialize current
    for (int i = 0; i < inst->nnodes; i++) {
        current[i] = i;
    }
    swap(current, 0, start);
    len++;
    for (int i = 1; i < inst->nnodes; i++) {
        for (int j = len; j < inst->nnodes; j++) {
            if (inst->cost[last * n + current[j]] < min3) {

                if (inst->cost[last * n + current[j]] < min2) {

                    if (inst->cost[last * n + current[j]] < min) {

                        min3 = min2;
                        best_pos3 = best_pos2;
                        min2 = min;
                        best_pos2 = best_pos;
                        min = inst->cost[last * n + current[j]];
                        best_pos = j;

                    }
                    else {

                        min3 = min2;
                        best_pos3 = best_pos2;
                        min2 = inst->cost[last * n + current[j]];
                        best_pos2 = j;

                    }
                }
                else {

                    min3 = inst->cost[last * n + current[j]];
                    best_pos3 = j;

                }
            }
        }
        double c = rand01();
        printf("random : %f\n", c);
        if (c < p3 && min3 < 3000) {
            printf("sono nel caso dove prendo terzo migliore\n");
            tourCost = tourCost + min3;
            swap(current, len, best_pos3);
            min = MAX;
            min2 = MAX;
            min3 = MAX;
            last = current[len];
            len++;
        }
        else if (p3 < c && c < p2 + p3 && min2 < 3000) {
            printf("sono nel caso dove prendo secondo migliore\n");
            tourCost = tourCost + min2;
            swap(current, len, best_pos2);
            min = MAX;
            min2 = MAX;
            min3 = MAX;
            last = current[len];
            len++;
        }
        else if (c > p3 + p2) {
            printf("sono nel caso dove prendo primo migliore\n");
            tourCost = tourCost + min;
            swap(current, len, best_pos);
            min = MAX;
            min2 = MAX;
            min3 = MAX;
            last = current[len];
            len++;
        }
    }
    tourCost += inst->cost[start * n + last];
    inst->bestsol = current;
    inst->zbest = tourCost;
}



