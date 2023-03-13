#include "TSP.h"
#include "utility.h"

int comparePoints(const Point* p1, const Point* p2) 
{
    //first compare x
    //note that with double coordinates it will always solve at x
    if (p1->x < p2->x)
        return -1;
    if (p2->x < p1->x)
        return 1;
    //then y
    if (p1->y < p2->y)
        return -1;
    if (p2->y < p1->y)
        return 1;
    return 0;
}

//check if the triangle made of a,b,c is counterclockwise oriented and not flat
//"p2 is under p3"
// figure and explaination at https://www.lri.fr/~marche/MPRI-2-36-1/2014/hull.pdf
bool ccw(const Point* a, const Point* b, const Point* c) 
{
    return (b->x - a->x) * (c->y - a->y) > (b->y - a->y) * (c->x - a->x);
}

void swap(int* array, int ind1, int ind2) {
    int temp;
    temp = array[ind1];
    array[ind1] = array[ind2];
    array[ind2] = temp;
}

// initialize the array of costs
void initialize_cost(Instance* inst) {
    int n = inst->nnodes;
    int size = n * n;
    inst->cost = (double*)calloc(size, sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inst->cost[i * n + j] = distance(&(inst->points[i]), &(inst->points[j]));
        }
    }
}

//comparison between double
bool is_equal_double(double a, double b)
{
    return fabs(a - b) < EPSILON;
}

//check if two points are equal
bool is_equal_points(Point* p1, Point* p2)
{
    return(is_equal_double(p1->x, p2->x) && is_equal_double(p1->y, p2->y));
}

double rand01()
{
    return (double)rand() / RAND_MAX;
}