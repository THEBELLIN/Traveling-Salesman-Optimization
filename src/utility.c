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