#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include <assert.h>

//returns the convex hull of points in instance
Point* convex_hull(Instance* inst, int* hsize) //mostly found online
{
	if (inst->nnodes < 1)
		return NULL;

    Point* hull;
    Point* points_copy;
	int capacity = 4;
    int size = 0;
    int n = inst->nnodes;
	hull = MALLOC(capacity, Point);
    points_copy = MALLOC(n, Point);
    for (int i = 0; i < n; i++)
        points_copy[i] = inst->points[i];

	qsort(points_copy, n, sizeof(Point), comparePoints);

    /* lower hull */
    for (int i = 0; i < n; ++i) 
    {
        while (size >= 2 && !ccw(&hull[size - 2], &hull[size - 1], &points_copy[i]))
            --size;
        if (size == capacity) 
        {
            capacity *= 2;
            hull = REALLOC(hull, capacity, Point);
        }
        assert(size >= 0 && size < capacity);
        hull[size++] = points_copy[i];
    }

    /* upper hull */
    int t = size + 1;
    for (int i = n - 1; i >= 0; i--) {
        while (size >= t && !ccw(&hull[size - 2], &hull[size - 1], &points_copy[i]))
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
    free(points_copy);
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

int* points_to_indexes(Instance* inst, Point* p, int n)
{
    int* indexes = MALLOC(n, int);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < inst->nnodes; j++)
        {
            if (is_equal_points(&p[i], &inst->points[j]))
            {
                indexes[i] = j;
                break;
            }
        }
    }
    return indexes;
}

void extra_mileage(Instance* inst, em_start start, em_options* options)
{
    if (options->opt == NORM)
        extra_mileage_det(inst, start);
    else if (options->opt == GRASP_2)
        extra_mileage_grasp2(inst, start, options->p1);
    else if (options->opt==GRASP_3)
        extra_mileage_grasp3(inst, start, options->p1, options->p2);
    else
        print_error("Wrong parameter for extra mileage algorithm");
}

void extra_mileage_det(Instance* inst, em_start start)
{
    //construct cost matrix if not already done
    if (inst->cost == NULL)
        print_error("Cost matrix not initialized");

    //select starting indexes
    int n = inst->nnodes;
    int* starting_points = NULL;
    int n_starting;

    if (start == RAND)
    {
        starting_points = MALLOC(2, int);
        n_starting = 2;
        starting_points[0] = rand() % n;
        starting_points[1] = starting_points[0];
        while(starting_points[1]==starting_points[0])
            starting_points[1] = rand() % n;
    }
    else if (start == MAX_DIST)
    {
        starting_points = MALLOC(2, int);
        n_starting = 2;
        int maxcost = 0;
        starting_points[0] = 0;
        starting_points[1] = 0;
        for (int i = 1; i < (n * n); i++)
        {
            if (inst->cost[i] > maxcost)
            {
                maxcost = inst->cost[i];
                starting_points[0] = i / n;
                starting_points[1] = i % n;
            }
        }
    }
    else if (start == CONV_HULL)
    {
        Point* conv = convex_hull(inst, &n_starting);
        starting_points = points_to_indexes(inst, conv, n_starting);
    }
    else
        print_error("Wrong starting flag in extra mileage");

    //extra mileage loop
    int current_nodes = n_starting;
    for (int i = 0; i < n; i++)
        inst->bestsol[i] = i;
    inst->bestsol[n] = starting_points[0];
    //set current tour
    for (int i = 0; i < n_starting; i++)
    {
        swap(inst->bestsol, i, starting_points[i]);
    }
    swap(inst->bestsol, n_starting, n);
    //untill all nodes are added
    while (current_nodes < n)
    {
        double extra_cost ;
        double min_extra_cost = INF_DOUBLE;
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
                extra_cost = - COST(e.from, e.to) + COST(e.from, inst->bestsol[j]) + COST(inst->bestsol[j], e.to);
                if (extra_cost < min_extra_cost)
                {
                    min_extra_cost = extra_cost;
                    new_point_index = j;
                    place = i + 1;
                }
            }
        }
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
}

void extra_mileage_grasp2(Instance* inst, em_start start, double p)
{
    //construct cost matrix if not already done
    if (inst->cost == NULL)
        print_error("Cost matrix not initialized");

    //select starting indexes
    int n = inst->nnodes;
    int* starting_points = NULL;
    int n_starting;

    if (start == RAND)
    {
        starting_points = MALLOC(2, int);
        n_starting = 2;
        starting_points[0] = rand() % n;
        starting_points[1] = starting_points[0];
        while (starting_points[1] == starting_points[0])
            starting_points[1] = rand() % n;
    }
    else if (start == MAX_DIST)
    {
        starting_points = MALLOC(2, int);
        n_starting = 2;
        int maxcost = 0;
        starting_points[0] = 0;
        starting_points[1] = 0;
        for (int i = 1; i < (n * n); i++)
        {
            if (inst->cost[i] > maxcost)
            {
                maxcost = inst->cost[i];
                starting_points[0] = i / n;
                starting_points[1] = i % n;
            }
        }
    }
    else if (start == CONV_HULL)
    {
        Point* conv = convex_hull(inst, &n_starting);
        starting_points = points_to_indexes(inst, conv, n_starting);
    }
    else
        print_error("Wrong starting flag in extra mileage");

    //extra mileage loop
    int current_nodes = n_starting;
    for (int i = 0; i < n; i++)
        inst->bestsol[i] = i;
    inst->bestsol[n] = starting_points[0];
    //set current tour
    for (int i = 0; i < n_starting; i++)
    {
        swap(inst->bestsol, i, starting_points[i]);
    }
    swap(inst->bestsol, n_starting, n);
    //untill all nodes are added
    while (current_nodes < n)
    {
        double draw = rand01();    // [0,p] -> first, [p, 1] -> second
        bool second = false;
        if (draw > p)
            second = true;
        double extra_cost;
        double min_extra_cost = INF_DOUBLE, s_min_extra_cost = INF_DOUBLE;
        int place = -1, s_place = -1;
        int new_point_index = -1, s_new_point_index = -1;
        //for every edge in current tour
        for (int i = 0; i < current_nodes; i++) //considering edge inst->bestsol[i]->inst->bestsol[i+1]
        {
            edge e;
            e.from = inst->bestsol[i];
            e.to = inst->bestsol[i + 1];
            //for all points not already considered
            for (int j = current_nodes + 1; j < n + 1; j++) //considering point inst->bestsol[j]
            {
                extra_cost = -COST(e.from, e.to) + COST(e.from, inst->bestsol[j]) + COST(inst->bestsol[j], e.to);
                if (extra_cost < min_extra_cost)
                {
                    s_min_extra_cost = min_extra_cost;
                    min_extra_cost = extra_cost;
                    s_new_point_index = new_point_index;
                    new_point_index = j;
                    s_place = place;
                    place = i + 1;
                }
                else if (extra_cost < s_min_extra_cost)
                {
                    s_min_extra_cost = extra_cost;
                    s_new_point_index = j;
                    s_place = i + 1;
                }
            }
        }
        current_nodes++;
        if(!second)
            add_in_position(new_point_index, place, inst->bestsol, current_nodes + 1);
        if(second)
            add_in_position(s_new_point_index, s_place, inst->bestsol, current_nodes + 1);
        //DEBUG
        /*for (int i = 0; i < current_nodes + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf(" || ");
        for (int i = current_nodes+1; i < n + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf("\n");*/
    }
}

void extra_mileage_grasp3(Instance* inst, em_start start, double p1, double p2)
{
    //construct cost matrix if not already done
    if (inst->cost == NULL)
        print_error("Cost matrix not initialized");

    //select starting indexes
    int n = inst->nnodes;
    int* starting_points = NULL;
    int n_starting;

    if (start == RAND)
    {
        starting_points = MALLOC(2, int);
        n_starting = 2;
        starting_points[0] = rand() % n;
        starting_points[1] = starting_points[0];
        while (starting_points[1] == starting_points[0])
            starting_points[1] = rand() % n;
    }
    else if (start == MAX_DIST)
    {
        starting_points = MALLOC(2, int);
        n_starting = 2;
        int maxcost = 0;
        starting_points[0] = 0;
        starting_points[1] = 0;
        for (int i = 1; i < (n * n); i++)
        {
            if (inst->cost[i] > maxcost)
            {
                maxcost = inst->cost[i];
                starting_points[0] = i / n;
                starting_points[1] = i % n;
            }
        }
    }
    else if (start == CONV_HULL)
    {
        Point* conv = convex_hull(inst, &n_starting);
        starting_points = points_to_indexes(inst, conv, n_starting);
    }
    else
        print_error("Wrong starting flag in extra mileage");

    //extra mileage loop
    int current_nodes = n_starting;
    for (int i = 0; i < n; i++)
        inst->bestsol[i] = i;
    inst->bestsol[n] = starting_points[0];
    //set current tour
    for (int i = 0; i < n_starting; i++)
    {
        swap(inst->bestsol, i, starting_points[i]);
    }
    swap(inst->bestsol, n_starting, n);
    //untill all nodes are added
    while (current_nodes < n)
    {
        double draw = rand01();    // [0,p1] -> first, [p1, p2] -> second, [p2, 1] -> third
        bool second = false;
        bool third = false;
        if (draw > p2)
            third = true;
        else if (draw > p1)
            second = true;
        double extra_cost;
        double min_extra_cost = INF_DOUBLE, s_min_extra_cost = INF_DOUBLE, t_min_extra_cost = INF_DOUBLE;
        int place = -1, s_place = -1, t_place = -1;
        int new_point_index = -1, s_new_point_index = -1, t_new_point_index = -1;
        //for every edge in current tour
        for (int i = 0; i < current_nodes; i++) //considering edge inst->bestsol[i]->inst->bestsol[i+1]
        {
            edge e;
            e.from = inst->bestsol[i];
            e.to = inst->bestsol[i + 1];
            //for all points not already considered
            for (int j = current_nodes + 1; j < n + 1; j++) //considering point inst->bestsol[j]
            {
                extra_cost = -COST(e.from, e.to) + COST(e.from, inst->bestsol[j]) + COST(inst->bestsol[j], e.to);
                if (extra_cost < min_extra_cost)
                {
                    t_min_extra_cost = s_min_extra_cost;
                    s_min_extra_cost = min_extra_cost;
                    min_extra_cost = extra_cost;
                    t_new_point_index = s_new_point_index;
                    s_new_point_index = new_point_index;
                    new_point_index = j;
                    t_place = s_place;
                    s_place = place;
                    place = i + 1;
                }
                else if (extra_cost < s_min_extra_cost)
                {
                    t_min_extra_cost = s_min_extra_cost;
                    s_min_extra_cost = extra_cost;
                    t_new_point_index = s_new_point_index;
                    s_new_point_index = j;
                    t_place = s_place;
                    s_place = i + 1;
                }
                else if (extra_cost < t_min_extra_cost)
                {
                    t_min_extra_cost = extra_cost;
                    t_new_point_index = j;
                    t_place = i + 1;
                }
            }
        }
        current_nodes++;
        if (second)
            add_in_position(s_new_point_index, s_place, inst->bestsol, current_nodes + 1);
        else if(third)
            add_in_position(t_new_point_index, t_place, inst->bestsol, current_nodes + 1);
        else
            add_in_position(new_point_index, place, inst->bestsol, current_nodes + 1);
        //DEBUG
        /*for (int i = 0; i < current_nodes + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf(" || ");
        for (int i = current_nodes+1; i < n + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf("\n");*/
    }
}