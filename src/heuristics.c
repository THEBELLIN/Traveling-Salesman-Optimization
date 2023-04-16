#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include <assert.h>


//returns the convex hull of points
//SORTS THE POINTS AND CHANGES INDEXES
Point* convex_hull(Point* p, int len, int* hsize) {
    if (len == 0 || p == NULL) {
        *hsize = 0;
        return NULL;
    }

    int i, size = 0, capacity = 4;
    Point* hull = MALLOC(capacity, Point);

    qsort(p, len, sizeof(Point), comparePoints);

    /* lower hull */
    for (i = 0; i < len; ++i) {
        while (size >= 2 && !ccw(&hull[size - 2], &hull[size - 1], &p[i]))
            --size;
        if (size == capacity) {
            capacity *= 2;
            hull = REALLOC(hull, capacity, Point);
        }
        assert(size >= 0 && size < capacity);
        hull[size] = p[i];
        size++;
    }

    /* upper hull */
    int t = size + 1;
    for (i = len - 1; i >= 0; i--) {
        while (size >= t && !ccw(&hull[size - 2], &hull[size - 1], &p[i]))
            --size;
        if (size == capacity) {
            capacity *= 2;
            hull = REALLOC(hull, capacity, Point);
        }
        assert(size >= 0 && size < capacity);
        hull[size] = p[i];
        size++;
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

int* points_to_indexes(Instance* inst, Point* p, int n, int* indexes)
{
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

void extra_mileage(Instance* inst,  em_options* options)
{
    //check for call parameters
    if (options->opt == NORM)
        extra_mileage_det(inst, options->start);
    else if (options->opt == GRASP_2)
        extra_mileage_grasp2(inst, options->start, options->p1);
    else if (options->opt==GRASP_3)
        extra_mileage_grasp3(inst, options->start, options->p1, options->p2);
    else
        print_error("Wrong parameter for extra mileage algorithm");

    //save the cost
    inst->currcost = get_cost(inst, inst->currsol);

    //save if solution found is better
    save_if_best(inst);
}

void extra_mileage_det(Instance* inst, em_start start)
{
    //construct cost matrix if not already done
    if (inst->cost == NULL)
        print_error("Cost matrix not initialized");

    //select starting indexes
    int n = inst->nnodes;
    int* starting_points = NULL;
    int n_starting = 0;

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
        Point* conv = convex_hull(inst->points, inst->nnodes, &n_starting);
        initialize_cost(inst);
        starting_points = MALLOC(n_starting, int);
        points_to_indexes(inst, conv, n_starting, starting_points);
        free(conv);
    }
    else
        print_error("Wrong starting flag in extra mileage");

    //extra mileage loop
    int current_nodes = n_starting;
    for (int i = 0; i < n; i++)
        inst->bestsol[i] = i;
    inst->bestsol[n] = *starting_points;
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
                //printf("\nconsidering edge %d -> %d against point %d", e.from, e.to, inst->bestsol[j]);
                extra_cost = -(COST(e.from, e.to)) + COST(e.from, inst->bestsol[j]) + COST(inst->bestsol[j], e.to);
                if (extra_cost < min_extra_cost)
                {
                    /*printf("removed cost: %f", COST(e.from, e.to));
                    printf("added costs cost: %f, %f", COST(e.from, inst->bestsol[j]), COST(inst->bestsol[j], e.to));
                    printf("\nmin extra cost found: %f, previous was %f. This is point %d between %d and %d", extra_cost, min_extra_cost, inst->bestsol[j], e.from, e.to);*/
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
    free(starting_points);
    //plot_generator(inst, inst->nnodes);
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
        Point* conv = convex_hull(inst->points, inst->nnodes, &n_starting);
        initialize_cost(inst);
        starting_points = MALLOC(n_starting, int);
        points_to_indexes(inst, conv, n_starting, starting_points);
        free(conv);
    }
    else
        print_error("Wrong starting flag in extra mileage");

    //extra mileage loop
    int current_nodes = n_starting;
    for (int i = 0; i < n; i++)
        inst->currsol[i] = i;
    inst->currsol[n] = starting_points[0];
    //set current tour
    for (int i = 0; i < n_starting; i++)
    {
        swap(inst->currsol, i, starting_points[i]);
    }
    swap(inst->currsol, n_starting, n);
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
            e.from = inst->currsol[i];
            e.to = inst->currsol[i + 1];
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
            add_in_position(new_point_index, place, inst->currsol, current_nodes + 1);
        if(second)
            add_in_position(s_new_point_index, s_place, inst->currsol, current_nodes + 1);
        //DEBUG
        /*for (int i = 0; i < current_nodes + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf(" || ");
        for (int i = current_nodes+1; i < n + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf("\n");*/
    }
    //free
    free(starting_points);
}

void extra_mileage_grasp3(Instance* inst, em_start start, double p1, double p2)
{
    //construct cost matrix if not already done
    if (inst->cost == NULL)
        print_error("Cost matrix not initialized");

    //select starting indexes
    int n = inst->nnodes;
    int* starting_points = NULL;
    int n_starting = 0;

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
        Point* conv = convex_hull(inst->points, inst->nnodes, &n_starting);
        initialize_cost(inst);
        starting_points = MALLOC(n_starting, int);
        points_to_indexes(inst, conv, n_starting, starting_points);
        free(conv);
    }
    else
        print_error("Wrong starting flag in extra mileage");

    //extra mileage loop
    int current_nodes = n_starting;
    for (int i = 0; i < n; i++)
        inst->bestsol[i] = i;
    if (starting_points) //!=NULL
        inst->bestsol[n] = *starting_points;
    else
        print_error("%s: Starting points is NULL pointer", __LINE__);
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
    free(starting_points);
}

//greedy NN given a starting point O(n^2)
void nearest_neighbor_det(Instance* inst, int start) 
{
    if (start < 0)
        print_error("Invalid choice of the start node");

    int n = inst->nnodes;
    int len = 0;
    double costo = 0;
    int last = start;
    double min = INF_DOUBLE;
    int best_pos;
    // initializes the array of costs 
    initialize_cost(inst);
    //initialize current with the trivial permutation 0-1-2....-n
    for (int i = 0; i < inst->nnodes; i++) 
    {
        inst->currsol[i] = i;
    }
    inst->currsol[n] = start;
    // moves at the beginnig the first covered node which is the starting node
    swap(inst->currsol, 0, start);
    len++;
    for (int i = 1; i < n; i++) 
    {
        // finds minimum cost node from the last visited node of the tour
        for (int j = len; j < n; j++) 
        {
            double c = COST(last, inst->currsol[j]);
            if (c < min) 
            {
                min = c;
                best_pos = j;
            }
        }
        // update the cost and rearranges the current array
        costo += min;
        swap(inst->currsol, len, best_pos);
        min = INF_DOUBLE;
        last = inst->currsol[len];
        len++;
    }
    costo += COST(start, last);
    inst->currcost = costo;
    save_if_best(inst);
}

// greedy NN with O(n^3) finds the best starting point and its solution
int nearest_neighbor_allstart(Instance* inst) 
{
    double best_current_cost = INF_DOUBLE;
    int best_start;
    for (int i = 0; i < inst->nnodes; i++) 
    {
        nearest_neighbor_det(inst, i);
        if (inst->bestcost < best_current_cost) {
            best_current_cost = inst->bestcost;
            best_start = i;
            /*DEBUG
            printf("The solution has been updated\n new solution: ");
            for (int j = 0; j < inst->nnodes; j++) 
            {
                printf("%d ", inst->bestsol[j]);
            }
            printf("\n The total cost is now: %f\n", inst->bestcost);*/
        }
    }
    return best_start;
}

void nearest_neighbor_grasp(Instance* inst, int start, double p2) 
{
    if (start < 0)
        print_error("Invalid choice of the start node");

    int n = inst->nnodes;
    int len = 0;
    double costo = 0;
    int last = start;
    double min = INF_DOUBLE;
    double min2 = INF_DOUBLE;
    int best_pos = -1;
    int best_pos2 = -1;
    initialize_cost(inst);
    //initialize current
    for (int i = 0; i < inst->nnodes; i++) 
    {
        inst->currsol[i] = i;
    }
    inst->currsol[n] = start;
    swap(inst->currsol, 0, start);
    len++;
    for (int i = 1; i < n; i++) 
    {
        for (int j = len; j < n; j++) 
        {
            double c = COST(last, inst->currsol[j]);
            if (c < min2) 
            {
                if (c < min) 
                {
                    min2 = min;
                    best_pos2 = best_pos;
                    min = c;
                    best_pos = j;
                }
                else 
                {
                    min2 = c;
                    best_pos2 = j;
                }
            }
        }
        double draw = rand01();
        if (draw < p2 ) 
        {
            costo += min2;
            swap(inst->currsol, len, best_pos2);
        }
        else
        {
            costo += min;
            swap(inst->currsol, len, best_pos);
        }
        min = INF_DOUBLE;
        min2 = INF_DOUBLE;
        last = inst->currsol[len];
        len++;
    }

    costo += COST(start, last);
    inst->currcost = costo;
    save_if_best(inst);
}

// grasp NN given a starting point and the 2 probabities
void nearest_neighbor_grasp2(Instance* inst, int start, double p2, double p3) 
{
    if (start < 0)
    {
        print_error("%s, Invalid choice of the start node", __LINE__);
    }

    int n = inst->nnodes;
    int len = 0;
    double costo = 0;
    int last = start;
    double min = INF_DOUBLE;
    double min2 = INF_DOUBLE;
    double min3 = INF_DOUBLE;
    int best_pos = -1;
    int best_pos2 = -1;
    int best_pos3 = -1;
    initialize_cost(inst);

    //initialize current
    for (int i = 0; i < inst->nnodes; i++) {
        inst->currsol[i] = i;
    }
    inst->currsol[n] = start;
    swap(inst->currsol, 0, start);
    len++;

    for (int i = 1; i < inst->nnodes; i++) 
    {
        for (int j = len; j < inst->nnodes; j++) 
        {
            double c = COST(last, inst->currsol[j]);
            if (c < min3) 
            {
                if (c < min2) 
                {
                    if (c < min) 
                    {
                        min3 = min2;
                        best_pos3 = best_pos2;
                        min2 = min;
                        best_pos2 = best_pos;
                        min = c;
                        best_pos = j;
                    }
                    else 
                    {
                        min3 = min2;
                        best_pos3 = best_pos2;
                        min2 = c;
                        best_pos2 = j;
                    }
                }
                else 
                {
                    min3 = c;
                    best_pos3 = j;
                }
            }
        }
        double draw = rand01();
        if (draw < p3 ) 
        {
            costo = costo + min3;
            swap(inst->currsol, len, best_pos3);
        }
        else if (p3 < draw && draw < p2 + p3 ) 
        {
            costo = costo + min2;
            swap(inst->currsol, len, best_pos2);
        }
        else if (draw > p3 + p2) 
        {
            costo = costo + min;
            swap(inst->bestsol, len, best_pos);
        }
        min = INF_DOUBLE;
        min2 = INF_DOUBLE;
        min3 = INF_DOUBLE;
        last = inst->bestsol[len];
        len++;
    }
    costo += COST(start, last);
    inst->currcost = costo;

    //save if newbest sol found
    save_if_best(inst);
}

// grasp where every 20 iteration it makes a random pick
void nearest_neighbor_grasp_random(Instance* inst, int start, double p2) 
{
    if (start < 0)
        print_error("Invalid choice of the start node");

    int n = inst->nnodes;
    int len = 0;
    double costo = 0;
    int last = start;
    double min = INF_DOUBLE;
    double min2 = INF_DOUBLE;
    int best_pos = -1;
    int best_pos2 = -1;
    initialize_cost(inst);

    //initialize current
    for (int i = 0; i < inst->nnodes; i++) {
        inst->currsol[i] = i;
    }
    inst->currsol[n] = start;
    swap(inst->currsol, 0, start);
    len++;

    for (int i = 1; i < inst->nnodes; i++) 
    {
        //every 20 iterations do a random add
        if (i % 20 == 0) 
        {
            int d = rand_int(len, n-1);
            costo += COST(last, inst->currsol[d]);
            swap(inst->currsol, len, d);
            last = inst->currsol[len];
            len++;
        }
        else
        {
            //compute min costs
            for (int j = len; j < inst->nnodes; j++) 
            {
                double c = COST(last, inst->currsol[j]);
                if (c < min2) 
                {
                    if (c < min) 
                    {
                        min2 = min;
                        best_pos2 = best_pos;
                        min = c;
                        best_pos = j;
                    }
                    else 
                    {
                        min2 = c;
                        best_pos2 = j;
                    }
                }

            }
            double draw = rand01();
            if (draw < p2) 
            {
                costo += min2;
                swap(inst->currsol, len, best_pos2);
            }
            else
            {
                costo += min;
                swap(inst->currsol, len, best_pos);
            }
            min = INF_DOUBLE;
            min2 = INF_DOUBLE;
            last = inst->currsol[len];
            len++;
        }
    }
    costo += COST(start, last);
    inst->currcost = costo;
    save_if_best(inst);
}


void next_bestsol(Instance* inst, int it) {
    //assume we have already found a solution
    //this funxtions returns the next best solution in the neighborhood that is not currsol and satisfies tabu list

    //===============2-OPT CHECK=============
    int best_arc = -1;
    int best_arc2 = -1;
    double best_delta = INF_DOUBLE;
    int n = inst->nnodes;

    //check for tabu tenure validity
    if (inst->tabu_tenure < 0)
        print_error("%s, Impossible value of tabu tenure", __LINE__);

    //check for the best 2-opt move possible
    for (int i = 0; i < n - 2; i++) 
    {
        //check for tabu nodes
        if (it - inst->tabu[inst->bestsol[i]] < inst->tabu_tenure || it - inst->tabu[inst->bestsol[i + 1]] < inst->tabu_tenure)
            continue;
        for (int j = i + 2; j < n - 1; j++) 
        {
            //check for tabu nodes
            if (it - inst->tabu[inst->bestsol[j]] < inst->tabu_tenure || it - inst->tabu[inst->bestsol[j + 1]] < inst->tabu_tenure)
                continue;
            // calculate the delta cost of swapping edges bestsol[(i,i+1)] and (j,j+1) for (i,j) and (i+1,j+1)
            double delta = COST(inst->bestsol[i], inst->bestsol[j]) + COST(inst->bestsol[i + 1], inst->bestsol[j + 1]) - COST(inst->bestsol[i], inst->bestsol[i + 1]) - COST(inst->bestsol[j], inst->bestsol[j + 1]);
            if (delta < best_delta) 
            {
                best_delta = delta;
                best_arc = i;
                best_arc2 = j;
            }
        }
    }
    if (best_delta == INF_DOUBLE)
        print_error("%s, Error in finding another solution", __LINE__);

    //put all 4 nodes swapped in tabu list
    inst->tabu[inst->currsol[best_arc]] = it;
    inst->tabu[inst->currsol[best_arc + 1]] = it;
    inst->tabu[inst->currsol[best_arc2]] = it;
    inst->tabu[inst->currsol[best_arc2 + 1]] = it;

    invert_nodes(inst->currsol, best_arc + 1, best_arc2); // reverse the sub-tour between best_arc and best_arc2
    inst->currcost += best_delta;
    if (inst->currcost < inst->bestcost) //if new solution is best found, save it
    {
        copy_array(inst->currsol, inst->bestsol, inst->nnodes + 1);
        inst->bestcost = inst->currcost;
        if (inst->verbose > 0)
            printf("New best solution found of cost: %f", inst->bestcost);
    }
}

void solve(Instance* inst, solve_options* options)
{
    inst->time_limit = options->timelimit;
    if (options->alg == EM)
    {
        extra_mileage(inst, options->em_opts);
    }
    else if (options->alg == NN)
    {
        nearest_neighbor(inst, 0); //TODO fix this
    }
    else if (options->alg == GEN)
    {
        genetic(inst); 
        return;
    }
    else
        print_error("%s, Error in setting algorithm options for solving", __LINE__);

    printf("starting to iterate");
    //iterate untill timelimit is exceeded
    int it = 0;
    int time_elapsed = time(NULL) - inst->tstart;
    while ( time_elapsed < options->timelimit)
    {
        next_bestsol(inst, it);
        it++;
        time_elapsed = time(NULL) - inst->tstart;
    }
    //plot best instance found
    //plot_generator(inst, inst->nnodes);
}


//plot cost of incumbent during solving (or cost in console but better to plot)
//tabulate formatted output

void VNS(Instance* inst) 
{
    // initialize solution
    nearest_neighbor_grasp_random(inst, 0, 0.5);
    int n = inst->nnodes;
    int i = 0;
    //copy bestsol to currsol since two_opt works on currsol
    copy_array(inst->bestsol, inst->currsol, inst->nnodes + 1);
    inst->currcost = inst->bestcost;
    two_opt(inst);

    while (time(NULL) - inst->tstart < inst->time_limit) 
    {
        kick(inst, inst->currsol);
        kick(inst, inst->currsol);
        two_opt(inst);
        save_if_best(inst);
    }
}

void two_opt(Instance* inst) 
{
    int n = inst->nnodes;
    double best_delta = - 1;
    while (best_delta < 0) 
    {
        best_delta = 0;
        int best_i = -1;
        int best_j = -1;
        for (int i = 0; i < n - 2; i++) 
        {
            for (int j = i + 2; j < n; j++) 
            {
                double dist_before = COST(inst->currsol[i], inst->currsol[i + 1]) + COST(inst->currsol[j], inst->currsol[j + 1]);
                double dist_after = COST(inst->currsol[i], inst->currsol[j]) + COST(inst->currsol[i + 1], inst->currsol[j + 1]);
                double delta = dist_after - dist_before;
                if (delta < best_delta) 
                {
                    best_i = i;
                    best_j = j;
                    best_delta = delta;
                }
            }
        }
        if (best_delta < 0) 
        {
            invert_nodes(inst->currsol, best_i + 1, best_j);
            inst->currcost += best_delta;
        }
        save_if_best(inst);
    }
}

void kick(Instance* inst, int* sol) {
    int n = inst->nnodes;

    int* randomEdges = (int*)calloc(3, sizeof(int));

    randomEdges[0] = rand() % (n - 1);
    // choses at random 3 nodes which represent the arcs to eliminate checking if they are equal and sorting them
    for (int i = 1; i < 3; i++) {
        double c = rand() % (n - 1);
        if (checkNonEqual(randomEdges, c, i)) {
            randomEdges[i] = c;
            for (int j = i; j > 0; j--) {
                if (randomEdges[j] < randomEdges[j - 1])
                    swap(randomEdges, j, j - 1);
            }
        }
        else
            i--;

    }

    int* solutionRearrenged = (int*)calloc(n + 1, sizeof(int));
    int i = 0;
    int len = 0;
    // copio fino al primo indice
    while (sol[i] != sol[randomEdges[0] + 1]) {
        solutionRearrenged[len] = sol[i];
        i++;
        len++;
    }
    // dopo sol[indice[0]] volgio sol[indice[1]+1]
    i = randomEdges[1] + 1;

    while (sol[i] != sol[randomEdges[2] + 1]) {
        solutionRearrenged[len] = sol[i];
        i++;
        len++;
    }
    // da sol[indice 2] voglio sol[[indice 0]+1]
    i = randomEdges[0] + 1;
    while (sol[i] != sol[randomEdges[1] + 1]) {
        solutionRearrenged[len] = sol[i];
        i++;
        len++;
    }
    //solutionRearrenged[len] = sol[indice[2] + 1];
   // len++;
    i = randomEdges[2] + 1;
    while (sol[i] != sol[0]) {
        solutionRearrenged[len] = sol[i];
        i++;
        len++;
    }
    solutionRearrenged[n] = sol[0];
    //memcpy(sol, solutionRearrenged, sizeof(solutionRearrenged));
    for (i = 0; i < n + 1; i++) {
        sol[i] = solutionRearrenged[i];
    }
    free(randomEdges);
    free(solutionRearrenged);
}


int checkNonEqual(int* array, int tocheck, int length) {
    for (int i = 0; i < length; i++) {
        if (array[i] == tocheck)
            return 0;
    }
    return 1;
}


