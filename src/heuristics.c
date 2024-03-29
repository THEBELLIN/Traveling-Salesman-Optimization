
#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include <assert.h>
// grasp iter nn tested it works also retested all nearest neighbor
// extra mileage corrected it didn't add the last node + added iterative grasp for extra mileage 
// simulated annealing tested and added in heuristics.c just for now
// two opt added as a method with its best solution initialization(nn grasp2)
// tabu seems to work but the cost is wrong so i added get cost in order to have the correct cost else it kept giving a final result different from the real one saved
// vns works





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
    for (int i = len; i > pos; i--)
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


void extra_mileage_det(Instance* inst)
{
    //construct cost matrix if not already done
    if (inst->cost == NULL)
        print_error("Cost matrix not initialized", __LINE__);

    //select starting indexes
    int n = inst->nnodes;
    int* starting_points = NULL;
    int n_starting = 0;
    int start = inst->solver.start;

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
        print_error("Wrong starting flag in extra mileage", __LINE__);

    //extra mileage loop
    int current_nodes = n_starting;
    int* utilized = CALLOC(n, int);
    //set the initial points in current solution
    for (int i = 0; i < n_starting; i++)
    {
        inst->currsol[i] = starting_points[i];
        utilized[starting_points[i]] = 1;
    }
    inst->currsol[n_starting] = starting_points[0];
    //set remaining points
    int k = 0;
    for (int i = n_starting + 1; i < n + 1; i++)
    {
        while (utilized[k] == 1)
            k++;
        inst->currsol[i] = k;
        utilized[k] = 1;
    }
    free(utilized);

    //untill all nodes are added
    while (current_nodes < n)
    {
        double extra_cost;
        double min_extra_cost = INF_DOUBLE;
        int place;
        int new_point_index;
        //for every edge in current tour
        for (int i = 0; i < current_nodes; i++) //considering edge inst->bestsol[i]->inst->bestsol[i+1]
        {
            myedge e;
            e.from = inst->currsol[i];
            e.to = inst->currsol[i + 1];
            //for all points not already considered
            for (int j = current_nodes + 1; j < n + 1; j++) //considering point inst->currsol[j]
            {
                //printf("\nconsidering edge %d -> %d against point %d", e.from, e.to, inst->currsol[j]);
                extra_cost = -(COST(e.from, e.to)) + COST(e.from, inst->currsol[j]) + COST(inst->currsol[j], e.to);
                if (extra_cost < min_extra_cost)
                {
                    /*printf("removed cost: %f", COST(e.from, e.to));
                    printf("added costs cost: %f, %f", COST(e.from, inst->currsol[j]), COST(inst->currsol[j], e.to));
                    printf("\nmin extra cost found: %f, previous was %f. This is point %d between %d and %d", extra_cost, min_extra_cost, inst->bestsol[j], e.from, e.to);*/
                    min_extra_cost = extra_cost;
                    new_point_index = j;
                    place = i + 1;
                }
            }
        }
        current_nodes++;
        add_in_position(new_point_index, place, inst->currsol, current_nodes + 1);
        //DEBUG
        /*for (int i = 0; i < current_nodes + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf(" || ");
        for (int i = current_nodes+1; i < n + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf("\n");*/
    }
    free(starting_points);

    //save the cost
    inst->currcost = get_cost(inst, inst->currsol);

    //save if solution found is better
    save_if_best(inst);
    //plot_generator(inst, inst->nnodes);
}





void extra_mileage(Instance* inst)
{
    //check for call parameters
    if (inst->solver.id == MILEAGE)
        extra_mileage_det(inst);
    else if (inst->solver.id == EM_GRASP2 || inst->solver.id == EM_GRASP3)
        grasp_extra_iter(inst);
    else
        print_error("Wrong parameter for extra mileage algorithm", __LINE__);


}



void extra_mileage_grasp2(Instance* inst)
{
    //construct cost matrix if not already done
    if (inst->cost == NULL)
        print_error("Cost matrix not initialized", __LINE__);

    //select starting indexes
    int n = inst->nnodes;
    int* starting_points = NULL;
    int n_starting;
    int start = inst->solver.start;
    int p = inst->solver.p1;

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
        print_error("Wrong starting flag in extra mileage", __LINE__);

    //extra mileage loop
    int current_nodes = n_starting;
    int* utilized = CALLOC(n, int);
    //set the initial points in current solution
    for (int i = 0; i < n_starting; i++)
    {
        inst->currsol[i] = starting_points[i];
        utilized[starting_points[i]] = 1;
    }
    inst->currsol[n_starting] = starting_points[0];
    //set remaining points
    int k = 0;
    for (int i = n_starting + 1; i < n + 1; i++)
    {
        while (utilized[k] == 1)
            k++;
        inst->currsol[i] = k;
        utilized[k] = 1;
    }
    free(utilized);
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
        for (int i = 0; i < current_nodes; i++) //considering edge inst->currsol[i]->inst->currsol[i+1]
        {
            myedge e;
            e.from = inst->currsol[i];
            e.to = inst->currsol[i + 1];
            //for all points not already considered
            for (int j = current_nodes + 1; j < n + 1; j++) //considering point inst->currsol[j]
            {
                extra_cost = -COST(e.from, e.to) + COST(e.from, inst->currsol[j]) + COST(inst->currsol[j], e.to);
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
        if (!second)
            add_in_position(new_point_index, place, inst->currsol, current_nodes + 1);
        if (second)
            add_in_position(s_new_point_index, s_place, inst->currsol, current_nodes + 1);
        //DEBUG
        /*for (int i = 0; i < current_nodes + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf(" || ");
        for (int i = current_nodes+1; i < n + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf("\n");*/
    }

    //save the cost
    inst->currcost = get_cost(inst, inst->currsol);

    //save if solution found is better
    save_if_best(inst);
    //free
    free(starting_points);
}

void extra_mileage_grasp3(Instance* inst)
{
    //construct cost matrix if not already done
    if (inst->cost == NULL)
        print_error("Cost matrix not initialized", __LINE__);

    //select starting indexes
    int n = inst->nnodes;
    int* starting_points = NULL;
    int n_starting = 0;
    int start = inst->solver.start;
    int p1 = inst->solver.p1;
    int p2 = inst->solver.p2;

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
        print_error("Wrong starting flag in extra mileage", __LINE__);

    //extra mileage loop
    int current_nodes = n_starting;
    int* utilized = CALLOC(n, int);
    //set the initial points in current solution
    for (int i = 0; i < n_starting; i++)
    {
        inst->currsol[i] = starting_points[i];
        utilized[starting_points[i]] = 1;
    }
    inst->currsol[n_starting] = starting_points[0];
    //set remaining points
    int k = 0;
    for (int i = n_starting + 1; i < n + 1; i++)
    {
        while (utilized[k] == 1)
            k++;
        inst->currsol[i] = k;
        utilized[k] = 1;
    }
    free(utilized);

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
        for (int i = 0; i < current_nodes; i++) //considering edge inst->currsol[i]->inst->currsol[i+1]
        {
            myedge e;
            e.from = inst->currsol[i];
            e.to = inst->currsol[i + 1];
            //for all points not already considered
            for (int j = current_nodes + 1; j < n + 1; j++) //considering point inst->currsol[j]
            {
                extra_cost = -COST(e.from, e.to) + COST(e.from, inst->currsol[j]) + COST(inst->currsol[j], e.to);
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
            add_in_position(s_new_point_index, s_place, inst->currsol, current_nodes + 1);
        else if (third)
            add_in_position(t_new_point_index, t_place, inst->currsol, current_nodes + 1);
        else
            add_in_position(new_point_index, place, inst->currsol, current_nodes + 1);
        //DEBUG
        /*for (int i = 0; i < current_nodes + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf(" || ");
        for (int i = current_nodes+1; i < n + 1; i++)
            printf("%d ", inst->bestsol[i]);
        printf("\n");*/
    }
    //save the cost
    inst->currcost = get_cost(inst, inst->currsol);

    //save if solution found is better
    save_if_best(inst);
    free(starting_points);
}


// iterative method 
void grasp_extra_iter(Instance* inst) {

    inst->tstart = time(NULL);
    int time_elapsed = time(NULL) - inst->tstart;
    if (inst->verbose > 1 && inst->solver.id == EM_GRASP2)
        printf("\nStarting iterations of extra milegae grasp2");
    else if (inst->verbose > 1 && inst->solver.id == EM_GRASP3)
        printf("\nStarting iterations of extra milegae grasp3");
    while (time_elapsed < inst->time_limit)
    {
        if (inst->solver.id == EM_GRASP2) {
            extra_mileage_grasp2(inst);
        }
        else if (inst->solver.id == EM_GRASP3) {
            extra_mileage_grasp3(inst);
        }

        time_elapsed = time(NULL) - inst->tstart;
    }
}


void nearest_neighbor(Instance* inst)
{
    //check for call parameters
    if (inst->solver.id == NN)
        nearest_neighbor_det(inst, rand_int(0, inst->nnodes));
    else if (inst->solver.id == GRASP2_NN || inst->solver.id == GRASP3_NN)
        grasp_iter(inst);
    else if (inst->solver.id == ALL_START_NN)
        nearest_neighbor_allstart(inst);
    else
        print_error("Wrong parameter for extra mileage algorithm", __LINE__);

    //save the cost
    inst->currcost = get_cost(inst, inst->currsol);

    //save if solution found is better
    save_if_best(inst);
}

//greedy NN given a starting point O(n^2)
void nearest_neighbor_det(Instance* inst, int start)
{
    if (start < 0)
        print_error("Invalid choice of the start node", __LINE__);

    int n = inst->nnodes;
    int len = 0;
    double costo = 0;
    int last = start;
    double min = INF_DOUBLE;
    int best_pos;
    // initializes the array of costs 
   //initialize_cost(inst);
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
    printf("\nStarting iterations for best nearest neighbor all start");
    double best_current_cost = INF_DOUBLE;
    int best_start = 0;
    for (int i = 0; i < inst->nnodes; i++)
    {
        nearest_neighbor_det(inst, i);
        save_if_best(inst);
    }
    return best_start;
}

void grasp_iter(Instance* inst)
{

    inst->tstart = time(NULL);
    int start = rand_int(0, inst->nnodes);
    int time_elapsed = time(NULL) - inst->tstart;
    if (inst->verbose > 1 && inst->solver.id == GRASP2_NN)
        printf("\nStarting iterations of nearest neighbor grasp2");
    else if (inst->verbose > 1 && inst->solver.id == GRASP3_NN)
        printf("\nStarting iterations of nearest neighbor grasp3");
    while (time_elapsed < inst->time_limit)
    {
        if (inst->solver.id == GRASP2_NN) {

            nearest_neighbor_grasp2(inst, start);
        }
        else if (inst->solver.id == GRASP3_NN) {

            nearest_neighbor_grasp3(inst, start);
        }

        // save_if_best(inst);
        time_elapsed = time(NULL) - inst->tstart;
    }


}

void nearest_neighbor_grasp2(Instance* inst, int start) {

    int n = inst->nnodes;
    int len = 0;
    double p2 = 0.1;
    double costo = 0;
    int last = start;
    double min = INF_DOUBLE;
    double min2 = INF_DOUBLE;
    int best_pos = -1;
    int best_pos2 = -1;
    //initialize_cost(inst);
    //initialize inst->currsol
    for (int i = 0; i < inst->nnodes; i++) {
        inst->currsol[i] = i;
    }
    inst->currsol[n] = start;
    swap(inst->currsol, 0, start);
    len++;
    for (int i = 1; i < inst->nnodes; i++) {
        for (int j = len; j < inst->nnodes; j++) {

            if (inst->cost[last * n + inst->currsol[j]] < min2) {

                if (inst->cost[last * n + inst->currsol[j]] < min) {
                    min2 = min;
                    best_pos2 = best_pos;
                    min = inst->cost[last * n + inst->currsol[j]];
                    best_pos = j;
                }
                else {
                    min2 = inst->cost[last * n + inst->currsol[j]];
                    best_pos2 = j;
                }
            }

        }
        double c = rand01();
        // printf("random : %f\n", c);
        if (c < p2 && min2 < INF_DOUBLE) {
            // printf("sono nel caso dove prendo secondo migliore\n");
            costo = costo + min2;
            swap(inst->currsol, len, best_pos2);

            last = inst->currsol[len];
            len++;
        }
        else
        {
            // printf("sono nel caso dove prendo primo migliore\n");
            costo = costo + min;
            swap(inst->currsol, len, best_pos);

            last = inst->currsol[len];
            len++;
        }
        min = INF_DOUBLE;
        min2 = INF_DOUBLE;
    }

    costo += inst->cost[start * n + last];
    inst->currcost = costo;
    save_if_best(inst);
}

void nearest_neighbor_grasp3(Instance* inst, int start) {
    if (start < 0)
        print_error("Invalid choice of the start node", __LINE__);
    int n = inst->nnodes;
    int len = 0;
    double costo = 0;
    int last = start;
    double p2 = inst->solver.p1;
    double p3 = inst->solver.p2;
    double min = INF_DOUBLE;
    double min2 = INF_DOUBLE;
    double min3 = INF_DOUBLE;
    int best_pos = -1;
    int best_pos2 = -1;
    int best_pos3 = -1;
    //initialize_cost(inst);
    //initialize inst->current 
    for (int i = 0; i < inst->nnodes; i++) {
        inst->currsol[i] = i;
    }
    inst->currsol[n] = start;
    swap(inst->currsol, 0, start);
    len++;
    for (int i = 1; i < inst->nnodes; i++) {
        for (int j = len; j < inst->nnodes; j++) {
            if (inst->cost[last * n + inst->currsol[j]] < min3) {

                if (inst->cost[last * n + inst->currsol[j]] < min2) {

                    if (inst->cost[last * n + inst->currsol[j]] < min) {

                        min3 = min2;
                        best_pos3 = best_pos2;
                        min2 = min;
                        best_pos2 = best_pos;
                        min = inst->cost[last * n + inst->currsol[j]];
                        best_pos = j;

                    }
                    else {

                        min3 = min2;
                        best_pos3 = best_pos2;
                        min2 = inst->cost[last * n + inst->currsol[j]];
                        best_pos2 = j;

                    }
                }
                else {

                    min3 = inst->cost[last * n + inst->currsol[j]];
                    best_pos3 = j;

                }
            }
        }
        double c = rand01();
        //  printf("random : %f\n", c);
        if (c < p3 && min3 < INF_DOUBLE) {
            // printf("sono nel caso dove prendo terzo migliore\n");
            costo = costo + min3;
            swap(inst->currsol, len, best_pos3);

            last = inst->currsol[len];
            len++;
        }
        else if (p3 <= c && c < p2 + p3 && min2 < INF_DOUBLE) {
            //printf("sono nel caso dove prendo secondo migliore\n");
            costo = costo + min2;
            swap(inst->currsol, len, best_pos2);

            last = inst->currsol[len];
            len++;
        }
        else {
            // printf("sono nel caso dove prendo primo migliore\n");
            costo = costo + min;
            swap(inst->currsol, len, best_pos);
            last = inst->currsol[len];
            len++;
        }
        min = INF_DOUBLE;
        min2 = INF_DOUBLE;
        min3 = INF_DOUBLE;
    }
    costo += inst->cost[start * n + last];
    inst->currcost = costo;
    save_if_best(inst);
}

// grasp where every 20 iteration it makes a random pick
void nearest_neighbor_grasp_random(Instance* inst, int start)
{
    double p2 = inst->solver.p1;
    if (start < 0)
        print_error("Invalid choice of the start node", __LINE__);

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
            int d = rand_int(len, n - 1);
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


void next_bestsol(Instance* inst, int it)
{

    //===============2-OPT CHECK=============
    int best_arc = -1;
    int best_arc2 = -1;
    double best_delta = INF_DOUBLE;
    int n = inst->nnodes;

    //check for tabu tenure validity
    if (inst->tabu_tenure < 0)
        print_error("%d, Impossible value of tabu tenure", __LINE__);

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
        print_error("%d, Error in finding another solution", __LINE__);

    //put all 4 nodes swapped in tabu list
    inst->tabu[inst->currsol[best_arc]] = it;
    inst->tabu[inst->currsol[best_arc + 1]] = it;
    inst->tabu[inst->currsol[best_arc2]] = it;
    inst->tabu[inst->currsol[best_arc2 + 1]] = it;

    invert_nodes(inst->currsol, best_arc + 1, best_arc2); // reverse the sub-tour between best_arc and best_arc2
    inst->currcost += best_delta;
   // inst->currcost = get_cost(inst, inst->currsol);// aggiunto questa riga per coerenza dei costi che sono sbagliati ma non capisco perch�
    // if (inst->currcost < inst->bestcost) //if new solution is best found, save it
     //{
     //    copy_array(inst->currsol, inst->bestsol, inst->nnodes + 1);
     //    inst->bestcost = inst->currcost;
      //   if (inst->verbose > 0)
       //      printf("\nNew best solution found of cost: %f", inst->bestcost);
    // }
    save_if_best(inst);
}
// cos� costi sono sbagliati no so perch�
void tabu_search(Instance* inst)
{ // inizializzo qua 
    nearest_neighbor_det(inst, 0);
    copy_array(inst->currsol, inst->bestsol, inst->nnodes + 1);
    inst->bestcost = inst->currcost;
    // two_opt(inst);
    printf("\n starting to iterate");
    //iterate untill timelimit is exceeded
    int it = 0;
    int time_elapsed = time(NULL) - inst->tstart;
    while (time_elapsed < inst->time_limit)
    {
        copy_array(inst->bestsol, inst->currsol, inst->nnodes + 1);
        inst->currcost = inst->bestcost;
        next_bestsol(inst, it);
        it++;
        time_elapsed = time(NULL) - inst->tstart;
    }
}



void two_opt(Instance* inst)
{
    int n = inst->nnodes;
    double best_delta = -1;
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

    }
    save_if_best(inst);
}



void simulated_annealing(Instance* inst) {
    int n = inst->nnodes;
    double T = INITIAL_TEMPERATURE; // initial temperature
    //printf("temp: %f", T);
    int i = 1;

    // initialize solution in currcost currcost � giusto in bestsol adesso c'e il two opt 
    nearest_neighbor_allstart(inst);
    copy_array(inst->bestsol, inst->currsol, inst->nnodes + 1);
    inst->currcost = inst->bestcost;
    two_opt(inst);
    inst->bestcost = get_cost(inst, inst->bestsol);

    double alpha = COOLING_FACTOR;// cooling factor
    time_t start_time = time(NULL);// inital time
    // int* newsol = (int*)calloc(n + 1, sizeof(int));// array for solution
    // double costnewsol = inst->bestcost;  // cost initialized 
   //  memcpy(newsol, inst->bestsol, (n + 1) * sizeof(int));//in new sol bestsol
    while (time(NULL) - start_time < inst->time_limit) {
        double costDiff = random_two_opt_move(inst, inst->currsol);// swap two cities   
        double deltaz = costDiff / ((inst->bestcost) / n);
        double c = rand01();

        double expo = exp(-deltaz / T);


        if ((costDiff < 0 || expo > c)) {
            printf("  costDiff: %f  ", costDiff);
            // Accept solution
          //  printf("accepted\n");
            memcpy(inst->bestsol, inst->currsol, (n + 1) * sizeof(int));
            inst->bestcost += costDiff;
            // inst->currcost = inst->bestcost;
            printf("\nSolution accepted: %f", inst->bestcost);

        }
        else {
            // printf("REJECTED");
             // Reject the new solution
            memcpy(inst->currsol, inst->bestsol, (n + 1) * sizeof(int));
            inst->currcost = inst->bestcost;
        }

        // Update the temperature
        T *= alpha;


    }
}


double random_two_opt_move(Instance* inst, int* sol) {
    int n = inst->nnodes;
    int second;
    int first;
    double delta;
    do {
        second = 2 + rand() % (inst->nnodes - 2);
        first = rand() % (second - 1);
        double dist_before = inst->cost[sol[first] * n + sol[first + 1]] + inst->cost[sol[second] * n + sol[second + 1]];
        double dist_after = inst->cost[sol[first] * n + sol[second]] + inst->cost[sol[first + 1] * n + sol[second + 1]];
        delta = dist_after - dist_before;
    } while (delta == 0);// this fixes a single case where 2 edges are consecutive

    invert_nodes(sol, first + 1, second);
    return delta;

}




// works better than 2 opt
void vns(Instance* inst) {
    nearest_neighbor_det(inst, 0);// initialize solution
    int start_time = time(NULL);
    int n = inst->nnodes;
    int i = 0;
    two_opt(inst);
    memcpy(inst->bestsol, inst->currsol, (n + 1) * sizeof(int));
    inst->bestcost = inst->currcost;
    while (time(NULL) - start_time < inst->time_limit) {
        memcpy(inst->currsol, inst->bestsol, (n + 1) * sizeof(int));
        kick(inst, inst->currsol);
        kick(inst, inst->currsol);
        inst->currcost = get_cost(inst, inst->currsol);
        two_opt(inst);
    }
}




int checkNonEqual(int* array, int tocheck, int length) {
    for (int i = 0; i < length; i++) {
        if (array[i] == tocheck)
            return 0;
    }
    return 1;
}




// gives the solution a kick
void kick(Instance* inst, int* sol) {
    int n = inst->nnodes;

    int* randomEdges = (int*)calloc(3, sizeof(int));

    randomEdges[0] = rand() % (n);
    // choses at random 3 nodes which represent the arcs to eliminate checking if they are equal and sorting them
    for (int i = 1; i < 3; i++) {
        double c = rand() % (n);
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