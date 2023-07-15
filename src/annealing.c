#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include "annealing.h"

void simulated_annealing(Instance* inst)
{
    int time_limit = inst->time_limit;
    int n = inst->nnodes;
    double T = 100.0; // initial temperature
    printf("temp: %f", T);
    int i = 1;
    double alpha = 0.99;// cooling factor
    time_t start_time = time(NULL);// inital time
    int* newsol = (int*)calloc(n + 1, sizeof(int));// array for solution
    double costnewsol = inst->bestcost;  // cost initialized 
    memcpy(newsol, inst->bestsol, (n + 1) * sizeof(int));//in new sol bestsol
    while (time(NULL) - start_time < time_limit) {
        random_two_opt_move(newsol, inst);// swap two cities
        double costDiff = get_cost(inst, newsol) - costnewsol;
        //printf("  costDiff= %f  ", costDiff);
        double deltaz = costDiff / ((inst->bestcost) / n);
        // printf("  deltaz = %f  ", deltaz);
        double c = rand01();
        // printf("  random: %f  ", c);
        double expo = exp(-deltaz / T);
        if (i == 1)
            printf(" exp=  %f\n  ", expo);


        i = 2;
        if (costDiff < 0 || expo > c) {
            //printf("ACCEPTED");
            // Accept solution
            memcpy(inst->bestsol, newsol, (n + 1) * sizeof(int));
            inst->bestcost += costDiff;
            costnewsol = inst->bestcost;
        }
        else {
            //printf("REJECTED");
            // Reject the new solution
            memcpy(newsol, inst->bestsol, (n + 1) * sizeof(int));
            costnewsol = inst->bestcost;
        }
        // Update the temperature
        T *= alpha;
    }
    free(newsol);
}


void random_two_opt_move(int* sol, Instance* inst) {
    int n = inst->nnodes;
    int first = rand() % (n - 4);
    int second = first + (rand() % (n - 2 - first + 1));
    invert_nodes(sol, first + 1, second);
}

