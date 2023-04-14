#ifndef GENETICS_H
#define GENETICS_H

#include "TSP.h"
#include "utility.h"
#include "heuristics.h"


//================HYPERPARAMETERS================
#define POP_SIZE 1000		//population size
#define AUTOMATIC_PASS 500  //number of best individuals that automatically pass the selection phase
#define N_PARENTS 500		//number of parents chosen at each iteration
#define N_CHILDREN 250		//number of children generated at each iteration
#define N_RAND 100			//number of randomly generated solutions at the start of the genetic algorithm
#define N_NN 400			//number of generated solutions with nearest neighbour at the start of the genetic algorithm
#define P_MUTATION 0.05		//probability of mutation

//================STRUCTS==================
typedef struct
{
	int* genes;
	double fitness; // 1000 / cost -> cost = 1000 / fitness
	double normalized_fitness;
}individual;

//================FUNCTIONS================
void genetic(Instance*, const int);
double get_fitness(Instance*, int*);
void generate_random_tour(int*, const int);
void choose_parents(individual*, int*, const int);
void normalize_fitness(individual*);
void generate_children(Instance*, individual*, int*, const int, const int, const int);
void crossover(Instance*, individual*, const int, const int, const int);
void mutation(Instance*, individual*, const int, const int);
void selection(individual*, const int, const int);

#endif
