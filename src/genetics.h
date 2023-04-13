#ifndef GENETICS_H
#define GENETICS_H

#include "TSP.h"
#include "utility.h"
#include "heuristics.h"


//================CONSTANTS================
#define POP_SIZE 1000

//================STRUCTS==================
typedef struct
{
	int* genes;
	double fitness; //1000/cost or something like that
	double normalized_fitness;
}individual;

//================FUNCTIONS================
void genetic(Instance*, const int);
double get_fitness(Instance*, int*);
void generate_random_tour(int*, const int);
void choose_parents(individual*, int*, const int);
void normalize_fitness(individual*);
void generate_children(individual*, int*, const int, const int);
void crossover(individual*, const int, const int, const int, const int);

#endif
