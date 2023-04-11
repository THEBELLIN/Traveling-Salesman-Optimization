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
	double fitness; //in our case this is actually the opposite. AKA the less fit, the better
	double normalized_fitness;
}individual;

//================FUNCTIONS================
void generate_random_tour(int*, const int);
void choose_parents(individual*, int*, const int);
void normalize_fitness(individual*);

#endif
