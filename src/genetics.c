#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include "genetics.h"

void generate_random_tour(int* genes, const int nnodes)
{
	//initialization is deterministic
	for (int i = 0; i < nnodes+1; i++)
	{
		genes[i] = i;
	}

	//then swap at random 2 nodes
	for (int i = 0; i < nnodes+1; i++)
	{
		int i1 = rand_int(0, nnodes);
		int i2 = rand_int(0, nnodes);
		swap(genes, i1, i2);
	}
}

void choose_parents(individual* population, int* parents, const int nparents)
{
	//sort population in increasing fitness (decreasing number)
	qsort(population, POP_SIZE, sizeof(individual), compareIndividual);

	//normalize fitness values
	normalize_fitness(population);

	//draw randomly between 0 and 1
	double draw = rand01();

	//seek for drawn parent

	//THIS ONLY WORKS IF FITNESS IS INCREASING IN BEST INDIVIDUALS

}

void normalize_fitness(individual* population)
{
	double fit_sum = 0;
	for (int i = 0; i < POP_SIZE; i++)
	{
		fit_sum += population[i].fitness;
	}
	//set cumulative normalized values
	double cum_sum = 0;
	for (int i = 0; i < POP_SIZE; i++)
	{
		cum_sum += (population[i].fitness / fit_sum);
		population[i].normalized_fitness = cum_sum;
	}
	if (fabs(cum_sum - 1) > EPSILON)
	{
		print_error("%sError in cumulative sum fitness calculation", __LINE__);
	}
}

