#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include "genetics.h"


void genetic(Instance* inst, const int niter)
{
	individual* population = MALLOC(POP_SIZE, individual);

	//add currsol and bestsol to current population
	population[0].genes = MALLOC(inst->nnodes + 1, int);
	population[1].genes = MALLOC(inst->nnodes + 1, int);
	for (int i = 0; i < inst->nnodes + 1; i++)
	{
		population[0].genes[i] = inst->currsol[i];
		population[1].genes[i] = inst->bestsol[i];
	}
	population[0].fitness = get_fitness(inst, population[0].genes);
	population[1].fitness = get_fitness(inst, population[1].genes);

	//populate 98 more at random
	for (int i = 2; i < 100; i++)
	{
		population[i].genes = MALLOC(inst->nnodes + 1, int);
		generate_random_tour(population[i].genes, inst->nnodes);
		population[i].fitness = get_fitness(inst, population[i].genes);
	}

	//populate other 900 using intelligient solutions
	for (int i = 100; i < 500; i++)
	{
		population[i].genes = MALLOC(inst->nnodes + 1, int);
		//TODO NN_GRASP2
		population[i].fitness = get_fitness(inst, population[i].genes);
	}

	for (int i = 500; i < POP_SIZE; i++)
	{
		population[i].genes = MALLOC(inst->nnodes + 1, int);
		extra_mileage_grasp2(inst, RAND, 0.1);
		population[i].fitness = get_fitness(inst, population[i].genes);
	}

	//free
	for (int i = 0; i < POP_SIZE; i++)
	{
		free(population[i].genes);
	}
	free(population);
 }

double get_fitness(Instance* inst, int* genes)
{
	double cost = get_cost(inst, genes);
	return 1000 / cost;
}

void generate_random_tour(int* genes, const int nnodes)
{
	//initialization is deterministic
	for (int i = 0; i < nnodes; i++)
	{
		genes[i] = i;
	}

	//then swap at random 2 nodes
	for (int i = 0; i < nnodes; i++)
	{
		int i1 = rand_int(0, nnodes-1);
		int i2 = rand_int(0, nnodes-1);
		swap(genes, i1, i2);
	}

	//close loop
	genes[nnodes] = genes[0];
}

//select parents based on roulette-wheel selection
void choose_parents(individual* population, int* parents, const int nparents)
{
	//sort population in increasing fitness (decreasing number)
	qsort(population, POP_SIZE, sizeof(individual), compareIndividual);

	//normalize fitness values
	normalize_fitness(population);

	int selected[POP_SIZE];
	memset(selected, 0, POP_SIZE * sizeof(int));
	int chosen_parents = 0;

	while (chosen_parents != nparents)
	{
		//draw randomly between 0 and 1
		double draw = rand01();

		//seek for drawn parent
		for (int i = 0; i < POP_SIZE; i++)
		{
			//first one to be higher gets chosen
			if (population[i].normalized_fitness >= draw)
			{
				if (selected[i] == 0)
				{
					selected[i] = 1;
					chosen_parents++;
				}
				break;
			}
		}
	}

	//copy selected parents to parents array
	int t = 0;
	for (int i = 0; i < POP_SIZE; i++)
	{
		if (selected[i])
		{
			parents[t] = i;
			t++;
		}
	}
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

void generate_children(individual* population, int* parents, const int nparents, const int nnodes)
{
	int n_child = 0;
	for (int i = 0; i < nparents; i = i + 2)
	{
		int j = (i + 1) % nparents; //%for odd number of parents
		crossover(population, parents[i], parents[j], n_child, nnodes);
		n_child++;
	}
}

void crossover(individual* population, const int parent1, const int parent2, const int n_child, const int nnodes)
{
	individual p1 = population[parent1];
	individual p2 = population[parent2];
	individual* child = population + POP_SIZE + n_child;

	//====================CUTOFF METHOD=====================
	int cutoff = rand_int(0, nnodes);
	int* visited = CALLOC(nnodes, int);
	int index = 0;

	for (int i = 0; i < nnodes; i++)
	{
		//before cutoff copy p1 tour
		if (i <= cutoff)
		{
			child->genes[i] = p1.genes[i];
			visited[p1.genes[i]] = 1;
		}

		//after cutoff, copy nonvisited nodes from p2
		else
		{
			if (visited[p2.genes[i]])
			{
				continue;
			}
			child->genes[index] = p2.genes[i];
			visited[p2.genes[index]] = 1;
		}
		index++;
	}

	//fill tour with missing nodes
	while (index < nnodes)
	{
		for (int i = 0; i <= cutoff; i++)
		{
			if (!visited[p2.genes[i]])
			{
				child->genes[index++] = p2.genes[i];	
			}
		}
	}
	free(visited); 

	//method2: add remaining nodes with extra mileage TODO
}

individual* get_champion(individual* population, const int size)
{
	qsort(population, size, sizeof(individual), compareIndividual);
	return &population[size - 1];
}

//MUTATION 1PARENT->1CHILD WITH KICK

//KILL TO POP_COUNT

//POPULATE START

//DEFINE HYPERAPRAMS