#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include "genetics.h"

void genetic(Instance* inst)
{
	//first 100 are population, 
	individual* population = MALLOC(POP_SIZE + N_CHILDREN, individual);
	int added = 0;

	if (inst->currcost < INF_INT)
	{
		//add currsol and bestsol to current population, if any
		population[0].genes = MALLOC(inst->nnodes + 1, int);
		population[1].genes = MALLOC(inst->nnodes + 1, int);
		for (int i = 0; i < inst->nnodes + 1; i++)
		{
			population[0].genes[i] = inst->currsol[i];
			population[1].genes[i] = inst->bestsol[i];
		}
		population[0].fitness = get_fitness(inst, population[0].genes);
		population[1].fitness = get_fitness(inst, population[1].genes);
		added = 2;
	}

	//populate 8 (or 10) more at random
	for (int i = added; i < N_RAND; i++)
	{
		population[i].genes = MALLOC(inst->nnodes + 1, int);
		generate_random_tour(population[i].genes, inst->nnodes);
		population[i].fitness = get_fitness(inst, population[i].genes);
	}

	//populate other 90 using intelligient solutions
	//40 nearest neighbor grasp 2
	for (int i = N_RAND; i < N_RAND + N_NN; i++)
	{
		population[i].genes = MALLOC(inst->nnodes + 1, int);
		int starting_node = rand_int(0, inst->nnodes - 1);
		inst->solver.p1 = 0.15;
		inst->solver.p2 = 0.05;
		nearest_neighbor_grasp2(inst, starting_node);
		copy_array(inst->currsol, population[i].genes, inst->nnodes + 1);
		population[i].fitness = get_fitness(inst, population[i].genes);
	}

	//50 extramileage grasp 2 with p=0.1
	for (int i = N_RAND + N_NN; i < POP_SIZE; i++)
	{
		population[i].genes = MALLOC(inst->nnodes + 1, int);
		inst->solver.id = EM_GRASP2;
		inst->solver.start = RAND;
		inst->solver.p1 = 0.1;
		inst->solver.p2 = 0;
		extra_mileage_grasp2(inst);

		copy_array(inst->currsol, population[i].genes, inst->nnodes + 1);
		population[i].fitness = get_fitness(inst, population[i].genes);
	}

	//allocate space for children
	for (int i = POP_SIZE; i < POP_SIZE + N_CHILDREN; i++)
	{
		population[i].genes = MALLOC(inst->nnodes + 1, int);
	}

	//genetic algorithm
	int it = 0;
	while ((time(NULL) - inst->tstart) < inst->time_limit)
	{
		//printf("\niteration %d", it);
		int* parents = MALLOC(N_PARENTS, int);
		//printf("\nchoosing parents");
		choose_parents(population, parents, N_PARENTS);
		//printf("\ngenerating children");
		generate_children(inst, population, parents, N_PARENTS, N_CHILDREN, inst->nnodes);
		//printf("\nselecting survivors");
		selection(population, POP_SIZE + N_CHILDREN, POP_SIZE);
		individual* champ = get_champion(population, POP_SIZE);
		printf("\niteration: %d, champion fitness: %f", it, 1000 / champ->fitness);
		copy_array(champ->genes, inst->currsol, inst->nnodes + 1);
		save_if_best(inst);
		free(parents);
		it++;
	}

	//free
	for (int i = 0; i < POP_SIZE + N_CHILDREN; i++)
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
		int i1 = rand_int(0, nnodes - 1);
		int i2 = rand_int(0, nnodes - 1);
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

//generates children after population in the same array
void generate_children(Instance* inst, individual* population, int* parents, const int nparents, const int nchild, const int nnodes)
{
	int i = 0;
	for (int c = 0; c < nchild; c++)
	{
		double p = rand01();
		if (p < P_MUTATION)
		{
			//printf("\nmutation");
			mutation(inst, population, i, c);
			i++;
		}
		else
		{
			int j = (i + 1) % nparents; //%for odd number of parents
			i = (i + 2) % nparents;
			//printf("\ncrossover between %d and %d", parents[i], parents[j]);
			crossover(inst, population, parents[i], parents[j], c);
		}
	}
}

void crossover(Instance* inst, individual* population, const int parent1, const int parent2, const int n_child)
{
	individual p1 = population[parent1];
	individual p2 = population[parent2];
	individual* child = population + POP_SIZE + n_child;

	//====================CUTOFF METHOD=====================
	int cutoff = rand_int(0, inst->nnodes);
	int* visited = CALLOC(inst->nnodes, int);
	int index = 0;

	for (int i = 0; i < inst->nnodes; i++)
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
	while (index < inst->nnodes)
	{
		for (int i = 0; i <= cutoff; i++)
		{
			if (!visited[p2.genes[i]])
			{
				child->genes[index++] = p2.genes[i];
			}
		}
	}
	//close the tour
	child->genes[inst->nnodes] = child->genes[0];
	child->fitness = get_fitness(inst, child->genes);

	free(visited);

	//method2: add remaining nodes with extra mileage
}

void mutation(Instance* inst, individual* population, const int parent, const int n_child)
{	//call with 39, 20 causes pop[9] to get fucked up
	copy_array(&population[parent].genes[0], &population[POP_SIZE + n_child].genes[0], inst->nnodes + 1);
	kick(inst, &(population[POP_SIZE + n_child].genes[0]));
	population[POP_SIZE + n_child].fitness = get_fitness(inst, &population[POP_SIZE + n_child].genes[0]);
}

individual* get_champion(individual* population, const int size)
{
	qsort(population, size, sizeof(individual), compareIndividual);
	return population + (size - 1);
}

void selection(individual* population, const int all_size, const int desired_size)
{
	if (all_size <= desired_size)
	{
		print_error("%s Error in resizing population in selection", __LINE__);
	}

	qsort(population, all_size, sizeof(individual), compareIndividual_rev);

	//first AUTOMATIC_PASS get automatically selected
	//other at random
	for (int i = AUTOMATIC_PASS; i < desired_size; i++)
	{
		int j = rand_int(i, all_size - 1);
		swap_individual(population, i, j);
	}
}