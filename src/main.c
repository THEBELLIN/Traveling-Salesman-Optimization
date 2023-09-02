#include "TSP.h"
#include "utility.h"
#include "heuristics.h"
#include "cplex_model.h"

#pragma comment(lib, "cplex2211.lib")
#pragma comment(lib, "ilocplex.lib")


int main(int argc, char** argv)
{
	Instance instance;
	initialize_instance(&instance);
	parse_args(&instance, argc, argv);
	srand(instance.randomseed);
	parse_TSPLIB(&instance);
	initialize_cost(&instance);
	print_params(&instance);
	solve(&instance);
	//plot_generator(&instance, instance.nnodes);
	free_instance(&instance);
	return 0;
}