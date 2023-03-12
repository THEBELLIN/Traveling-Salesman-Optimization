#include "TSP.h"
#include "utility.h"
#include "heuristics.h"

int main(int argc, char** argv)
{
	Instance instance;
	initialize_instance(&instance);
	parse_args(&instance, argc, argv);
	srand(instance.randomseed);
	parse_TSPLIB(&instance);
	initialize_cost(&instance);
	em_options o = NORMAL;
	extra_mileage(&instance, MAX_DIST, &o);
	system("cd");
	plot_generator(&instance);
	free_instance(&instance);
	return 0;
}