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
	em_options o = NORMAL;					//deterministic extra_mileage
	extra_mileage(&instance, RAND, &o);		//random starting points
	extra_mileage(&instance, MAX_DIST, &o); //max distance starting points
	plot_generator(&instance);
	free_instance(&instance);
	return 0;
}