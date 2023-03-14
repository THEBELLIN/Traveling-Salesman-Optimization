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
	em_options o = NORMAL;							//NORMAL, GRASP2(p), GRASP3(p1, p2)
	extra_mileage(&instance, CONV_HULL, &o);		//RAND, MAX_DIST, CONV_HULL
	plot_generator(&instance);
	free_instance(&instance);
	return 0;
}