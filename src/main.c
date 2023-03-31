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

	//solve options
	em_options oo = NORMAL;				//NORMAL, GRASP2(p), GRASP3(p1, p2)
	oo.start = CONV_HULL;
	solve_options o =
	{
		o.alg = EM,
		o.timelimit = 300, //5 mins
		o.em_opts = &oo
	};
	//extra_mileage(&instance, &o);		//RAND, MAX_DIST, CONV_HULL
	//plot_generator(&instance);
	solve(&instance, &o);
	free_instance(&instance);
	return 0;
}