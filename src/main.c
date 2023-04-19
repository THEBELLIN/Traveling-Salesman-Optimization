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
	
	/*
	//solve options
	em_options oo = NORMAL;				//NORMAL, GRASP2(p), GRASP3(p1, p2)
	oo.start = CONV_HULL;
	solve_options o =
	{
		.alg = GEN,
		.timelimit = 10, 
		.em_opts = &oo
	};
	//extra_mileage(&instance, &o);		//RAND, MAX_DIST, CONV_HULL
	solve(&instance, &o);
	*/
	TSPopt(&instance);
	//plot_generator(&instance, instance.nnodes);
	free_instance(&instance);
	return 0;
}