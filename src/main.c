#include "TSP.h"

int main(int argc, char** argv)
{
	Instance instance;
	initialize_instance(&instance);
	parse_args(&instance, argc, argv);
	parse_TSPLIB(&instance);
	print_points(&instance);
	instance.points = generate_random_points(instance.nnodes); //0 to 10k
	print_points(&instance);
	instance.points=generate_random_points_range(instance.nnodes, 0.0, 20000.0); //0 to 20k
	print_points(&instance);
	return 0;
}