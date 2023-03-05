#include "TSP.h"

int main(int argc, char** argv)
{
	Instance instance;
	initialize_instance(&instance);
	parse_args(&instance, argc, argv);
	parse_TSPLIB(&instance);
	print_points(&instance);
	return 0;
}