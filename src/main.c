#include "TSP.h"

int main(int argc, char** argv)
{
	Instance instance;
	instance.inputfile = "../data/att48.tsp";
	parse_TSPLIB(&instance);
	printf("%d", instance.nnodes);
	return 0;
}