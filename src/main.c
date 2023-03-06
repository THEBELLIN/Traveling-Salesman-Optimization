#include "TSP.h"

int main(int argc, char** argv)
{
	Instance insta;
	Point p[6];
	int sol[] = { 0, 2, 4, 5,3,1 };
	int x[] = { 1, 7, 5, 6,11, 13 };
	int y[] = { 5, 4, 10, 2,8, 3 };

	Instance instance;
	initialize_instance(&instance);
	parse_args(&instance, argc, argv);
	parse_TSPLIB(&instance);
	print_points(&instance);
	instance.points = generate_random_points(instance.nnodes); //0 to 10k
	print_points(&instance);
	instance.points=generate_random_points_range(instance.nnodes, 0.0, 20000.0); //0 to 20k
	print_points(&instance);
	
	insta.bestsol = sol;
	insta.nnodes = 6;
	insta.randomseed = 72;
	for (int i = 0; i < 6; i++) {
		p[i].x = x[i];
		p[i].y = y[i];
	}
	insta.points = p;
	plot_generator(&insta);
	return 0;
}