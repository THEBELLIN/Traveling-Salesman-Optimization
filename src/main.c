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
	/*em_options o = GRASP3(0.7, 0.2);							//NORMAL, GRASP2(p), GRASP3(p1, p2)
	extra_mileage(&instance, CONV_HULL, &o);		//RAND, MAX_DIST, CONV_HULL
	plot_generator(&instance);*/
	int hsize;
	Point* hull = convex_hull(&instance, &hsize);
	/*FILE* out_lines;
	out_lines = fopen("../../Traveling-Salesman-Optimization/data/data_lines.dat", "w");*/
	for (int i = 0; i < hsize-1; i++)
	{
		printf("%f %f\n", instance.points[instance.bestsol[i]].x, instance.points[instance.bestsol[i]].y);
		printf("%f %f\n\n\n", instance.points[instance.bestsol[i + 1]].x, instance.points[instance.bestsol[i + 1]].y);
		//fprintf(out_lines, "%f %f\n", instance.points[instance.bestsol[i]].x, instance.points[instance.bestsol[i]].y);
		//fprintf(out_lines, "%f %f\n\n\n", instance.points[instance.bestsol[i + 1]].x, instance.points[instance.bestsol[i + 1]].y);
	}
	printf("%f %f\n", instance.points[instance.bestsol[hsize - 1]].x, instance.points[instance.bestsol[hsize - 1]].y);
	printf("%f %f\n\n\n", instance.points[instance.bestsol[0]].x, instance.points[instance.bestsol[0]].y);
	/*fprintf(out_lines, "%f %f\n", instance.points[instance.bestsol[hsize - 1]].x, instance.points[instance.bestsol[hsize - 1]].y);
	fprintf(out_lines, "%f %f\n\n\n", instance.points[instance.bestsol[0]].x, instance.points[instance.bestsol[0]].y);
	chdir("../../Traveling-Salesman-Optimization/plot");
	system("gnuplot -persistent gp_lines.gp");
	fclose(out_lines);*/
	free_instance(&instance);
	return 0;
}