#include "TSP.h"  
#include "utility.h"
#include <time.h>

void initialize_instance(Instance* inst)
{
	inst->nnodes = -1;
	inst->bestsol = NULL;
	inst->cost = NULL;
	inst->randomseed = 1337;
}

void parse_args(Instance* inst, int argc, char** argv)
{
	if (argc < 2)
		print_error("Not enough arguments in the command line while launching the program");
	int arg = 2;
	for (int i = 1; i < argc; i++)
	{
		if (strncmp(argv[i], "-f", 2) == 0) //input file
		{
			strcpy(inst->inputfile, argv[++i]);
			continue;
		}
		if (strncmp(argv[i], "-v", 2) == 0) //verbose level
		{
			strcpy(inst->verbose, atoi(argv[++i]));
			continue;
		}
		if (strncmp(argv[i], "-seed", 5) == 0) //verbose level
		{
			inst->randomseed = atoi(argv[++i]);
			continue;
		}
	}
}

void parse_TSPLIB(Instance* inst)
{
	FILE* fin = fopen(inst->inputfile, "r");
	printf("%s", inst->inputfile);
	if (fin == NULL)
		print_error("input file not opened correctly");

	char line[180];

	//read dimension and set nnnodes in instance
	while (fgets(line, sizeof(line), fin) != NULL)
	{
		char* par = strtok(line, " :");
		if (strncmp(par, "NODE_COORD_SECTION", 18) == 0)
		{
			if (inst->nnodes < 0)
				print_error("Starting points coordinate before having the dimension");
			break;
		}
		if (strcmp(par, "DIMENSION", 9) == 0 )
		{
			if (inst->nnodes > 0)
				print_error("Dimension of dataset already set");
			int n = atoi(strtok(NULL, " :"));
			inst->nnodes = n; 
			continue;
		}
	}

	//allocate memory
	inst->points = CALLOC(inst->nnodes, Point);
	inst->cost = CALLOC(inst->nnodes * inst->nnodes, double);
	inst->bestsol = CALLOC(inst->nnodes + 1, int);

	//reading coordinates
	for (int i = 0; i < inst->nnodes; i++)
	{
		fgets(line, sizeof(line), fin);
		char* token = strtok(line, " ");	//skip index
		token = strtok(NULL, " ");
		double x = atof(token); 
		token = strtok(NULL, " ");
		double y = atof(token);
		Point p;
		p.x = x;
		p.y = y;
		inst->points[i] = p;
	}
}

void print_points(Instance* inst)
{
	printf("\nPOINTS:\n");
	for (int i = 0; i < inst->nnodes; i++)
		printf("%f, %f \n", inst->points[i].x, inst->points[i].y);
}

void print_error(const char* msg)
{
	printf("ERROR: %s", msg);
	exit(1);
}

Point* generate_random_points(int n)
{
	//allocate memory
	Point* points = CALLOC(n, Point);

	//generation
	double range_min = 0;
	double range_max = 10000;
	for (int i = 0; i < n; i++)
	{
		Point p;
		p.x = ((double)rand() / RAND_MAX) * (range_max - range_min) + range_min;
		p.y = ((double)rand() / RAND_MAX) * (range_max - range_min) + range_min;
		points[i] = p;
	}
	return points;
}

Point* generate_random_points_range(int n, double range_min, double range_max)
{
	//allocate memory
	Point* points = CALLOC(n, Point);

	//generation
	for (int i = 0; i < n; i++)
	{
		Point p;
		p.x = ((double)rand() / RAND_MAX) * (range_max - range_min) + range_min;
		p.y = ((double)rand() / RAND_MAX) * (range_max - range_min) + range_min;
		points[i] = p;
	}
	return points;
}

// creates a graph of the cycle on gnuplot
void plot_generator(Instance* inst) {
	FILE* out_lines;
	FILE* out_points;
	out_lines = fopen("../../Traveling-Salesman-Optimization/data/data_lines.dat", "w");
	out_points = fopen("../../Traveling-Salesman-Optimization/data/data_points.dat", "w");
	if (out_lines == NULL || out_points == NULL)
	{
		fclose(out_lines);
		fclose(out_points);
		print_error("Error in opening output data file");
	}

	for (int i = 0; i < inst->nnodes - 1; i++) {
		fprintf(out_points, "%f %f\n", inst->points[i].x, inst->points[i].y);
		if (inst->bestsol != NULL)
		{
			fprintf(out_lines, "%f %f\n", inst->points[inst->bestsol[i]].x, inst->points[inst->bestsol[i]].y);
			fprintf(out_lines, "%f %f\n\n\n", inst->points[inst->bestsol[i + 1]].x, inst->points[inst->bestsol[i + 1]].y);
		}
	}
	fprintf(out_points, "%f %f\n", inst->points[inst->nnodes-1].x, inst->points[inst->nnodes-1].y);
	if (!inst->bestsol == NULL)
	{
		fprintf(out_lines, "%f %f\n", inst->points[inst->bestsol[inst->nnodes - 1]].x, inst->points[inst->bestsol[inst->nnodes - 1]].y);
		fprintf(out_lines, "%f %f\n", inst->points[inst->bestsol[inst->nnodes]].x, inst->points[inst->bestsol[inst->nnodes]].y);
	}
	fclose(out_lines);
	fclose(out_points);

	chdir("../../Traveling-Salesman-Optimization/plot");
	if(inst->bestsol != NULL)
		system("gnuplot -persistent gp_lines.gp");
	else
		system("gnuplot -persistent gp_points.gp");
}

void free_instance(Instance* inst)
{
	free(inst->bestsol);
	free(inst->cost);
	free(inst->points);
}

double distance(Point* p1, Point* p2)
{
	double d = 0;
	d = sqrt(pow(p1->x - p2->x, 2) + pow(p1->y - p2->y, 2));
	return d;
}



//get_cost(i, j, inst)
//compute or return matrix value if already there