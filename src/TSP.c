#include "TSP.h"
#include <time.h>


void parse_TSPLIB(Instance*);
void print_error(char*);
void initialize_instance(Instance*);
void parse_args(Instance*, int, char**);
void print_points(Instance*);
Point* generate_random_points(int);
Point* generate_random_points_range(int, double, double);

void initialize_instance(Instance* inst)
{
	inst->nnodes = -1;
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
	inst->points = (Point*) calloc(inst->nnodes, sizeof(Point));

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

void print_error(char* msg)
{
	printf("ERROR: %s", msg);
	exit(1);
}

Point* generate_random_points(int n)
{
	//allocate memory
	Point* points = (Point*) calloc(n, sizeof(Point));

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
	Point* points = (Point*)calloc(n, sizeof(Point));

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