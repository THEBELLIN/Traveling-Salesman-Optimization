#include "TSP.h"
#include <time.h>


void parse_TSPLIB(Instance*);
void print_error(char*);


void parse_TSPLIB(Instance* inst)
{
	FILE* fin = fopen(inst->inputfile, "r");
	if (fin == NULL)
		printf("File not found");

	char line[180];

	//read dimension and set nnnodes in instance
	while (fgets(line, sizeof(line), fin) != NULL && inst->nnodes<0)
	{
		char* cm = strtok(line, " :");
		if (strcmp(cm, "DIMENSION") == 0 )
		{
			if (inst->nnodes > 0)
				print_error("Dimension of dataset already set");
			int n = atoi(strtok(NULL, " :"));
			inst->nnodes = n; 
		}
	}

	//allocate memory
	inst->points = (Point*) calloc(inst->nnodes, sizeof(Point));

	//reading coordinates
	for (int i = 0; i < inst->nnodes; i++)
	{
		fgets(line, sizeof(line), fin);
		int x = strtok(line, " ");
		x = strtok(NULL, " "); //skip index
		int y = strtok(NULL, " ");
		Point p;
		p.x = x;
		p.y = y;
		inst->points[i] = p;
	}
}

void print_error(char* msg)
{
	printf("ERROR: %s", msg);
	exit(1);
}
