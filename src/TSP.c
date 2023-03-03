#include "TSP.h"
#include <time.h>


void parse_TSPLIB(Instance*);


void parse_TSPLIB(Instance* inst)
{
	FILE* fin = fopen(inst->inputfile, "r");
	if (fin == NULL)
		printf("File not found");

	char line[180];

	while (fgets(line, sizeof(line), fin) != NULL)
	{
		char* cm = strtok(line, " :");
		if (strcmp(cm, "DIMENSION") == 0)
		{
			int n = atoi(strtok(NULL, " :"));
			inst->nnodes = n; 
			printf("%d",n);
		}
	}
}


