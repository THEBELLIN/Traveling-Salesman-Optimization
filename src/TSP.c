#include "TSP.h"  
#include "utility.h"
#include "heuristics.h"
#include "annealing.h"
#include "cplex_model.h"
#include "genetics.h"
#include "local_branching.h"
#include "hard_fixing.h"
#include <time.h>

//check feasibility of a solution with n nodes
bool check_feasibility(int* sol, int n)
{
	//create and populate counters
	int* counters = CALLOC(n, int);
	for (int i = 0; i < n; i++)
	{
		counters[sol[i]]++;
	}

	//check that it is a cycle
	if (sol[0] != sol[n])
	{
		free(counters);
		return 0;
	}

	//check all counters are 1
	for (int i = 0; i < n; i++)
	{
		if (*(counters+i) != 1)
		{
			free(counters);
			return 0;
		}
	}
	free(counters);
	return 1;
}

void initialize_instance(Instance* inst)
{
	inst->nnodes = -1;
	inst->bestsol = NULL;
	inst->currsol = NULL;
	inst->cost = NULL;
	inst->currcost = INF_INT;
	inst->bestcost = INF_INT;
	inst->randomseed = 1337;
	inst->tstart = time(NULL);
	inst->tabu_tenure = -1;
}

void parse_args(Instance* inst, int argc, char** argv)
{
	//TODO check all fundamental arguments are filled
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
			inst->verbose = atoi(argv[++i]);
			continue;
		}
		if (strncmp(argv[i], "-seed", 5) == 0) //verbose level
		{
			inst->randomseed = atoi(argv[++i]);
			continue;
		}
		if (strncmp(argv[i], "-tenure", 7) == 0) //tabu tenure level
		{
			inst->tabu_tenure = atoi(argv[++i]);
			continue;
		}
		if (strncmp(argv[i], "-time", 5) == 0) //time limit
		{
			inst->time_limit = atoi(argv[++i]);
			continue;
		}
	}
}

void parse_TSPLIB(Instance* inst)
{
	FILE* fin = fopen(inst->inputfile, "r");
	printf("%d", inst->inputfile);
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
	inst->points = MALLOC(inst->nnodes, Point);
	inst->cost = MALLOC(inst->nnodes * inst->nnodes, double);
	inst->bestsol = MALLOC(inst->nnodes + 1, int);
	inst->currsol = MALLOC(inst->nnodes + 1, int);
	inst->tabu = MALLOC(inst->nnodes, int);

	//initialize tbau list
	for (int i = 0; i < inst->nnodes; i++)
	{
		inst->tabu[i] = -INF_INT;
	}

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
	printf("ERROR: %d", msg);
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
/*void plot_generator(Instance* inst) {
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
}*/

void plot_generator(Instance* inst, int n_edges) {
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
	//print points
	for (int i = 0; i < inst->nnodes - 1; i++) {
		fprintf(out_points, "%f %f \"%d\"\n", inst->points[i].x, inst->points[i].y, i);
	}
	fprintf(out_points, "%f %f \"%d\"\n", inst->points[inst->nnodes - 1].x, inst->points[inst->nnodes - 1].y, inst->nnodes - 1);
	
	//print lines
	for (int i = 0; i < n_edges; i++)
	{
		fprintf(out_lines, "%f %f\n", inst->points[inst->bestsol[i]].x, inst->points[inst->bestsol[i]].y);
		fprintf(out_lines, "%f %f\n\n\n", inst->points[inst->bestsol[i + 1]].x, inst->points[inst->bestsol[i + 1]].y);
		printf("\n%d: plotting line from %d to %d", i, inst->bestsol[i], inst->bestsol[i + 1]);
	}
	
	fclose(out_lines);
	fclose(out_points);

	chdir("../../Traveling-Salesman-Optimization/plot");
	system("gnuplot -persistent gp_points_and_lines.gp");
}

void free_instance(Instance* inst)
{
	free(inst->bestsol);
	free(inst->cost);
	free(inst->points);
	free(inst->tabu);
	free(inst->currsol);
}

double distance(Point* p1, Point* p2)
{
	double d = 0;
	d = sqrt((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y));
	return d;
}

void print_points_file(Point* points, int n, FILE* out)
{
	fprintf(out, "DIMENSION : %d\n", n);
	for (int i = 0; i < n; i++)
		fprintf(out, "%d %f %f \n", i + 1, points[i].x, points[i].y);
}

Instance* generate_test_bed(int seed, int n_instances, int n_points) {
	FILE* out = fopen("../../Traveling-Salesman-Optimization/data/test_bed.txt", "w");
	if (out == NULL)
		print_error("Error in opening output data file");
	Instance* inst_set = CALLOC(n_instances, Instance);
	srand(seed);
	for (int i = 0; i < n_instances; i++) {
		inst_set[i].randomseed = seed;
		inst_set[i].nnodes = n_points;
		inst_set[i].points = generate_random_points(n_points);
		strcpy(inst_set[i].inputfile, "test_bed.txt");
		print_points_file(inst_set[i].points, n_points, out);
	}
	fprintf(out, "EOF");
	fclose(out);
	return inst_set;
}

double get_cost(Instance* inst, int* sol)
{
	double cost = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		cost += distance(&inst->points[sol[i]], &inst->points[sol[i + 1]]);
	}
	return cost;
}

//checks if currsol is better than bestsol and saves it if that's the case
//verbose > 1 print every best solution
void save_if_best(Instance* inst)
{
	if (!check_feasibility(inst->currsol, inst->nnodes))
	{
		print_error("%d, Unfeasible solution found", __LINE__);
	}
	if (inst->currcost < inst->bestcost)
	{
		copy_array(inst->currsol, inst->bestsol, inst->nnodes + 1);
		inst->bestcost = inst->currcost;
		if (inst->verbose > 1)
		{
			printf("\nNew best solution found of cost: %f", inst->bestcost);
		}
	}
}

void file_perf_prof(int n_instances, int n_points, int seed) {
	// generate n_instances random instances
	FILE* out = fopen("../data/PerfProf.txt", "w");
	if (out == NULL)
		print_error("Error in opening output data file");
	Instance* set = generate_test_bed(seed, n_instances, n_points);
	fprintf(out, "%d nearest_neighbor allstart\n", 2);
	for (int i = 0; i < n_instances; i++) {
		nearest_neighbor(&set[i], 0);
		fprintf(out, "Instance%d %f ", i + 1, set[i].bestcost);
		nearest_neighbor_allstart(&set[i]);
		fprintf(out, "%f\n", set[i].bestcost);
	}
}

void sol_to_tsp(int* sol, double* cplex, Instance* inst) {
	for (int i = 0; i < inst->ncols; i++) {
		cplex[i] = 0.0;
	}
	for (int i = 0; i < inst->nnodes; i++) {
		cplex[xpos(sol[i], sol[i + 1], inst)] = 1.0;
	}

}
// works
// transforms a solution in the format of successors in a permutation and saves it in best sol
void transform_in_perm_and_save(int* succ, Instance* inst) {
	inst->bestsol[0] = 0;
	int curr = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		inst->bestsol[i + 1] = succ[curr];
		curr = succ[curr];
	}
}

//TODO check all parameters needed are set in parseargs
void solve(Instance* inst)
{
	solver_id sid = inst->solver.id;
	if (sid == EM || sid == EM_GRASP2 || sid == EM_GRASP3)
	{
		extra_mileage(inst);
	}
	else if (sid == NN || sid == NN_GRASP2 || sid == NN_GRASP3)
	{
		nearest_neighbor(inst);
	}
	else if (sid == GEN)
	{
		genetic(inst);
		return;
	}
	else if (sid == TABU)
	{
		tabu_search(inst);
	}
	else if (sid == VNS)
	{
		vns(inst);
	}
	else if (sid == SIMANN)
	{
		simulated_annealing(inst); 
	}
	else if (sid == CPLEX_BENDERS || sid == CPLEX_CALLBACK)
	{
		TSPopt(inst);
	}
	else if (sid == LOCAL_BRANCHING)
	{
		local_branching(inst);
	}
	else if (sid == HARD_FIXING )
	{
		hard_fixing(inst);
	}
	else
		print_error("%d, Error in setting algorithm options for solving", __LINE__);


	//plot best instance found
	//plot_generator(inst, inst->nnodes);
}

//plot cost of incumbent during solving (or cost in console but better to plot)
//tabulate formatted output