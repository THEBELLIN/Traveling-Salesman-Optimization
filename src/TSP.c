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
	inst->solver.id = -1;
	inst->solver.p1 = -1;
	inst->solver.p2 = -1;
	inst-> solver.start = -1;
	inst->verbose = 0;
	inst->time_limit = -1;
}

void parse_args(Instance* inst, int argc, char** argv)
{
	if (argc < 2)
		print_error("Not enough arguments in the command line while launching the program", __LINE__);
	int arg = 2;
	for (int i = 1; i < argc; i++)
	{
		if (strncmp(argv[i], "-help", 5) == 0) //parameters help and after it finishes the program
		{
			//print instructions
			printf("Usage: %s -f <inputfile> -v <verbose> -seed <randomseed> -tenure <tabu_tenure> -time <time_limit> -solver <solver> -p1 <p1> -p2 <p2> -em_start <em_start>\n", argv[0]);

			//print solvers
			printf("Available solvers: NN, GRASP2_NN, GRASP3_NN, MILEAGE, EM_GRASP2, EM_GRASP3, GEN, VNS, SIMANN, TABU, CPLEX_BENDERS, CPLEX_CALLBACK, LOCAL_BRANCHING, HARD_FIXING\n");

			//print available em_start
			printf("Available em_starts: RAND, MAX_DIST, CONV_HULL\n");

			//print paramenters needed
			printf("Solvers that need p1: GRASP2_NN, GRASP3_NN, EM_GRASP2, EM_GRASP3, VNS, LOCAL_BRANCHING, HARD_FIXING\n");
			printf("Solvers that also need p2: GRASP3_NN, EM_GRASP3, HARD_FIXING\n");
			printf("Solvers that need em_start: MILEAGE, EM_GRASP2, EM_GRASP3\n");

			//print default values
			printf("Default values: verbose: 0, seed: 1337.\n");
			break;
		}
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
		if (strncmp(argv[i], "-seed", 5) == 0) //seed
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
		if (strncmp(argv[i], "-p1", 3) == 0) //p1
		{
			inst->solver.p1 = atof(argv[++i]);
			//check that p1 is in the range [0,1]
			if (inst->solver.p1 < 0 || inst->solver.p1 > 1)
				print_error("%d, p1 must be in the range [0,1]", __LINE__);
			continue;
		}
		if (strncmp(argv[i], "-p2", 3) == 0) //p2
		{
			inst->solver.p2 = atof(argv[++i]);
			//check that p2 is in the range [0,1]
			if (inst->solver.p2 < 0 || inst->solver.p1 > 1)
				print_error("%d, p1 must be in the range [0,1]", __LINE__);
			continue;
		}
		if (strncmp(argv[i], "-em_start", 9) == 0) //em_start
		{
			inst->solver.start = atoi(argv[++i]);
			//check that em_start is one of the following values: 0, 1, 2 
			if (inst->solver.start < 0 || inst->solver.start > 2)
				print_error("%d, em_start must be one of the following values: RAND, MAX_DIST, CONV_HULL", __LINE__);

			continue;
		}
		if (strncmp(argv[i], "-solver", 7) == 0) //solver
		{
			if (strncmp(argv[++i], "nn", 2) == 0)
				inst->solver.id = NN;
			else if (strncmp(argv[i], "grasp2_nn", 9) == 0)
				inst->solver.id = GRASP2_NN;
			else if (strncmp(argv[i], "grasp3_nn", 9) == 0)
				inst->solver.id = GRASP3_NN;
			else if (strncmp(argv[i], "mileage", 7) == 0)
				inst->solver.id = MILEAGE;
			else if (strncmp(argv[i], "em_grasp2", 9) == 0)
				inst->solver.id = EM_GRASP2;
			else if (strncmp(argv[i], "em_grasp3", 9) == 0)
				inst->solver.id = EM_GRASP3;
			else if (strncmp(argv[i], "two_opt", 7) == 0)
				inst->solver.id = TWO_OPT;
			else if (strncmp(argv[i], "all_start_nn", 12) == 0)
				inst->solver.id = ALL_START_NN;
			else if (strncmp(argv[i], "vns", 3) == 0) {
				inst->solver.id = VNS;
				//printf("fatto");
			}
			else if (strncmp(argv[i], "genetic", 7) == 0)
				inst->solver.id = GEN;
			else if (strncmp(argv[i], "simann", 6) == 0)
				inst->solver.id = SIMANN;
			else if (strncmp(argv[i], "tabu", 4) == 0)
				inst->solver.id = TABU;
			else if (strncmp(argv[i], "cplex_benders", 13) == 0)
				inst->solver.id = CPLEX_BENDERS;
			else if (strncmp(argv[i], "cplex_callback", 14) == 0)
				inst->solver.id = CPLEX_CALLBACK;
			else if (strncmp(argv[i], "local_branching", 15) == 0)
				inst->solver.id = LOCAL_BRANCHING;
			else if (strncmp(argv[i], "hard_fixing", 11) == 0)
				inst->solver.id = HARD_FIXING;
			else
				print_error("solver not recognized", __LINE__);
			continue;
		}
	}
}


void print_params(Instance* inst) {
	printf("\ninput file:%s", inst->inputfile);
	printf("\nnnodes:%d", inst->nnodes);
	printf("\ntabu_tenure:%d", inst->tabu_tenure);
	printf("\ntime_limit:%f", inst->time_limit);
	printf("\nverbose: %d", inst->verbose);
	printf("\np2: %f", inst->solver.p1);
	printf("\np3: %f", inst->solver.p2);
}

void check_params(Instance* inst)
{
	//check that input file is set
	if (inst->inputfile[0] == '\0')
		print_error("%d, Input file not set", __LINE__);

	//check that time limit is set
	if (inst->time_limit < 0)
		print_error("%d, Time limit not set", __LINE__);

	//check that solver is set
	if (inst->solver.id == NULL)
		print_error("%d, Solver not set", __LINE__);

	//check that p1 is set
	int sid = inst->solver.id;
	if ((sid == GRASP2_NN || sid == GRASP3_NN || sid == EM_GRASP2 || sid == EM_GRASP3 || sid == VNS || sid == LOCAL_BRANCHING || sid == HARD_FIXING) && inst->solver.p1 < 0)
		print_error("%d, p1 not set", __LINE__);

	//check that p2 is set
	if ((sid == GRASP3_NN || sid == EM_GRASP3 || sid == HARD_FIXING) && inst->solver.p2 < 0)
		print_error("%d, p2 not set", __LINE__);

	//check that em_start is set
	if ((sid == MILEAGE || sid == EM_GRASP2 || sid == EM_GRASP3) && inst->solver.start < 0)
		print_error("%d, em_start not set", __LINE__);

	//check tenure is set
	if (sid == TABU && inst->tabu_tenure < 0)
		print_error("%d, tabu tenure not set", __LINE__);
}


void parse_TSPLIB(Instance* inst)
{
	FILE* fin = fopen(inst->inputfile, "r");
	printf("%s", inst->inputfile);
	if (fin == NULL)
		print_error("input file not opened correctly", __LINE__);

	char line[180];

	//read dimension and set nnnodes in instance
	while (fgets(line, sizeof(line), fin) != NULL)
	{
		char* par = strtok(line, " :");
		if (strncmp(par, "NODE_COORD_SECTION", 18) == 0)
		{
			if (inst->nnodes < 0)
				print_error("Starting points coordinate before having the dimension", __LINE__);
			break;
		}
		if (strcmp(par, "DIMENSION", 9) == 0 )
		{
			if (inst->nnodes > 0)
				print_error("Dimension of dataset already set", __LINE__);
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

void print_error(const char* msg, const int line)
{
	printf("\n===========================================================");
	printf("\nERROR at line %d: %s", line, msg);
	printf("\n===========================================================\n");
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
		print_error("Error in opening output data file", __LINE__);
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
		print_error("Error in opening output data file", __LINE__); 
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
		print_error("Error in opening output data file", __LINE__);
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
		print_error("Error in opening output data file", __LINE__);
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
	if (sid == MILEAGE || sid == EM_GRASP2 || sid == EM_GRASP3)
	{
		extra_mileage(inst);
	}

	else if (sid == NN || sid == GRASP2_NN || sid == GRASP3_NN || sid == ALL_START_NN)
	{
		nearest_neighbor(inst);
	}
	else if (sid == TABU)
	{
		tabu_search(inst);
	}
	else if (sid == VNS)
	{
		vns(inst);
	}
	else if (sid == TWO_OPT)
	{
		inst->solver.id = GRASP2_NN;
		grasp_iter(inst);
		// copy bestsol which has the solution for best start in currsol since two opt works on currsol
		copy_array(inst->bestsol, inst->currsol, inst->nnodes + 1);
		inst->currcost = inst->bestcost;
		if(inst->verbose>2)
			printf("\nStarting two optimization");
		two_opt(inst);
	}
	else if (sid == SIMANN)
	{
		simulated_annealing(inst);
	}
	else if (sid == GEN)
	{
		genetic(inst);
		return;
	}
	else if (sid == CPLEX_BENDERS || sid == CPLEX_CALLBACK)
	{
		TSPopt(inst);
	}
	else if (sid == LOCAL_BRANCHING)
	{
		//local_branching(inst);
	}
	else if (sid == HARD_FIXING)
	{
		//hard_fixing(inst);
	}
	else
		print_error("%d, Error in setting algorithm options for solving", __LINE__);


	//plot best instance found
	plot_generator(inst, inst->nnodes);
}

//plot cost of incumbent during solving (or cost in console but better to plot)
//tabulate formatt