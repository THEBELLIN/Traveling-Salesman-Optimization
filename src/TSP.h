#ifndef VRP_H_ //avoid defining twice if included twice

#define VRP_H_

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//==================enums==================
typedef enum
{
	NN,
	NN_GRASP2,
	NN_GRASP3,
	EM,			//always select the best option at each step -> deterministic	
	EM_GRASP2,	//use the second best option with probability 1-p
	EM_GRASP3,	//use the 1st option with probability p1, the 2nd with p2, 3rd with 1-p1-p2
	GEN, 
	VNS, 
	SIMANN, 
	TABU,
	CPLEX_BENDERS, 
	CPLEX_CALLBACK,
	LOCAL_BRANCHING,
	HARD_FIXING
}solver_id;

//option to select the method for extra mileage algorithm  start
typedef enum
{
	RAND,		//select 2 nodes at random for the strating set
	MAX_DIST,	//select the 2 most distant nodes for the starting set
	CONV_HULL	//start with the convex hull
}em_start;

//===============structs===============
typedef struct
{
	double x;
	double y;
}Point;

typedef struct
{
	int from;
	int to;
}myedge;

typedef struct
{
	solver_id id;
	double p1, p2; //probabilities for GRASP
	em_start start;
}solver;

typedef struct
{
	Point* points;
	int nnodes;
	double* cost;
	solver solver;
	int ncols; //ncols of CPLEX model

	char inputfile[1000];
	char outputfile[1000];
	int randomseed;

	double tstart;
	double time_limit;
	double bestcost, tbest, bestlb;		// incumbent solution
	int* bestsol;						// incumbent solution
	double currcost;					// current solution	
	int* currsol;						// current solution
	int verbose;						//  0, 10, 20, 30, 40, 50

	int* tabu;
	int tabu_tenure;

}Instance;

//=====================functions=======================================
void parse_TSPLIB(Instance*);
void print_error(const char*);
void initialize_instance(Instance*);
bool check_feasibility(int*, int);
void parse_args(Instance*, int, char**);
void print_points(Instance*);
Point* generate_random_points(int);
Point* generate_random_points_range(int, double, double);
void free_instance(Instance*);
double distance(Point*, Point*);
//void plot_generator(Instance*);
void plot_generator(Instance* inst, int n_edges);
void print_points_file(Point*, int, FILE*);
Instance* generate_test_bed(int, int, int);
double get_cost(Instance*, int*);
void save_if_best(Instance*);
void file_perf_prof(int, int, int);
void sol_to_tsp(int*, double*, Instance*);
void transform_in_perm_and_save(int*, Instance*);
void solve(Instance*); 

#endif
