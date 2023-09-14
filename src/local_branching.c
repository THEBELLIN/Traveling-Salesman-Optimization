#include "local_branching.h"

#include "local_branching.h"

void local_branching(Instance* inst) {
	int k = 20;
	// define the fraction of time to spend for every call
	double time_per_call = inst->time_limit / 5.0;
	//double start = time(NULL);
	//printf("time start %f", start);
	inst->time_limit = 50;

	
	// initialize with an heuristic method x_0 first solution
	// a questo punto bestsol is initialized and has a best cost initialized

	vns(inst);
	double start = time(NULL);
	printf("\nheuristic solution cost %f", inst->bestcost);
	inst->time_limit = time_per_call;
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error)
	{
		print_error("%d, CPXopenCPLEX() error", __LINE__);
	}
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model hard_fixing");
	if (error)
	{
		print_error("%d, CPXcreateprob() error", __LINE__);
	}
	// fine apertura cplex
	if (inst->verbose > 10)
	{
		CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	}
	// build the degree model
	build_model(inst, env, lp);
	//----------------------
	//printf("fine building model\n");
	//set time limit
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit - (time(NULL) - start));
	//----------------------

	// allocate memory for local branching constraint
	int* index = CALLOC(inst->nnodes, int);
	double* value = CALLOC(inst->nnodes, double);
	char** cname = CALLOC(1, char*);
	cname[0] = CALLOC(100, char);
	sprintf(cname[0], "local branching");


	//int k = 20; // number of neighbourhood;

	inst->ncols = CPXgetnumcols(env, lp);
	//while cycle untill i have time
	int iter = 0;
	while (time(NULL) < start + inst->time_limit) {
		// post a solution
		add_mip_start(inst, env, lp, inst->bestsol);
		// add branching constraint
		double rhs = inst->nnodes - k;
		char sense = 'G';
		int nnz = 0;

		// find variables which are set to 1
		for (int i = 0; i < inst->nnodes; i++) {
			index[nnz] = xpos(inst->bestsol[i], inst->bestsol[i + 1], inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		int nrows = CPXgetnumrows(env, lp);// to get the index of the local branch constraint
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0])) print_error("Error in adding row for local branchin", __LINE__);

		//set the time to spend getting a better solution
		if (inst->time_limit - (time(NULL) - start) < time_per_call) {
			CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit - (time(NULL) - start));
		}
		else
			CPXsetdblparam(env, CPX_PARAM_TILIM, time_per_call);
		inst->tstart = time(NULL);
		//CALL CALLBACK FUNCTION TO FIND AN OPTIMIZED SOLUTION solition saved in bestsol if better than the previous solution found
		benders_loop2(inst, env, lp);
		//remove the local branching constraint added
		if (CPXdelrows(env, lp, nrows, nrows)) print_error("error in deleating a constraint", __LINE__);
		printf("Solution at time %f has cost %f\n", time(NULL) - start, inst->bestcost);
	}// while

	free(cname[0]);
	free(cname);
	free(index);
	free(value);

	// close enviroment 
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

}