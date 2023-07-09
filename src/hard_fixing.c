#include "hard_fixing.h"

void hard_fixing(Instance* inst, double percentage, double timelimit) {
	// define the fraction of time to spend for every call
	double time_per_call = timelimit / 10.0;
	double start = time(NULL);
	printf("time start %f", start);
	inst->time_limit = time_per_call;
	printf("time limit %f ", inst->time_limit);

	// initialize with an heuristic method x_0 first solution
	// a questo punto bestsol is initialized and has a best cost initialized
	nearest_neighbor_grasp_random(inst, 0, 0.5);
	two_opt(inst, inst->bestsol);
	inst->bestcost = calculateCost(inst, inst->bestsol);
	printf("heuristic solution cost %f", inst->bestcost);

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
	printf("fine building model\n");
	//set time limit
	CPXsetdblparam(env, CPX_PARAM_TILIM, timelimit - (time(NULL) - start));
	//----------------------

	//while cycle untill i have time
	int iter = 0;
	while (time(NULL) < start + timelimit) {
		//add the best soultion found as a start to cplex
		add_mip_start(inst, env, lp, inst->bestsol);

		//fix random edges to 1 for hard fixing
		double one = 1.0;
		char lb = 'L';
		for (int i = 0; i < inst->nnodes; i++) {
			if (rand01() < percentage) {
				int pos = xpos(inst->bestsol[i], inst->bestsol[i + 1], inst);
				if (CPXchgbds(env, lp, 1, &pos, &lb, &one))
					print_error("error in changing the bound");
			}
		}

		//set the time to spend getting a better solution
		if (timelimit - (time(NULL) - start) < time_per_call) {
			CPXsetdblparam(env, CPX_PARAM_TILIM, timelimit - (time(NULL) - start));
		}
		else
			CPXsetdblparam(env, CPX_PARAM_TILIM, time_per_call);
		inst->tstart = time(NULL);
		//CALL CALLBACK FUNCTION TO FIND AN OPTIMIZED SOLUTION solition saved in bestsol if better than the previous solution found
		benders_loop(inst, env, lp);
		//unfix all previous edges
		double zero = 0.0;
		for (int i = 0; i < inst->ncols; i++) {
			int pos = i;
			if (CPXchgbds(env, lp, 1, &pos, &lb, &zero))
				print_error("error in changing the bound");
		}
		iter++;
		printf("Solution at time %f has cost %f\n", time(NULL) - start, inst->bestcost);
	}// while
	// close enviroment 
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

}


// function that adds a solution in the node format as a start in cplex
void add_mip_start(Instance* inst, CPXENVptr env, CPXLPptr lp) {
	double* xheu = (double*)calloc(inst->ncols, sizeof(double));
	sol_to_tsp(inst->bestsol, xheu, inst);
	int* ind = (int*)malloc(inst->ncols * sizeof(int));
	for (int i = 0; i < inst->ncols; i++) {
		ind[i] = i;// values that i want to set in this case all of them
	}
	int effortlevel = CPX_MIPSTART_NOCHECK;
	int beg = 0;
	if (CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, ind, xheu, &effortlevel, NULL))
		print_error("CPXaddmipstarts() error");
	free(ind);
	free(xheu);

}