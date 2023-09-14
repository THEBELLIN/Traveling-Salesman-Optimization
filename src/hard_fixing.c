#include "hard_fixing.h"

void hard_fixing(Instance* inst)
{
	double percentage = 0.6;
	double timelimit = inst->time_limit;
	// define the fraction of time to spend for every call
	double time_per_call = timelimit / 5.0;
	//double start = time(NULL);
	//printf("time start %f", start);
	inst->time_limit = 100;

	double start = time(NULL);
	// initialize with an heuristic method x_0 first solution
	// a questo punto bestsol is initialized and has a best cost initialized

	vns(inst);
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
	printf("fine building model\n");
	//set time limit
	CPXsetdblparam(env, CPX_PARAM_TILIM, timelimit - (time(NULL) - start));
	//----------------------
	inst->ncols = CPXgetnumcols(env, lp);
	//while cycle untill i have time
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
					print_error("error in changing the bound", __LINE__);
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
		callback_solution(inst, env, lp);
		//unfix all previous edges
		double zero = 0.0;
		for (int i = 0; i < inst->ncols; i++) {
			int pos = i;
			if (CPXchgbds(env, lp, 1, &pos, &lb, &zero))
				print_error("error in changing the bound", __LINE__);
		}
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
		print_error("CPXaddmipstarts() error", __LINE__);
	free(ind);
	free(xheu);

}

void benders_loop2(Instance* inst, CPXENVptr env, CPXLPptr lp)
{
	//time elapsed
	double time_elapsed = time(NULL) - inst->tstart;
	//allocate mamory for results
	int* succ = MALLOC(inst->nnodes, int);
	int* comp = MALLOC(inst->nnodes, int);
	int ncomp = 2;
	//Benders' loop
	while (time_elapsed < inst->time_limit && ncomp>1)
	{
		//set time limit
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit - time_elapsed);

		//solve the problem
		int error = CPXmipopt(env, lp);
		if (error)
		{
			printf("CPX error code %d\n", error);
			print_error("%d, CPXmipopt() error", __LINE__);
		}

		//get results
		int ncols = CPXgetnumcols(env, lp);

		double* xstar = CALLOC(ncols, double);
		if (CPXgetx(env, lp, xstar, 0, ncols - 1))
		{
			print_error("%d, CPXgetx() error", __LINE__);
		}
		build_sol(xstar, inst, succ, comp, &ncomp);

		//check number of connected components
		if (ncomp > 1)
		{
			add_SEC(inst, ncomp, comp, env, lp);
		}
		time_elapsed = time(NULL) - inst->tstart;
		free(xstar);
	}
	if (ncomp == 1 && inst->bestcost > mip_value(env, lp)) {
		inst->bestcost = mip_value(env, lp);
		transform_in_perm_and_save(succ, inst);
		printf("Updated ");
	}
	//free
	free(comp);
	free(succ);
}