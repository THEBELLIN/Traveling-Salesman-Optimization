#include "cplex_model.h"

int TSPopt(Instance* inst)
{

	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) 
	{
		print_error("%d, CPXopenCPLEX() error", __LINE__);
	}
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error)
	{
		print_error("%d, CPXcreateprob() error", __LINE__);
	}

	build_model(inst, env, lp);

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	if (inst->verbose > 10)
	{
		CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	}
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->randomseed);
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit);

	//TODO maybe add mipstart

	callback_solution(inst, env, lp);

	// use the optimal solution found by CPLEX
	/*
	
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			if (xstar[xpos(i, j, inst)] > 0.5)
			{
				printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1);
			}
		}
	}
	*/
	//TODO divide in separate functions with wrapper to choose the solver method
	//benders_loop(inst, env, lp);
	//just print one line with lb at every iteration

	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
	
	if (inst->verbose >= 0)
	{
		CPXwriteprob(env, lp, "../../Traveling-Salesman-Optimization/data/model.lp", NULL);
	}

	return 0; 
}

int xpos(int i, int j, Instance* inst)                                           
{
	if (i == j)
	{
		print_error("%d, i == j in xpos", __LINE__);
	}
	//i should be < than j
	if (i > j)
	{
		return xpos(j, i, inst);
	}
	int pos = i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;  //triangular matrix position
	return pos;
}

void build_model(Instance* inst, CPXENVptr env, CPXLPptr lp)
{
	double zero = 0;
	char binary = 'B';

	char** cname = CALLOC(1, char*);
	cname[0] = CALLOC(100, char);

	// add binary var.s x(i,j) for i < j  
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);  	
			double obj = distance(&inst->points[i], &inst->points[j]); 
			double lb = 0.0;
			double ub = 1.0;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
			{
				print_error("%d, wrong CPXnewcols on x var.s", __LINE__);
			}
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst))
			{
				print_error("%d, wrong position for x var.s", __LINE__);
			}
		}
	}

	// add the degree constraints 
	int* index = CALLOC(inst->nnodes, int);
	double* value = CALLOC(inst->nnodes, double);

	for (int h = 0; h < inst->nnodes; h++)  		// add the degree constraint on node h
	{
		double rhs = 2;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h + 1);
		int nnz = 0;
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h) 
				continue;
			index[nnz] = xpos(i, h, inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]))
		{
			print_error("%d, CPXaddrows(): error 1", __LINE__);
		}
	}

	free(value);
	free(index);
	free(cname[0]);
	free(cname);

}

// build succ() and comp() wrt xstar()...
void build_sol(const double* xstar, Instance* inst, int* succ, int* comp, int* ncomp) 
{
	*ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		succ[i] = -1;
		comp[i] = -1;
	}

	for (int start = 0; start < inst->nnodes; start++)
	{
		if (comp[start] >= 0) 
			continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while (!done)  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for (int j = 0; j < inst->nnodes; j++)
			{
				if (i != j && xstar[xpos(i, j, inst)] > 0.5 && comp[j] == -1) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}
		succ[i] = start;  // last arc to close the cycle
	}
}

double mip_value(CPXENVptr env, CPXLPptr lp)
{
	double zz;
	if (CPXgetobjval(env, lp, &zz))
		zz = CPX_INFBOUND;
	return zz;
}

void mip_delete_all_mipstarts(CPXENVptr env, CPXLPptr lp)
{
	int nmipstart = CPXgetnummipstarts(env, lp);
	if (nmipstart > 0 && CPXdelmipstarts(env, lp, 0, nmipstart - 1)) 
		print_error("%d, CPXdelmipstarts error", __LINE__);
}

int mip_solution_available(CPXENVptr env, CPXLPptr lp)
{
	double zz;
	if (CPXgetobjval(env, lp, &zz)) 
		return 0;
	return 1;
}

int mip_solved_to_optimality(CPXENVptr env, CPXLPptr lp)
{
	int lpstat = CPXgetstat(env, lp);
	int solved = (lpstat == CPXMIP_OPTIMAL) || (lpstat == CPXMIP_OPTIMAL_INFEAS) ||	(lpstat == CPXMIP_OPTIMAL_TOL);
	return solved;
}

int mip_infeasible(CPXENVptr env, CPXLPptr lp)
{
	int lpstat = CPXgetstat(env, lp);
	int infeas = (lpstat == CPXMIP_INFEASIBLE);
	return infeas;
}

void benders_loop(Instance* inst, CPXENVptr env, CPXLPptr lp)
{
	//time elapsed
	double time_elapsed = time(NULL) - inst->tstart; 

	//allocate memory for results
	int* succ = MALLOC(inst->nnodes, int);
	int* comp = MALLOC(inst->nnodes, int);
	int ncomp = 2;

	//Benders' loop
	while (time_elapsed < inst->time_limit && ncomp > 1)
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
			double patched_cost, prev_ub;
			patched_cost = patching(inst, ncomp, comp, succ); 

			//set CPLEX new cutoff (Upper Bound) to prune early edges that are worse than this value
			CPXgetdblparam(env, CPXPARAM_MIP_Tolerances_UpperCutoff, &prev_ub); 
			if (patched_cost < prev_ub)
			{
				CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_UpperCutoff, patched_cost); 
				printf("previous upperbound was: %f, new computed is %f\n", prev_ub, patched_cost); 
			}	  
		}
		time_elapsed = time(NULL) - inst->tstart;
		free(xstar);
	}

	//free
	free(comp);
	free(succ);
}

void add_SEC(Instance* inst, const int ncomp, int* comp, CPXENVptr env, CPXLPptr lp)
{
	// add the degree constraints 
	int ncols = CPXgetnumcols(env, lp);
	int* index = CALLOC(ncols, int);
	double* value = CALLOC(ncols, double);
	char** cname = CALLOC(1, char*);
	cname[0] = CALLOC(100, char);

	//for every component
	for (int c = 1; c <= ncomp; c++)  		
	{
		double rhs = 0;
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (comp[i] == c)
			{
				rhs++;
			}
		}
		//rhs to |S|-1
		rhs--;

		char sense = 'L';                            
		int nnz = 0;
		sprintf(cname[0], "SEC");	
		for (int i = 0; i < inst->nnodes - 1; i++)
		{
			if (comp[i] != c)
				continue;
			for (int j = i + 1; j < inst->nnodes; j++)
			{
				if (comp[j] != c)
					continue;
				index[nnz] = xpos(i, j, inst);
				value[nnz] = 1;
				nnz++;
			}
		}
		int izero = 0;
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]))
		{
			print_error("%d, CPXaddrows(): error 1", __LINE__);
		}
	}

	//free
	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}

//remember to allocate index, value, cname, cname[0] in advance
void get_SEC(Instance* inst, const int ncomp, int* comp, int* index, double* value, char** cname, int* rhs, char* sense, int* nnz, int c)
{
	// add the degree constraints 
	int ncols = inst->ncols;

	//for only component c
    *rhs = 0; 
    for (int i = 0; i < inst->nnodes; i++)
    {
        if (comp[i] == c)
        {
            (*rhs)++; 
        }
    }
    //rhs to |S|-1
    (*rhs)--; 

    *sense = 'L'; 
    *nnz = 0;
    sprintf(cname[0], "SEC");
    for (int i = 0; i < inst->nnodes - 1; i++)
    {
        if (comp[i] != c)
            continue;
        for (int j = i + 1; j < inst->nnodes; j++)
        {
            if (comp[j] != c)
                continue;
            index[*nnz] = xpos(i, j, inst);
            value[*nnz] = 1;
            (*nnz)++;
        }
    }
    int izero = 0;
	
	//free
	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}

double patching(Instance* inst, int ncomp, int* comp, int* succ)
{
	while (ncomp > 1)
	{
		int first_comp = 1;
		int second_comp = ncomp;
		int first_node, second_node, node_a, node_b;

		for (int i = 0; i < inst->nnodes; i++)
		{
			if (comp[i] == first_comp)
				first_node = i;
			else if (comp[i] == second_comp)
				second_node = i;
		}

		node_a = first_node;
		node_b = second_node;

		double min_delta_cost = INF_DOUBLE;
		int min_a, min_b;

		//scan every pair of edges (i,j) and (k,l)
		do
		{
			//consider the two edges (a, succ[a]) and (b, succ[b])
			edge e1, e2;
			e1.from = node_a;
			e1.to = succ[node_a];
			do
			{
				e2.from = node_b;
				e2.to = succ[node_b];

				//check the delta cost of swapping edge (a, succ[a]) with (a, succ[b]) and edge (b, succ[b]) with (b, succ[a])
				double delta_cost = COST(e1.from, e2.to) + COST(e2.from, e1.to) - COST(e1.from, e1.to) - COST(e2.from, e2.to);

				//if the delta cost is minimum, update the minimum delta cost and save the two edges' indices
				if (delta_cost < min_delta_cost)
				{
					min_delta_cost = delta_cost;
					min_a = e1.from;
					min_b = e2.from;
				}

				//update node b
				node_b = succ[node_b];

			} while (node_b != second_node);
			
			//update node a
			node_a = succ[node_a];

		} while (node_a != first_node);

		//swap the two edges
		int tmp = succ[min_a];
		succ[min_a] = succ[min_b];
		succ[min_b] = tmp;

		//update comp[i] for every i in the second component
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (comp[i] == second_comp)
				comp[i] = first_comp;
		}

		//update the number of connected components
		ncomp--;
	}

	//return the cost of the new patched solution
	double cost = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		cost += COST(i, succ[i]);
	}
	return cost;

	//TODO maybe add 2-opt before return to make better incumbent
}

void callback_solution(Instance* inst, CPXENVptr env, CPXLPptr lp)
{
	//get the number of columns
	int ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;

	//install the callback function
	if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION, my_callback, inst))
		print_error("%d, CPXcallbacksetfunc() error", __LINE__);

	//solve the problem with CPLEX
	int status = CPXmipopt(env, lp);
}

static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle)
{
	Instance* inst = (Instance*)userhandle;
	double* xstar = MALLOC(inst->ncols, double);
	double objval = CPX_INFBOUND;

	//check if we have integer solution -- LAZY CUT
	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
	{
		if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval))
			print_error("CPXcallbackgetcandidatepoint error");

		//allocate memory for results
		int* succ = MALLOC(inst->nnodes, int);
		int* comp = MALLOC(inst->nnodes, int);
		int ncomp;
		build_sol(xstar, inst, succ, comp, &ncomp);

		//check if the solution is feasible
		if (ncomp == 1)
			return;

		//allocate memory for results
		int* index = CALLOC(inst->ncols, int); 
		double* value = CALLOC(inst->ncols, double); 
		char** cname = CALLOC(1, char*); 
		cname[0] = CALLOC(100, char); 
		int* rhs;
		char* sense; 
		int* nnz;

		for (int c = 1; c <= ncomp; c++)
		{
			//maybe use add_SEC with a check for context in order to use addrow or reject sol
			get_SEC(inst, ncomp, comp, index, value, cname, rhs, sense, nnz, c);

			if (*nnz > 0) // means that the solution is infeasible and a violated cut has been found
			{
				int izero = 0;
				if (CPXcallbackrejectcandidate(context, 1, nnz, rhs, sense, &izero, index, value))
					print_error("CPXcallbackrejectcandidate() error"); // reject the solution and adds one cut 
			}

		}
		 
		double objheu = patching(inst, succ, comp, ncomp);

		int* ind = (int*)malloc(inst->ncols * sizeof(int));
		double* xheu = (double*)calloc(inst->ncols, sizeof(double));  // all zeros, initially
		succ_to_xheu(inst, succ, xheu);
		for (int j = 0; j < inst->ncols; j++) 
			ind[j] = j;
		//TODO maybe add only if it's better than incumbent
		if (CPXcallbackpostheursoln(context, inst->ncols, ind, xheu, objheu, CPXCALLBACKSOLUTION_NOCHECK)) 
			print_error("CPXcallbackpostheursoln() error");
		free(ind); 

		//free allocated memory
		free(xheu);
		free(index); 
		free(value); 
		free(cname[0]); 
		free(cname); 
		free(xstar);
	}

	//check if we have a fractional solution -- USER CUT
	if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
	{
		if (CPXcallbackgetrelaxationpoint(context, xstar, 0, inst->ncols - 1, &objval))
			print_error("CPXcallbackgetrelaxationpoint error");

		int local = 0;
		int purgeable = CPX_USECUT_FILTER;

		//TODO fractional cuts
		//CPXcallbackaddusercuts(context, 1, nnz, &rhs, &sense, &izero, cutind, cutval, &purgeable, &local)
	}

	return 0; 
} 

void succ_to_xheu(Instance* inst, int* succ, double* xheu)
{
	for (int i = 0; i < inst->nnodes; i++) 
		xheu[xpos(i, succ[i], inst)] = 1.0;
}