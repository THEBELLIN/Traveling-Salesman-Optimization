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
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (inst->verbose > 10)
	{
		CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	}
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->randomseed);
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit);

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
	benders_loop(inst, env, lp);
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

	//allocate mamory for results
	int* succ = MALLOC(inst->nnodes, int);
	int* comp = MALLOC(inst->nnodes, int);

	//Benders' loop
	while (time_elapsed < inst->time_limit)
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
		int ncomp;
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
		sprintf(cname[0], "SEC(%d)", c);	//TODO fix progressive number
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