#ifndef H_CPLEX_MODEL
#define H_CPLEX_MODEL

#include "cplex.h"
#include "TSP.h"
#include "utility.h"
#include "heuristics.h"

//==================structs=================
typedef struct
{
	CPXCALLBACKCONTEXTptr context;
	Instance* inst;

}relaxation_callback_params;

//================functions=================
int TSPopt(Instance*);
int xpos(int, int, Instance*);
void build_model(Instance*, CPXENVptr, CPXLPptr);
void build_sol(const double*, Instance*, int*, int*, int*);
double mip_value(CPXENVptr, CPXLPptr);
void mip_delete_all_mipstarts(CPXENVptr, CPXLPptr);
int mip_solution_available(CPXENVptr, CPXLPptr);
int mip_solved_to_optimality(CPXENVptr, CPXLPptr);
int mip_infeasible(CPXENVptr, CPXLPptr);
void benders_loop(Instance*, CPXENVptr, CPXLPptr);
void add_SEC(Instance*, const int, int*, CPXENVptr, CPXLPptr);
void get_SEC(Instance*, const int, int*, int*, double*, char**, int*, char*, int*);
double patching(Instance*, int, int*, int*);
void callback_solution(Instance*, CPXENVptr, CPXLPptr);
static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr, CPXLONG, void*);
void succ_to_xheu(Instance*, int*, double*);
static int violated_cuts_callback(double, int, int*, void*);
static int add_SEC_cuts_fractional(Instance*, CPXCALLBACKCONTEXTptr, int, int*, int*, double*);
#endif
