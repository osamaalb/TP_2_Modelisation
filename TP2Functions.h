#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include<ilcplex/cplex.h>
#include <assert.h>




//IP Mproblem structure
typedef struct IP_problem
{
        //Internal structures
        CPXENVptr env;
        CPXLPptr lp;
	//number of variables
	int nv;
	//number of constraints
	int nc;

        //Output solution
        double *x;
        //Costs array
        double *cost;
        //TYpe of the variables (binary or continuous)
        char *c_type;
        //Bounds over the variables
        double *up_bound;
        double *low_bound;
	//Names of the variables
	char** var_name;

        //Right hand section of the constraints
        double *rhs;
        //Sense of the constraints
        char *sense;
        //Left hand section for data constraints
        int *rmatbeg;
        int *rmatind;
        double *rmatval;
	int nz;
	//Names of the constraints 
	char** const_name;
        //Solution status
        int solstat;
        double objval;
} IP_problem;

typedef struct dataSet 
{
	//Attributes of the instance
	//Nombre d'objets
	int n;
	//Capacite b
	int b;
	//Capacite g 
	int g;

	//Tableau d'entiers de taille n contenant la valeur de chacun des objets
	int*c;
	//Tableau d'entiers de taille n contenant le poids de chacun des objets
	int*a;
	//Tableau d'entiers de taille n contenant le poids #2 de chacun des objets
	int*f;

	//CPLEX problem data
	IP_problem master;

} dataSet;


int read_TP2_instance(FILE*fin,dataSet* dsptr);
int generate_TP2_instance(dataSet* dsptr, int n, int b, int g, int max_c, int max_a, int max_f);
int TP2_solve_exact(dataSet* dsptr);
int TP2_solve_exact_1D(dataSet* dsptr);
int TP2_relax_lagr(dataSet* d, float tolerance);


