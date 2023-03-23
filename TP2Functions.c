#include "TP2Functions.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <stdio.h>
#include <ilcplex/cplex.h>

int read_TP2_instance(FILE*fin,dataSet* dsptr)
{
	int rval = 0;

	//capacites b et g
	int b;
	int g;
	//Nombre d'objets
	int n;
	rval = fscanf(fin,"%d,%d,%d\n",&n,&b,&g);
	dsptr->b = b;
	dsptr->g = g;

	dsptr->n = n;
	dsptr->c = (int*)malloc(sizeof(int)*n);
	dsptr->a = (int*)malloc(sizeof(int)*n);
	dsptr->f = (int*)malloc(sizeof(int)*n);


	int i;
	for( i = 0 ; i < n ; i++)
		rval = fscanf(fin,"%d,%d,%d\n",&(dsptr->c[i]),&(dsptr->a[i]),&(dsptr->f[i]));

	fprintf(stderr,"\nInstance file read, we have capacites [b,g] = [%d,%d] and there is %d items of values/weights:\n",
			b,g,n);
	fprintf(stderr,"i\tci\tai\tfi\n");
	fprintf(stderr,"--------------------------\n");


	for( i = 0 ; i < n ; i++)
		fprintf(stderr,"%d\t%d\t%d\t%d\n",i,dsptr->c[i],dsptr->a[i],dsptr->f[i]);
	fprintf(stderr,"\n");

	return rval;
}


// une fonction pour générer un problème aléatoirement
// n est le nombre des objets
// b est le volume maximum
// g est le volume #2 maximum
// max_c est la valeur maximum de la valeur des objets
// max_a est la valeur maximum du poids des objets
// max_f est la valeur maximum du poids #2 des objets
int generate_TP2_instance(dataSet* dsptr, int n, int b, int g, int max_c, int max_a, int max_f) {
	// initialiser les variables
	dsptr->n = n;
	dsptr->b = b;
	dsptr->g = g;
	dsptr->c = (int*)malloc(sizeof(int)*n);
	dsptr->a = (int*)malloc(sizeof(int)*n);
	dsptr->f = (int*)malloc(sizeof(int)*n);

	// détérminer le seed
	srand(time(NULL));

	// générer n objets alétoirement
	for (int i = 0; i < n; i++) {
		dsptr->c[i] = rand() % max_c + 1;
		dsptr->a[i] = rand() % max_a + 1;
		dsptr->f[i] = rand() % max_f + 1;
	}

	return 0;
}


int TP2_solve_exact(dataSet* dsptr)
{
	int rval = 0;

	IP_problem* ip_prob_ptr = &(dsptr->master);
	ip_prob_ptr->env = NULL;
	ip_prob_ptr->lp = NULL;
	ip_prob_ptr->env = CPXopenCPLEX (&rval);
	if(rval) fprintf(stderr,"ERROR WHILE CALLING CPXopenCPLEX\n");
	if ( ip_prob_ptr->env == NULL ) 
	{
		char  errmsg[1024];
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (ip_prob_ptr->env, rval, errmsg);
		fprintf (stderr, "%s", errmsg);
		exit(0);	
	}

	//We create the MIP problem
	ip_prob_ptr->lp = CPXcreateprob (ip_prob_ptr->env, &rval, "TP2");
	if(rval) fprintf(stderr,"ERROR WHILE CALLING CPXcreateprob\n");

	//verify the data passed as argument and shows the intermediate operations while doing the optimisation
	rval = CPXsetintparam (ip_prob_ptr->env, CPX_PARAM_DATACHECK, CPX_ON); 
	rval = CPXsetintparam (ip_prob_ptr->env, CPX_PARAM_SCRIND, CPX_ON);

	int n = dsptr->n;
	int nv = n;
	int* a = dsptr->a;
	int* f = dsptr->f;
	int* c = dsptr->c;
	int b = dsptr->b;
	int g = dsptr->g;

	//We fill our arrays
	//Memory
	ip_prob_ptr->nv = nv;
	// where to store the solution
	ip_prob_ptr->x = (double*)malloc(sizeof(double)*nv);
	// the costs
	ip_prob_ptr->cost = (double*)malloc(sizeof(double)*nv);
	// the type of each variable: binary, integer, continuous, etc.
	ip_prob_ptr->c_type = (char*)malloc(sizeof(char)*nv);
	// bounds of each variable value
	ip_prob_ptr->up_bound = (double*)malloc(sizeof(double)*nv);
	ip_prob_ptr->low_bound = (double*)malloc(sizeof(double)*nv);
	// the name of each variable
	ip_prob_ptr->var_name = (char**)malloc(sizeof(char*)*nv);

	int j,id = 0;
	//Structures keeping the index of each variable
	int*id_x_j = (int*)malloc(sizeof(int)*n);

	//variables xi definition
	for( j = 0 ; j < n ; j++)
	{
		//We keep the id
		id_x_j[j] = id;

		//We generate the variable attributes
		ip_prob_ptr->x[id] = 0;
		ip_prob_ptr->cost[id] = c[j];
		ip_prob_ptr->c_type[id] = 'B'; // B: binary, C: continuous, I: integer
		ip_prob_ptr->up_bound[id] = 1;
		ip_prob_ptr->low_bound[id] = 0;
		ip_prob_ptr->var_name[id] = (char*)malloc(sizeof(char)*1024);
	        snprintf(       ip_prob_ptr->var_name[id],
        	                1024,
                	        "x_j%d",
                        	j+1);
		id++;
	}


	rval = CPXnewcols( ip_prob_ptr->env, ip_prob_ptr->lp, 
			nv, 
			ip_prob_ptr->cost, 
			ip_prob_ptr->low_bound,
			ip_prob_ptr->up_bound,
			ip_prob_ptr->c_type,
			ip_prob_ptr->var_name);
	if(rval)
		fprintf(stderr,"CPXnewcols returned errcode %d\n",rval);



	//Constraints part
	ip_prob_ptr->rhs = (double*)malloc(sizeof(double));
	ip_prob_ptr->sense = (char*)malloc(sizeof(char));
	ip_prob_ptr->rmatbeg = (int*)malloc(sizeof(int));
	ip_prob_ptr->nz = n;


	ip_prob_ptr->rmatind = (int*)malloc(sizeof(int)*nv);
	ip_prob_ptr->rmatval = (double*)malloc(sizeof(double)*nv);
	ip_prob_ptr->const_name = (char**)malloc(sizeof(char*));
	ip_prob_ptr->const_name[0] = (char*)malloc(sizeof(char)*1024);

	//We fill what we can 
	ip_prob_ptr->rmatbeg[0] = 0;

	//We generate and add each constraint to the model

	//capacity constraint #1
	ip_prob_ptr->rhs[0] = b;
	ip_prob_ptr->sense[0] = 'L'; // L: less than or equal, G greater than or equal, E equal to
	//Constraint name
	snprintf(       ip_prob_ptr->const_name[0],
			1024,
			"capacity_1"
			);
	id=0;
	//variables x_j coefficients
	for( j = 0 ; j < n ; j++)
	{
		ip_prob_ptr->rmatind[id] = id_x_j[j];
		ip_prob_ptr->rmatval[id] =  a[j];
		id++;
	}
	rval = CPXaddrows( ip_prob_ptr->env, ip_prob_ptr->lp, 
			0,//No new column
			1,//One new row
			n,//Number of nonzero coefficients
			ip_prob_ptr->rhs, 
			ip_prob_ptr->sense, 
			ip_prob_ptr->rmatbeg, 
			ip_prob_ptr->rmatind, 
			ip_prob_ptr->rmatval,
			NULL,//No new column
			ip_prob_ptr->const_name );
	if(rval)
		fprintf(stderr,"CPXaddrows returned errcode %d\n",rval);

	//capacity constraint 2#
	ip_prob_ptr->rhs[0] = g;
	ip_prob_ptr->sense[0] = 'L';
	//Constraint name
	snprintf(       ip_prob_ptr->const_name[0],
			1024,
			"capacity_2"
			);
	id=0;
	//variables x_j coefficients
	for( j = 0 ; j < n ; j++)
	{
		ip_prob_ptr->rmatind[id] = id_x_j[j];
		ip_prob_ptr->rmatval[id] =  f[j];
		id++;
	}
	rval = CPXaddrows( ip_prob_ptr->env, ip_prob_ptr->lp, 
			0,//No new column
			1,//One new row
			n,//Number of nonzero coefficients
			ip_prob_ptr->rhs, 
			ip_prob_ptr->sense, 
			ip_prob_ptr->rmatbeg, 
			ip_prob_ptr->rmatind, 
			ip_prob_ptr->rmatval,
			NULL,//No new column
			ip_prob_ptr->const_name );
	if(rval)
		fprintf(stderr,"CPXaddrows returned errcode %d\n",rval);

	//We switch to maximization
	rval = CPXchgobjsen( ip_prob_ptr->env, ip_prob_ptr->lp, CPX_MAX );


	//We write the problem for debugging purposes, can be commented afterwards
	rval = CPXwriteprob (ip_prob_ptr->env, ip_prob_ptr->lp, "multiDKnapsack.lp", NULL);
	if(rval)
		fprintf(stderr,"CPXwriteprob returned errcode %d\n",rval);

	//We solve the model
	rval = CPXmipopt (ip_prob_ptr->env, ip_prob_ptr->lp); // solve the problem with integers ; CPXlpopt solves using continuous numbers. If the variables are declared as integers, this will not work
	if(rval)
		fprintf(stderr,"CPXmipopt returned errcode %d\n",rval);

	rval = CPXsolwrite( ip_prob_ptr->env, ip_prob_ptr->lp, "multiDKnapsack.sol" );
	if(rval)
		fprintf(stderr,"CPXsolwrite returned errcode %d\n",rval);

	//We get the objective value
	rval = CPXgetobjval( ip_prob_ptr->env, ip_prob_ptr->lp, &(ip_prob_ptr->objval) );
	if(rval)
		fprintf(stderr,"CPXgetobjval returned errcode %d\n",rval);

	//We get the best solution found 
	rval = CPXgetobjval( ip_prob_ptr->env, ip_prob_ptr->lp, &(ip_prob_ptr->objval) );
	rval = CPXgetx( ip_prob_ptr->env, ip_prob_ptr->lp, ip_prob_ptr->x, 0, nv-1 );
	if(rval)
		fprintf(stderr,"CPXgetx returned errcode %d\n",rval);

	//We display the solution
	double tolerance = 0.0001;

	/**************** FILL HERE ***************/



	return rval;
}


int TP2_solve_exact_1D(dataSet* dsptr)
{
	int rval = 0;

	IP_problem* ip_prob_ptr = &(dsptr->master);
	ip_prob_ptr->env = NULL;
	ip_prob_ptr->lp = NULL;
	ip_prob_ptr->env = CPXopenCPLEX (&rval);
	if(rval) fprintf(stderr,"ERROR WHILE CALLING CPXopenCPLEX\n");
	if ( ip_prob_ptr->env == NULL ) 
	{
		char  errmsg[1024];
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (ip_prob_ptr->env, rval, errmsg);
		fprintf (stderr, "%s", errmsg);
		exit(0);	
	}

	//We create the MIP problem
	ip_prob_ptr->lp = CPXcreateprob (ip_prob_ptr->env, &rval, "TP2");
	if(rval) fprintf(stderr,"ERROR WHILE CALLING CPXcreateprob\n");

	//verify the data passed as argument and shows the intermediate operations while doing the optimisation
	rval = CPXsetintparam (ip_prob_ptr->env, CPX_PARAM_DATACHECK, CPX_ON); 
	rval = CPXsetintparam (ip_prob_ptr->env, CPX_PARAM_SCRIND, CPX_ON);

	int n = dsptr->n;
	int nv = n;
	int* a = dsptr->a;
	int* c = dsptr->c;
	int b = dsptr->b;

	//We fill our arrays
	//Memory
	ip_prob_ptr->nv = nv;
	// where to store the solution
	ip_prob_ptr->x = (double*)malloc(sizeof(double)*nv);
	// the costs
	ip_prob_ptr->cost = (double*)malloc(sizeof(double)*nv);
	// the type of each variable: binary, integer, continuous, etc.
	ip_prob_ptr->c_type = (char*)malloc(sizeof(char)*nv);
	// bounds of each variable value
	ip_prob_ptr->up_bound = (double*)malloc(sizeof(double)*nv);
	ip_prob_ptr->low_bound = (double*)malloc(sizeof(double)*nv);
	// the name of each variable
	ip_prob_ptr->var_name = (char**)malloc(sizeof(char*)*nv);

	int j,id = 0;
	//Structures keeping the index of each variable
	int*id_x_j = (int*)malloc(sizeof(int)*n);

	//variables xi definition
	for( j = 0 ; j < n ; j++)
	{
		//We keep the id
		id_x_j[j] = id;

		//We generate the variable attributes
		ip_prob_ptr->x[id] = 0;
		ip_prob_ptr->cost[id] = c[j];
		ip_prob_ptr->c_type[id] = 'B'; // B: binary, C: continuous, I: integer
		ip_prob_ptr->up_bound[id] = 1;
		ip_prob_ptr->low_bound[id] = 0;
		ip_prob_ptr->var_name[id] = (char*)malloc(sizeof(char)*1024);
	        snprintf(       ip_prob_ptr->var_name[id],
        	                1024,
                	        "x_j%d",
                        	j+1);
		id++;
	}


	rval = CPXnewcols( ip_prob_ptr->env, ip_prob_ptr->lp, 
			nv, 
			ip_prob_ptr->cost, 
			ip_prob_ptr->low_bound,
			ip_prob_ptr->up_bound,
			ip_prob_ptr->c_type,
			ip_prob_ptr->var_name);
	if(rval)
		fprintf(stderr,"CPXnewcols returned errcode %d\n",rval);



	//Constraints part
	ip_prob_ptr->rhs = (double*)malloc(sizeof(double));
	ip_prob_ptr->sense = (char*)malloc(sizeof(char));
	ip_prob_ptr->rmatbeg = (int*)malloc(sizeof(int));
	ip_prob_ptr->nz = n;


	ip_prob_ptr->rmatind = (int*)malloc(sizeof(int)*nv);
	ip_prob_ptr->rmatval = (double*)malloc(sizeof(double)*nv);
	ip_prob_ptr->const_name = (char**)malloc(sizeof(char*));
	ip_prob_ptr->const_name[0] = (char*)malloc(sizeof(char)*1024);

	//We fill what we can 
	ip_prob_ptr->rmatbeg[0] = 0;

	//We generate and add each constraint to the model

	//capacity constraint #1
	ip_prob_ptr->rhs[0] = b;
	ip_prob_ptr->sense[0] = 'L'; // L: less than or equal, G greater than or equal, E equal to
	//Constraint name
	snprintf(       ip_prob_ptr->const_name[0],
			1024,
			"capacity_1"
			);
	id=0;
	//variables x_j coefficients
	for( j = 0 ; j < n ; j++)
	{
		ip_prob_ptr->rmatind[id] = id_x_j[j];
		ip_prob_ptr->rmatval[id] =  a[j];
		id++;
	}
	rval = CPXaddrows( ip_prob_ptr->env, ip_prob_ptr->lp, 
			0,//No new column
			1,//One new row
			n,//Number of nonzero coefficients
			ip_prob_ptr->rhs, 
			ip_prob_ptr->sense, 
			ip_prob_ptr->rmatbeg, 
			ip_prob_ptr->rmatind, 
			ip_prob_ptr->rmatval,
			NULL,//No new column
			ip_prob_ptr->const_name );
	if(rval)
		fprintf(stderr,"CPXaddrows returned errcode %d\n",rval);

	//We switch to maximization
	rval = CPXchgobjsen( ip_prob_ptr->env, ip_prob_ptr->lp, CPX_MAX );


	//We write the problem for debugging purposes, can be commented afterwards
	rval = CPXwriteprob (ip_prob_ptr->env, ip_prob_ptr->lp, "1DKnapsack.lp", NULL);
	if(rval)
		fprintf(stderr,"CPXwriteprob returned errcode %d\n",rval);

	//We solve the model
	rval = CPXmipopt (ip_prob_ptr->env, ip_prob_ptr->lp); // solve the problem with integers ; CPXlpopt solves using continuous numbers. If the variables are declared as integers, this will not work
	if(rval)
		fprintf(stderr,"CPXmipopt returned errcode %d\n",rval);

	rval = CPXsolwrite( ip_prob_ptr->env, ip_prob_ptr->lp, "1DKnapsack.sol" );
	if(rval)
		fprintf(stderr,"CPXsolwrite returned errcode %d\n",rval);

	//We get the objective value
	rval = CPXgetobjval( ip_prob_ptr->env, ip_prob_ptr->lp, &(ip_prob_ptr->objval) );
	if(rval)
		fprintf(stderr,"CPXgetobjval returned errcode %d\n",rval);

	//We get the best solution found 
	rval = CPXgetobjval( ip_prob_ptr->env, ip_prob_ptr->lp, &(ip_prob_ptr->objval) );
	rval = CPXgetx( ip_prob_ptr->env, ip_prob_ptr->lp, ip_prob_ptr->x, 0, nv-1 );
	if(rval)
		fprintf(stderr,"CPXgetx returned errcode %d\n",rval);

	//We display the solution
	double tolerance = 0.0001;

	return rval;
}

dataSet TP2_solve_Pk(dataSet* d, double l) {
	// copier le dataset avec un seul modification c[j] devient c[j] - (l * f[j])
	dataSet d_k;
	d_k.n = d->n;
	d_k.b = d->b;
	d_k.c = (int*)malloc(sizeof(int)*d->n);
	d_k.a = (int*)malloc(sizeof(int)*d->n);

	for (int i = 0; i < d->n; i++) {
		d_k.c[i] = d->c[i] - (l * d->f[i]);
		d_k.a[i] = d->a[i];
	}

	TP2_solve_exact_1D(&d_k);

	d_k.master.objval += (l * d->g);

	return d_k;
}

double max(double a, double b) {
	if (a > b) {
		return a;
	}
	return b;
}

int TP2_relax_lagr(dataSet* d, float tolerance) {
	double* x_ = (double*)malloc(sizeof(double)*d->n);
	for(int i = 0; i < d->n; i++) {
		x_[i] = 0;
	}

	double z_ = INT_MAX;

	double z_x_ = 0;

	int k = 0;

	// pas
	double p = 1; // on a choisit la mis à jour de pas : p = p / 4

	// lambda
	double l = 0;

	double z_k = 0;

	while (p > tolerance || (z_ - z_k) > tolerance) {
		// résoudre le problem Pk
		dataSet d_k = TP2_solve_Pk(d, l);
		z_k = d_k.master.objval;
		
		// le volume utilisé par la solution de Pk pour la contarainte 2
		int used_vol_2 = 0;
		for (int i = 0; i < d_k.n; i++) {
			if (d_k.master.x[i] == 1) {
				used_vol_2 += d->f[i];
			}
		}

		if (used_vol_2 <= d->g && z_k > z_x_) {
			for (int i = 0; i < d->n; i++) {
				x_[i] = d_k.master.x[i];
			}
			z_x_ = z_k;
		}

		if (z_k < z_) {
			z_ = z_k;
		}

		double gamma = d->g - used_vol_2;

		l = max(0, l - p*gamma);

		k++;
		p = p / 4;
	}

	fprintf(stderr,"z_ relax. lagrangienne = %f\n", z_);

	return 0;
}