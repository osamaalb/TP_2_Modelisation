#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "TP2Functions.h"
#include<ilcplex/cplex.h>


int main(int argc, char **argv)
{
	int rval =0;	
	//File instance name
	//-F option
        char instance_file[1024];
        snprintf(       instance_file,
                        1024,
                        "%s",
                        "instance.csv");

	char c;
        while ((c=getopt (argc, argv,"F:h")) != EOF)
	{
		switch(c)
		{
                        case 'F':
				snprintf(       instance_file,
						1024,
						"%s",
						optarg);
				break;

			case 'h':
				fprintf(stderr,"Usage: ./TP2 [options]\nOptions:\n\n");
				fprintf(stderr,"******** INSTANCE DATA ********\n");
				fprintf(stderr,"\t-F Instance file name to load............................(default %s).\n",instance_file);
				break;
			default:
				exit(0);
		}
	}

	dataSet data;

	//Open the instance file
	FILE* fin = fopen(instance_file,"r");
	read_TP2_instance(fin,&data);
	fclose(fin);

	//Generate a random instance
	generate_TP2_instance(&data, 10, 50, 40, 15, 15, 12);

	fprintf(stderr,"%d,%d,%d\n",data.n,data.b,data.g);
	fprintf(stderr,"------------------------\n");
	for( int i = 0 ; i < data.n ; i++)
		fprintf(stderr,"%d,%d,%d\n",data.c[i],data.a[i],data.f[i]);

	fprintf(stderr,"\n\n");
	//execute your solution methods on the instance you just read
	//Exact solution
	TP2_solve_exact(&data);

	//Exact solution 1D
	TP2_solve_exact_1D(&data);

	return rval;
}

