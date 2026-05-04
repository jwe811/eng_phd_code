/* This program will count the total knottypes of a series of knot files	*/


#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>


/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/

#define	L 4
#define	M 1
#define	totalspan 150
#define	runnum 1

#define	maxfactors 15
//maximum number of factors in poly


/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/
int main(void)
{
	system("date");  /* prints the date and time */
	clock();	/* Start clock to time program */

	char line[100];
	char filename[100];

	sprintf(filename, "%dx%dhs%dr%d-id_master.txt", L, M, totalspan, runnum); //name of file
	FILE* fp = fopen(filename, "r");	//open "filename" (a knot file)

	if(fp != NULL){
		printf("file %s opened successfully\n", filename);
	}
	else{
		printf("file pointer is pointing to NULL. filename was %s. Exitting\n", filename);
		exit(1);
	}

	int i;
	int factors[maxfactors];
	for(i=0; i<maxfactors; i++){
		factors[i]=0;
	}
	int num_fact;//number of factors in line
	int blanklines=0;
	int bigger=0;

	while( fgets(line, 100, fp)!=NULL ){	//while there's a line to read
		i=0;
		num_fact=0;
		while( line[i] != '\n'){	//while not at the end of the line
			if(line[i]=='#' || line[i]=='c'){
				num_fact++;
			}
			i++;
		}
		if(i>0){
		  if(num_fact<maxfactors){
			factors[num_fact]++;
		  }
		  else{
		    bigger++;
		  }
		}
		else{
			blanklines++;
		}
	}

	for(i=0; i<maxfactors; i++){
		printf("%d factors: %d\n", i+1, factors[i]);
	}
	printf("There were %d polys with more than %d factors\n", bigger, maxfactors);
	printf("There were %d blank lines\n", blanklines);



	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	system("date");  /* prints the date and time */

	return 0;
}				/* end of main */
/***************************************************************************/


