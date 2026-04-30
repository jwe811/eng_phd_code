/* This program will count the total knottypes of a series of knot files	*/


#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>


/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDYIED **/

#define	L 2				/* number of horizontal edges  (x-direction)	*/
#define	M 1				/* number of vertical edges (y-direction)	*/
#define	totalspan 9			/* span of generated polygon (z-direction)	*/
#define	tot_num_files 1

#define	vM (M+1)				/* number of vertices in vertical direction   = M + 1      */
#define	vL (L+1)				/* number of vertices in horizontal direction = L + 1      */


/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/
main(void)
{
	system("date");  /* prints the date and time */
	clock();	/* Start clock to time program */

	int i, j, k;
	int unknot=0;
	int knot=0;
	int unknown=0;
	char line[8];
	char dummy[8];

	int totalfiles=tot_num_files;
	int filenum=1;
	char filename[60];	//this will hold the filename of the current text file of knottypes that will be counted


	for(filenum=1; filenum<=totalfiles; filenum++){	//for each knot file
		sprintf(filename, "ImpCreatorTwoStringsPrime_PKshrunkL%dM%dspan%dnum%d.txt_shrunkID", L, M, totalspan, filenum); //name of knot file
		FILE* fp = fopen(filename, "r");	//open "filename" (a knot file)

		if(fp != NULL){
			printf("file %s opened successfully\n", filename);
		}
		else{
			printf("file pointer is pointing to NULL\n");
			break;
		}

		while( fgets(line, 8, fp)!=NULL ){	//while there's a line to read		fgets(line, 7, fp);
			if( strcmp(line, "unknot\n")==0 ){	//we have an unknot
				unknot++;
			}
			else if ( strcmp(line, "knot\n")==0 ){	//we have a positive trefoil
				knot++;
			}
			else{	//unknown knot type
				printf("unknown found, it is: %s\n", line);
				scanf("%s", dummy);
				unknown++;
			}
		}



	}

	printf("unknots: %d\n", unknot);
	printf("knots: %d\n", knot);
	printf("unknown: %d\n", unknown);


	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	system("date");  /* prints the date and time */
}				/* end of main */
/***************************************************************************/


