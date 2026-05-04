/*This program will take in a UofS polygons file, and will count how many of the polygons in the file(s) have a certain number
of edges. This is being used to find the number of compact polygons in a group of polygons that are in a tube */

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>


main(void)
{
	printf("\nmain starting\n\n");

	const int L = 2;	//CHANGEABLE
	const int M = 1;	//CHANGEABLE
	const int span = 4;	//CHANGEABLE
	const int maxfiles = 11;

	int filenum=1;

	int i, j, k;
	int NumEdgeInterest = (L+1)*(M+1)*(span+1);	//Compact check

	int edgecount=0;	//will hold current number of edges in curpoly
	int curpoly=0;		//will count number of polys in file
	int InterestCount=0;	//will count number of polys in the file with NumEdgeInterest edges

	char line[7];		//line that is read in from file
	int x, y, z;
	int direction=0;

	char filename[45];

	for(filenum=1; filenum<=maxfiles; filenum++){
		sprintf(filename, "CreatorAllPolysL%dM%dspan%dnum%d.txt", L, M, span,filenum);	//CHANGEABLE
		FILE* fp = fopen(filename, "r");
		if(fp != NULL){
			printf("Opened '%s'\n", filename);
		}
		else{
			printf("Couldn't open '%s', fp is pointing at NULL. Stopping here!\n", filename);
			break;
		}

		//read in the next polys
		while( fgets(line, 7, fp)!=NULL ){
			sscanf(line, "%d %d %d", &x, &y, &z);
			if(x==-999){
				break;
			}
			curpoly++;

			while(1){
				fgets(line, 7, fp);
				sscanf(line, "%d", &direction);
				if(direction==-111){
					break;
				}
				else{
					edgecount++;
				}
			}
			//polygon is read

			if(edgecount==NumEdgeInterest){	//if polygon had the number of edges of interest
				InterestCount++;			//count it
			}
			edgecount=0;	//reset for next polygon.
		}
		fclose(fp);
		x=0;
	}

	printf("Number of polygons in all files: %d\n",curpoly);
	printf("Number of polygons with %d edges is %d\n",NumEdgeInterest,InterestCount);


}











































