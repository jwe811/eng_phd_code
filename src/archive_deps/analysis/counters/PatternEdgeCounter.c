//Takes a polygon file (or set of files) and for each polygon, counts how many edges are in the "pattern" contained in each poly.

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>



main(void)
{
	printf("\nmain starting\n\n");

	const int L = 3;	//CHANGEABLE
	const int M = 1;	//CHANGEABLE
	const int span = 5;	//CHANGEABLE
	const int totfiles = 12; //CHANGEABLE

	int i, j, k;
	int curfile;

	char FileName[45];	//files with the polygons (patterns) of interest
	char line[7];		//line in file

	int x, y, z;
	int polygon[4+(L+1)*(M+1)*(span+1)];	//first 3 entries (0,1,2) will be the starting xyz
	int curpoly=0;
	int direction=0;
	int curindex;
	int counts[(L+1)*(M+1)*(span+1)];
	int length;

	for(i=0; i<=(L+1)*(M+1)*(span+1)-1; i++){
		counts[i]=0;
	}

	for(curfile=1; curfile<=totfiles; curfile++){
		sprintf(FileName, "Global_InterestPKshrunkL%dM%dspan%dnum%d_v3.txt", L, M, span, curfile);	//CHANGEABLE
		FILE* fp = fopen(FileName, "r");
		if(fp != NULL){
			printf("Opened '%s'\n", FileName);
		}
		else{
			printf("Couldn't open '%s', fp is pointing at NULL. QUITTING!\n", FileName);
			exit(1);
		}

		//get rid of UofS line
		fgets(line, 7, fp);

		//Get polygon
		//read in the next polygon
		while( fgets(line, 7, fp)!=NULL ){
			sscanf(line, "%d %d %d", &x, &y, &z);
			if(x==-999){
				break;
			}
			curpoly++;

			//intialize polygon[] to zeros
			for(i=0; i<=(4+(L+1)*(M+1)*(span+1)-1); i++){
				polygon[i]=0;
			}
			polygon[0] = x;
			polygon[1] = y;
			polygon[2] = z;

			direction=0;
			curindex=3;
			while(1){
				fgets(line, 7, fp);
				sscanf(line, "%d", &direction);
				if(direction==-111){
					break;
				}
				else{
					polygon[curindex]=direction;
					curindex++;
				}
			}
			//polygon[] is set.
//			printf("polygon is set, entire length=%d\n", curindex-3);

			//count how many edges are in the pattern (if z=0 or z=span and direc is in x or y direction, then don't count.
			//x,y,z is polygon[0,1,2]
//			printf("z=%d\n", z);
			i=3;
			length=0;
//			printf("\n");
			while(polygon[i]>0){
//				printf("direc=%d, z=%d\n", polygon[i], z);
				if(z>0 && z<span){	//if not in start/finish hinge, then count the edge.
					length++;
//					printf("counted\n");
				}
				else if(polygon[i]==5 || polygon[i]==6){	//else we are in either the start/finish hinge. only count if it's a z-direction.
					length++;
//					printf("counted\n");
				}
				if(polygon[i]==5 ){
					z++;
				}
				else if(polygon[i]==6){
					z--;
				}
				i++;
			}
			
//			printf("pattern length=%d\n", length);
			counts[length]++;

			if((curindex-3)-length<2){	//check
				printf("diff <2, quitting\n");
				exit(1);
			}
		}

		fclose(fp);
	}

	printf("Polygons looked at: %d\n", curpoly);

	for(i=4; i<=(L+1)*(M+1)*(span+1)-1; i++){
		if(counts[i]!=0){
			printf("Patterns with %d edges: %d\n", i, counts[i]);
		}		
	}

	printf("\nmain ending\n\n");

}




