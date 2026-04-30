//This program takes a set of polygon files (i.e. "foo_num1.txt" to "foo_num23.txt") and will knot ID using bfacf_shrink_v2.c
//It will output a text file "foo_num1.txt_shrunkID" with either "unknot" or "knot"


#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>


int bfacf_shrink(int x, int y, int z, int polygon[], int length);

main(void)
{
	printf("\nmain starting\n\n");

	const int L = 2;	//CHANGEABLE
	const int M = 1;	//CHANGEABLE
	const int span = 8;	//CHANGEABLE
	const int num_files = 1;	//CHANGEABLE (number of files)

	int i, j, k;
	int curfile;
	int curpoly=0;	//counts how many polys we shrunkID'd

	for(curfile=1; curfile<=num_files; curfile++){
		//file with the polygons
		char filename[100];
		sprintf(filename, "ArcCheckTwoStringsPrime_PKshrunkL%dM%dspan%dnum%d.txt", L, M, span, curfile);	//CHANGEABLE
		FILE* fp = fopen(filename, "r");
		if(fp != NULL){
			printf("Opened '%s'\n", filename);
		}
		else{
			printf("Couldn't open '%s', fp is pointing at NULL. QUITTING!\n", filename);
			exit(1);
		}

		//file that will contain with knot ID's (knot or unknot)
		char filenametwo[100];
		sprintf(filenametwo, "ArcCheckTwoStringsPrime_PKshrunkL%dM%dspan%dnum%d.txt_shrunkID", L, M, span, curfile);	//CHANGEABLE
		FILE* fptwo = fopen(filenametwo, "w");
		if(fp != NULL){
			printf("Opened '%s'\n", filenametwo);
		}
		else{
			printf("Couldn't open '%s', fp is pointing at NULL. QUITTING!\n", filenametwo);
			exit(1);
		}


		char line[7];
		int x, y, z;
		int polygon[(L+1)*(M+1)*(span+1)+10];	//make a bit bigger to allow for +2 moves
		int direction, curindex;

		//get rid of UofS line
		fgets(line, 7, fp);

		//read in the next polygon
		while( fgets(line, 7, fp)!=NULL ){
			sscanf(line, "%d %d %d", &x, &y, &z);
			if(x==-999){
				break;
			}
			curpoly++;

			//intialize polygon[] to zeros
			for(i=0; i<=((L+1)*(M+1)*(span+1)+10-1); i++){
				polygon[i]=0;
			}

			curindex=0;
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
			//length is curindex
//			printf("length is %d\n", curindex);

			if( bfacf_shrink(x,y,z,polygon,curindex) < 24 ){ //unknotted
				fprintf(fptwo, "unknot\n");
			}
			else{ //knotted
				fprintf(fptwo, "knot\n");
			}
			



		}

		printf("KnotID'd %d polys\n", curpoly);
		fclose(fp);
		printf("file %s closed\n", filename);
		fclose(fptwo);
		printf("file %s closed\n", filenametwo);
		



	}




}

#include "../generators/bfacf_shrink_v2.c"
