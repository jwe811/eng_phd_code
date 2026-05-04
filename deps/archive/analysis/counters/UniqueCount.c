// Counts number of unique saps in a set of files
// (assumes some rule for startpoint orientation, so just looks at directions)

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>



int main(void) {
	printf("\nmain starting\n\n");

	const int L = 2;	//CHANGEABLE
	const int M = 1;	//CHANGEABLE
	const int span = 4;	//CHANGEABLE

	const int max_unique = 64558;	//CHANGEABLE
//	const int max_unique = 4580;	//CHANGEABLE
//	const int max_unique = 7631;	//CHANGEABLE

	int i, j, k;

	//count how many polygons are in the file
	char filename[100];
	sprintf(filename, "MCpolysHamL%dM%dspan%dnum1.txt", L, M, span);	//CHANGEABLE
	FILE* fpOriginal = fopen(filename, "r");
	if(fpOriginal != NULL){
		printf("Opened '%s'\n", filename);
	}
	else{
		printf("Couldn't open '%s', fpOriginal is pointing at NULL. QUITTING!\n", filename);
		exit(1);
	}

	char line[7];
	int totpolys=0;
	int x, y, z;
	while( fgets(line, 7, fpOriginal)!=NULL ){
		sscanf(line, "%d %d %d", &x, &y, &z);
		if(x==-111){
			totpolys++;
		}
	}
	printf("There are %d polygons in the file %s.\n", totpolys, filename);

	//read in and store unique polygons in array polygon[][]
//	printf("totpolys=%d, size=%d\n", totpolys, 3+(L+1)*(M+1)*(span+1));

	int polygon[max_unique][3+(L+1)*(M+1)*(span+1)];	//first 3 entries (0,1,2) will be the starting xyz (will hold unique polys)
	int curpolygon[3+(L+1)*(M+1)*(span+1)];
	int counts[max_unique];	//counts how many time the unique poly occurs

//	printf("polygon[][] created\n");
	//intialize polygon[][] to all -1
	//initialize counts[] to all zeros
	for(i=0; i<=max_unique-1; i++){
		counts[i]=0;
		for(j=0; j<=3+(L+1)*(M+1)*(span+1)-1; j++){
			polygon[i][j]=-1;
		}
	}

//	printf("polygon[][] initialized to zeros\n");


	int direction=0;
	int curpoly=-1;
	int curindex;

	int num_duplicates=0;
	int num_unique=0;

	//read in polygons
	fseek(fpOriginal, 0, SEEK_SET);
	fgets(line, 7, fpOriginal);	//get rid of UofS line
	
	while( fgets(line, 7, fpOriginal)!=NULL ){
		sscanf(line, "%d %d %d", &x, &y, &z);
		if(x==-999){
			break;
		}
		curpoly++;
		for(j=0; j<=3+(L+1)*(M+1)*(span+1)-1; j++){
			curpolygon[j]=-1;
		}
		curpolygon[0] = x;
		curpolygon[1] = y;
		curpolygon[2] = z;

		direction=0;
		curindex=3;
		while(1){
			fgets(line, 7, fpOriginal);
			sscanf(line, "%d", &direction);
			if(direction==-111){
				break;
			}
			else{
				curpolygon[curindex]=direction;
				curindex++;
			}
		}

//		printf("Checking poly %d: ", i);
//		for(k=3; k<=(3+(L+1)*(M+1)*(span+1)-1); k++){
//			printf("%d ", curpolygon[k]);
//		}
//		printf("\n");

		//check if curpolygon is unique
		int match=0;	//0 no match, 1 match

		//compare curpolygon with all polygons in polygon[][]
		for(j=0; j<=max_unique-1; j++){	//for each polygon in polygon[][]
			if(polygon[j][0]==-1){
				break;
			}

			match=1;	//assume match
			for(k=0; k<=(3+(L+1)*(M+1)*(span+1)-1); k++){	//for each coordinate/direction
				if(curpolygon[k] != polygon[j][k]){
					match=0;
					break;
				}
			}
			if(match==1){
				break;
			}
		}
		if(match==1){
			num_duplicates++;
			counts[j]++;
		}
		else{
			num_unique++;
			j=0;
			while(polygon[j][0]!=-1){
				j++;
			}
			for(k=0; k<=(3+(L+1)*(M+1)*(span+1)-1); k++){	//for each coordinate/direction
				polygon[j][k] = curpolygon[k];
			}
			counts[j]++;
		}
	}

	printf("%d polys in file\n", totpolys);
	printf("%d duplicates\n", num_duplicates);
	printf("%d unique\n", num_unique);

	printf("COUNTS:\n");
	for(i=0; i<=max_unique-1; i++){
		printf("counts[%d]=%d\n", i, counts[i]);
	}



	printf("\nmain finished\n\n");
	return(0);
}



































