//Takes a polygon file (with the patterns of interest) and the corresponding knot ID file (from the arcs from the patterns of interest)
//will split the polygons in the polygon file into two files, one with locally knotted polygons (patterns), and one with globally knotted polygons (patterns)

//changeskies if you used shrinkknotID.c instead of normal knot ID

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
	const int span = 8;	//CHANGEABLE
	const int totfiles = 1; //CHANGEABLE

	int i, j, k;
	int curfile;

	char InterestName[100];	//files with the polygons (patterns) of interest
	char KnotIDName[45];	//files with the knot ID's corresponding to arcs from polygons (patterns) of interest
	char LocalName[45];		//files where local knots will go
	char GlobalName[45];	//files where global knots will go

		char line[7];		//line in Interest
		char IDone[10];	//first of pair in knotID
		char IDtwo[10];	//second of pair in knotID

		int x, y, z;
		int polygon[4+(L+1)*(M+1)*(span+1)];	//first 3 entries (0,1,2) will be the starting xyz
		int curpoly=0;
		int direction=0;
		int curindex;
		int localcount=0;
		int globalcount=0;

	for(curfile=1; curfile<=totfiles; curfile++){
		//file with original polygons
		sprintf(InterestName, "ImpCreatorTwoStringsPrime_PKshrunkL%dM%dspan%dnum%d.txt", L, M, span, curfile);	//CHANGEABLE
		FILE* fpInterest = fopen(InterestName, "r");
		if(fpInterest != NULL){
			printf("Opened '%s'\n", InterestName);
		}
		else{
			printf("Couldn't open '%s', fpInterest is pointing at NULL. QUITTING!\n", InterestName);
			exit(1);
		}

		//CHANGESKIES if used shrunk_ID
		//file with knot IDs of arcs
		sprintf(KnotIDName, "ArcCheckTwoStringsPrime_PKshrunkL%dM%dspan%dnum%d.txt_shrunkID", L, M, span, curfile);	//CHANGEABLE
		FILE* fpKnotID = fopen(KnotIDName, "r");
		if(fpKnotID != NULL){
			printf("Opened '%s'\n", KnotIDName);
		}
		else{
			printf("Couldn't open '%s', fpKnotID is pointing at NULL. QUITTING!\n", KnotIDName);
			exit(1);
		}

		//file where local knots will go
		sprintf(LocalName, "LocalTwoStringsPrime_PKshrunkL%dM%dspan%dnum%d.txt", L, M, span, curfile);	//CHANGEABLE
		FILE* fpLocal = fopen(LocalName, "w");
		if(fpLocal != NULL){
			printf("Created '%s'\n", LocalName);
		}
		else{
			printf("Couldn't open '%s', fpLocal is pointing at NULL. QUITTING!\n", LocalName);
			exit(1);
		}
		fprintf(fpLocal, "UofS\n");



		//file where global knots will go
		sprintf(GlobalName, "GlobalTwoStringsPrime_PKshrunkL%dM%dspan%dnum%d.txt", L, M, span, curfile);	//CHANGEABLE
		FILE* fpGlobal = fopen(GlobalName, "w");
		if(fpGlobal != NULL){
			printf("Created '%s'\n", GlobalName);
		}
		else{
			printf("Couldn't open '%s', fpGlobal is pointing at NULL. QUITTING!\n", GlobalName);
			exit(1);
		}
		fprintf(fpGlobal, "UofS\n");




		//get rid of UofS line
		fgets(line, 7, fpInterest);


		//Get polygon from Interest

		//read in the next polygon
		while( fgets(line, 7, fpInterest)!=NULL ){
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
				fgets(line, 7, fpInterest);
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

			//Determine which file to print polygon[] to
			fgets(IDone, 10, fpKnotID);
			fgets(IDtwo, 10, fpKnotID);

			//if both say unknot, put in global... else put in local
			if(strcmp(IDone, "unknot\n")==0 && strcmp(IDtwo, "unknot\n")==0){					//CHANGESKIES
	//			printf("Both arcs are unknots, this polygon is globally knotted\n");
				globalcount++;
				//print polygon[] into fpGlobal
				fprintf(fpGlobal, "%d %d %d\n", polygon[0], polygon[1], polygon[2]);
				i=3;
				while(polygon[i] != 0){
					fprintf(fpGlobal, "%d\n", polygon[i]);
					i++;
				}
				fprintf(fpGlobal, "-111\n");
			}
			else{
	//			printf("This polygon is locally knotted\n");
				localcount++;
				//print polygon[] into fpLocal
				fprintf(fpLocal, "%d %d %d\n", polygon[0], polygon[1], polygon[2]);
				i=3;
				while(polygon[i] != 0){
					fprintf(fpLocal, "%d\n", polygon[i]);
					i++;
				}
				fprintf(fpLocal, "-111\n");
			}

		}







		fprintf(fpLocal, "-999");
		fprintf(fpGlobal, "-999");


		fclose(fpInterest);
		fclose(fpKnotID);
		fclose(fpLocal);
		fclose(fpGlobal);
	}

	printf("globalcount=%d, localcount=%d\n", globalcount, localcount);
	printf("\nmain ending\n\n");

}




