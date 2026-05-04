//This program will create a file "Interest_Locally_KnottedL3M1span5.txt" which will contain the polygons which are considered locally knotted from the knotted ID program.



#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>



main(void)
{
	printf("\nmain starting\n\n");

	int i;

	//file with the polygons of interest
	char filenameInterest[45];
	sprintf(filenameInterest, "InterestL3M1span5num1_3_1.txt");	//CHANGEABLE
	FILE* fpInterest = fopen(filenameInterest, "r");
	if(fpInterest != NULL){
		printf("Opened '%s'\n", filenameInterest);
	}
	else{
		printf("Couldn't open '%s', fpInterest is pointing at NULL. QUITTING!\n", filenameInterest);
		exit(1);
	}

	//file with knot ID's of arcs
	char filenameArcCheckID[45];
	sprintf(filenameArcCheckID, "ArcCheckL3M1span5num1_3_1.txt_knot");	//CHANGEABLE
	FILE* fpArcCheckID = fopen(filenameArcCheckID, "r");
	if(fpArcCheckID != NULL){
		printf("Opened '%s'\n", filenameArcCheckID);
	}
	else{
		printf("Couldn't open '%s', fpArcCheckID is pointing at NULL. QUITTING!\n", filenameArcCheckID);
		exit(1);
	}

	//file to write to
	char filename[45];
	sprintf(filename, "Interest_Locally_KnottedL3M1span5.txt");	//CHANGEABLE
	FILE* fp = fopen(filename, "w");
	if(fp != NULL){
		printf("Opened '%s'\n", filename);
	}
	else{
		printf("Couldn't open '%s', fp is pointing at NULL. QUITTING!\n", filename);
		exit(1);
	}
	fprintf(fp, "UofS\n");




	//ALGORITHM
	//read in polygon of interest
	//read in two knot types
	//if one of them is 3_1, print corresponding polygon of interest

	int polygon[4+(3+1)*(1+1)*(5+1)];	//first 3 entries (0,1,2) will be the starting xyz
	int polynum=0;
	int arcnum=0;
	int x, y, z;
	int direction;
	int curindex;
	int KIone;
	int KItwo;

	char lineInterest[7];	//line from fpInterest
	char lineIDone[7];		//1st line from fpArcCheckID
	char lineIDtwo[7];		//2nd line from fpArcCheckID

	fgets(lineInterest, 7, fpInterest);	//get UofS out of the way.

	while( fgets(lineIDone, 7, fpArcCheckID)!=NULL ){	//first arc knot type
		arcnum++;
		fgets(lineIDtwo, 7, fpArcCheckID);		//second arc knot type
		arcnum++;
		
		//read in polygon
		fgets(lineInterest, 7, fpInterest);	//first line (startpoint)
		sscanf(lineInterest, "%d %d %d", &x, &y, &z);
		if(x==-999){
			printf("Problem, reached end of interst file before the end of the knotID file. Quitting\n");
			printf("polynum=%d, arcnum=%d\n", polynum, arcnum);
			exit(1);
		}
		polynum++;

		//intialize polygon[] to zeros
		for(i=0; i<=(4+(3+1)*(1+1)*(5+1)-1); i++){
			polygon[i]=0;
		}
		polygon[0] = x;
		polygon[1] = y;
		polygon[2] = z;

		direction=0;
		curindex=3;
		while(1){
			fgets(lineInterest, 7, fpInterest);
			sscanf(lineInterest, "%d", &direction);
			if(direction==-111){
				break;
			}
			else{
				polygon[curindex]=direction;
				curindex++;
			}
		}
		//polygon[] is set.

		sscanf(lineIDone, "%d_1", &KIone);
		sscanf(lineIDtwo, "%d_1", &KItwo);
		if(KIone==3 || KItwo==3){	//then print polygon[] to file
			fprintf(fp, "%d %d %d\n", polygon[0], polygon[1], polygon[2]);
			i=3;
			while(polygon[i]!=0){
				fprintf(fp, "%d\n", polygon[i]);
				i++;
			}
			fprintf(fp, "-111\n");
		}

	}



	printf("polynum=%d, arcnum=%d\n", polynum, arcnum);
	fprintf(fp, "-999\n");
	fclose(fp);

}



























