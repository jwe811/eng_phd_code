//Takes a polygon file and the corresponding knot ID file. Will create a new file with only the knotted polys

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

	const int L = 5;	//CHANGEABLE
	const int M = 1;	//CHANGEABLE
	const int span = 50;	//CHANGEABLE
	const int totfiles = 1; //CHANGEABLE

	int i, j, k;
	int curfile;

	char Filename[100];		//files with the polygons
	char KnotIDName[100];	//files with the knot IDs
	char Knottedname[100];	//file where the knotted polygons will go

	char line[8];		//line in poly file
	char ID[8];		//line in knotID file

	int x, y, z;
	int polygon[4+(L+1)*(M+1)*(span+1)];	//first 3 entries (0,1,2) will be the starting xyz
	int curpoly=0;
	int direction=0;
	int curindex;
	int knottedcount=0;
	int unknottedcount=0;

	FILE* fpPolyfile;
	FILE* fpKnotIDfile;
	FILE* fpKnottedfile;

	for(curfile=1; curfile<=totfiles; curfile++){
		//file with original polygons
		sprintf(Filename, "MCpolysHamL%dM%dspan%drun1num%d.txt", L, M, span, curfile);	//CHANGEABLE
		fpPolyfile = fopen(Filename, "r");
		if(fpPolyfile != NULL){
			printf("Opened '%s'\n", Filename);
		}
		else{
			printf("Couldn't open '%s', fpPolyfile is pointing at NULL. QUITTING!\n", Filename);
			exit(1);
		}

		//CHANGESKIES if used shrunk_ID
		//file with knot IDs of polys
		sprintf(KnotIDName, "MCpolysHamL%dM%dspan%drun1num%d.txt_knot", L, M, span, curfile);	//CHANGEABLE
		FILE* fpKnotIDfile = fopen(KnotIDName, "r");
		if(fpKnotIDfile != NULL){
			printf("Opened '%s'\n", KnotIDName);
		}
		else{
			printf("Couldn't open '%s', fpKnotIDfile is pointing at NULL. QUITTING!\n", KnotIDName);
			exit(1);
		}

		//file where knotted polys will go
		sprintf(Knottedname, "MCpolysHam_Knotted_L%dM%dspan%dnum%d.txt", L, M, span, curfile);	//CHANGEABLE
		FILE* fpKnottedfile = fopen(Knottedname, "w");
		if(fpKnottedfile != NULL){
			printf("Created '%s'\n", Knottedname);
		}
		else{
			printf("Couldn't open '%s', fpKnottedfile is pointing at NULL. QUITTING!\n", Knottedname);
			exit(1);
		}
		fprintf(fpKnottedfile, "UofS\n");




		//get rid of UofS line
		fgets(line, 7, fpPolyfile);

		//Get polygon

		//read in the next polygon
		while( fgets(line, 7, fpPolyfile)!=NULL ){
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
				fgets(line, 7, fpPolyfile);
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
			fgets(ID, 8, fpKnotIDfile);

			//if both say unknot, put in global... else put in local
			if(strcmp(ID, "0_1\n")!=0){					//CHANGESKIES
	//			printf("Polygon is knotted\n");
				knottedcount++;
				//print polygon[] into fpKnottedfile
				fprintf(fpKnottedfile, "%d %d %d\n", polygon[0], polygon[1], polygon[2]);
				i=3;
				while(polygon[i] != 0){
					fprintf(fpKnottedfile, "%d\n", polygon[i]);
					i++;
				}
				fprintf(fpKnottedfile, "-111\n");
			}
			else{
				unknottedcount++;
			}
		}







		fprintf(fpKnottedfile, "-999");


		fclose(fpPolyfile);
		fclose(fpKnotIDfile);
		fclose(fpKnottedfile);
	}

	printf("knottedcount=%d\n", knottedcount);
	printf("unknottedcount=%d\n", unknottedcount);
	printf("\nmain ending\n\n");

}




