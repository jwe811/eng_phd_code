//Takes a polygon file and the corresponding knot ID file. Will create 2 new files, one with 3_1 and the other with 3_1s

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

	const int L = 3;	//CHANGEABLE
	const int M = 1;	//CHANGEABLE
	const int span = 3;	//CHANGEABLE
	const int totfiles = 1; //CHANGEABLE

	int i, j, k;
	int curfile;

	char Filename[100];		//files with the polygons
	char KnotIDName[100];	//files with the knot IDs
	char PosTrefname[100];	//file where the knotted polygons will go
	char NegTrefname[100];	//file where the knotted polygons will go
	//used to be Knottedname

	char line[8];		//line in poly file
	char ID[8];		//line in knotID file

	int x, y, z;
	int polygon[4+(L+1)*(M+1)*(span+1)];	//first 3 entries (0,1,2) will be the starting xyz
	int curpoly=0;
	int direction=0;
	int curindex;
	int unknottedcount=0;
	int PosTrefcount=0;
	int NegTrefcount=0;

	FILE* fpPolyfile;
	FILE* fpKnotIDfile;
	FILE* fp3_1file;
	FILE* fp3_1sfile;
	//used to be fpKnottedfile

	for(curfile=1; curfile<=totfiles; curfile++){
		//file with original polygons
		sprintf(Filename, "ImpCreatorPrimePolysL%dM%dspan%dnum%d.txt", L, M, span, curfile);	//CHANGEABLE
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
		sprintf(KnotIDName, "ImpCreatorPrimePolysL%dM%dspan%dnum%d.txt_knot", L, M, span, curfile);	//CHANGEABLE
		fpKnotIDfile = fopen(KnotIDName, "r");
		if(fpKnotIDfile != NULL){
			printf("Opened '%s'\n", KnotIDName);
		}
		else{
			printf("Couldn't open '%s', fpKnotIDfile is pointing at NULL. QUITTING!\n", KnotIDName);
			exit(1);
		}

		//file where 3_1 polys will go
		sprintf(PosTrefname, "ImpCreatorPrimePolysL%dM%dspan%dnum%d_3_1.txt", L, M, span, curfile);	//CHANGEABLE
		fp3_1file = fopen(PosTrefname, "w");
		if(fp3_1file != NULL){
			printf("Created '%s'\n", PosTrefname);
		}
		else{
			printf("Couldn't open '%s', fp3_1file is pointing at NULL. QUITTING!\n", PosTrefname);
			exit(1);
		}
		fprintf(fp3_1file, "UofS\n");

		//file where 3_1s polys will go
		sprintf(NegTrefname, "ImpCreatorPrimePolysL%dM%dspan%dnum%d_3_1s.txt", L, M, span, curfile);	//CHANGEABLE
		fp3_1sfile = fopen(NegTrefname, "w");
		if(fp3_1sfile != NULL){
			printf("Created '%s'\n", NegTrefname);
		}
		else{
			printf("Couldn't open '%s', fp3_1sfile is pointing at NULL. QUITTING!\n", NegTrefname);
			exit(1);
		}
		fprintf(fp3_1sfile, "UofS\n");



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

			//if it says 3_1, put in fp3_1file
			if(strcmp(ID, "3_1\n")==0){					//CHANGESKIES
	//			printf("This is a 3_1\n");
				PosTrefcount++;
				//print polygon[] into fp3_1file
				fprintf(fp3_1file, "%d %d %d\n", polygon[0], polygon[1], polygon[2]);
				i=3;
				while(polygon[i] != 0){
					fprintf(fp3_1file, "%d\n", polygon[i]);
					i++;
				}
				fprintf(fp3_1file, "-111\n");
			}
			else if(strcmp(ID, "3_1s\n")==0){					//CHANGESKIES
	//			printf("This is a 3_1s\n");
				NegTrefcount++;
				//print polygon[] into fp3_1sfile
				fprintf(fp3_1sfile, "%d %d %d\n", polygon[0], polygon[1], polygon[2]);
				i=3;
				while(polygon[i] != 0){
					fprintf(fp3_1sfile, "%d\n", polygon[i]);
					i++;
				}
				fprintf(fp3_1sfile, "-111\n");
			}
			else{
				unknottedcount++;
			}
		}







		fprintf(fp3_1file, "-999");
		fprintf(fp3_1sfile, "-999");


		fclose(fpPolyfile);
		fclose(fpKnotIDfile);
		fclose(fp3_1file);
		fclose(fp3_1sfile);
	}

	printf("PosTrefcount=%d\n", PosTrefcount);
	printf("NegTrefcount=%d\n", NegTrefcount);
	printf("unknottedcount=%d\n", unknottedcount);
	printf("\nmain ending\n\n");

}




