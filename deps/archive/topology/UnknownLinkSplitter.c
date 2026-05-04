//Takes a Link file and the corresponding Link ID file. Will create a new file with only the "unknown" polys

//test with "MC2SAPsHamL2M1span90run1num*_LP.txt"
//ID files in folder "knot_lists/MC2SAPsHamL2M1span90run1num*_LP.txt_knot_new"

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
	const int span = 90;	//CHANGEABLE
	const int totfiles = 100; //CHANGEABLE

	unsigned long int LineLength = (L+1)*(M+1)*(span+1)*6+1;

	int curfile;

	char Filename[100];		//files with the links
	char LinkIDname[100];	//files with the link IDs
	char Unknownname[100];	//file where the unknown links will go

	char* line;		//poly in LP format
	line = (char*)malloc(sizeof(char)*(LineLength));
	if(line==NULL){
		fprintf(stderr, "LP line too big. Exitting\n");
		exit(1);
	}
	char ID[20];		//line in linkID file
	char ID2[20];		//next line in linkID file

	int foundUnknown;
	int IDcounter=0;
	int UnknownCounter=0;
	unsigned long int linecounter=0;

	FILE* fpLinkfile;
	FILE* fpLinkIDfile;
	FILE* fpUnknownfile;

	for(curfile=1; curfile<=totfiles; curfile++){
		//file with original links
		sprintf(Filename, "MC2SAPsHamL%dM%dspan%drun1num%d_LP.txt", L, M, span, curfile);	//CHANGEABLE
		fpLinkfile = fopen(Filename, "r");
		if(fpLinkfile != NULL){
			printf("Opened '%s'\n", Filename);
		}
		else{
			printf("Couldn't open '%s', fpLinkfile is pointing at NULL. QUITTING!\n", Filename);
			exit(1);
		}

		//file with knot IDs of polys
		sprintf(LinkIDname, "knot_lists/MC2SAPsHamL%dM%dspan%drun1num%d_LP.txt_knot_new", L, M, span, curfile);	//CHANGEABLE
		fpLinkIDfile = fopen(LinkIDname, "r");
		if(fpLinkIDfile != NULL){
			printf("Opened '%s'\n", LinkIDname);
		}
		else{
			printf("Couldn't open '%s', fpLinkIDfile is pointing at NULL. QUITTING!\n", LinkIDname);
			exit(1);
		}

		//file where knotted polys will go
		sprintf(Unknownname, "MC2SAPsHamL%dM%dspan%drun1num%d_LP_unknown.txt", L, M, span, curfile);	//CHANGEABLE
		fpUnknownfile = fopen(Unknownname, "w");
		if(fpUnknownfile != NULL){
			printf("Created '%s'\n", Unknownname);
		}
		else{
			printf("Couldn't open '%s', fpUnknownfile is pointing at NULL. QUITTING!\n", Unknownname);
			exit(1);
		}

		//Read in ID and ID2
		fgets(ID, 20, fpLinkIDfile);
		while( fgets(ID2, 20, fpLinkIDfile)!=NULL ){
			IDcounter++;
			foundUnknown=0;
//			printf("\nID=%s", ID);
			if(strcmp(ID, "unknown\n")==0){
				foundUnknown=1;
				UnknownCounter++;
//				printf("unknown %d found in 2SAP %d\n", UnknownCounter, IDcounter);
			}

			//read in 2 polys. if unknown, then print them.
			fgets(line, LineLength, fpLinkfile);
//			printf("Read in a poly: %s", line);
			linecounter++;
			if(foundUnknown==1){
//				printf("printing line %lu\n", linecounter);
				fprintf(fpUnknownfile, "%s", line);
			}
			fgets(line, LineLength, fpLinkfile);
//			printf("Read in a poly: %s", line);
			linecounter++;
			if(foundUnknown==1){
//				printf("printing line %lu\n", linecounter);
				fprintf(fpUnknownfile, "%s", line);
			}

			//copy ID2 into ID
			strcpy(ID, ID2);
		}
				



		fclose(fpLinkfile);
		fclose(fpLinkIDfile);
		fclose(fpUnknownfile);
	}

	printf("IDcounter=%d\n", IDcounter);
	printf("UnknownCounter=%d\n", UnknownCounter);
	printf("\nmain ending\n\n");

	return 0;

}




