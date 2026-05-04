//Reads in knot ID file of Unknown links, and linking number file of unknown links. Prints Summary?

//test with "MC2SAPsHamL2M1span90run1num1_LP_unknown.txt_knot"
//test with "MC2SAPsHamL2M1span90run1num1_LP_unknown_LN.txt"

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>


#define L 2
#define M 1
#define span 90
#define totfiles 100

int main(void) {
	printf("\nmain starting\n\n");

	int curfile;

	char IDname[100];		//files with the Knot IDs
	char LinkNumname[100];	//files with Linking numbers

	FILE* fpID;
	FILE* fpLinkNum;

	char ID1[100];		//line with ID1
	char ID2[100];		//line with ID2
	char LinkNum[100];	//line with LinkingNumber
	int LinkNumInt;

	int LinkedUnknots=0;
	int MaybeUnlinkedUnknots=0;
	int LinkedOneKnot=0;
	int LinkedTwoKnots=0;
	int MaybeUnlinkedOneKnot=0;
	int MaybeUnlinkedTwoKnots=0;
	int UnknownCounter=0;


	for(curfile=1; curfile<=totfiles; curfile++){
		//file with Knot IDs
		sprintf(IDname, "knot_lists/MC2SAPsHamL%dM%dspan%drun1num%d_LP_unknown.txt_knot", L, M, span, curfile);	//CHANGEABLE
		fpID = fopen(IDname, "r");
		if(fpID != NULL){
			printf("Opened '%s'\n", IDname);
		}
		else{
			printf("Couldn't open '%s', fpID is pointing at NULL. QUITTING!\n", IDname);
			exit(1);
		}

		//file with LinkingNums
		sprintf(LinkNumname, "MC2SAPsHamL%dM%dspan%drun1num%d_LP_unknown_LN.txt", L, M, span, curfile);	//CHANGEABLE
		fpLinkNum = fopen(LinkNumname, "r");
		if(fpLinkNum != NULL){
			printf("Opened '%s'\n", LinkNumname);
		}
		else{
			printf("Couldn't open '%s', fpLinkNum is pointing at NULL. QUITTING!\n", LinkNumname);
			exit(1);
		}

		while( fgets(ID1, 100, fpID)!=NULL ){
			fgets(ID2, 100, fpID);
			fgets(LinkNum, 100, fpLinkNum);
			UnknownCounter++;

//			printf("ID1=%s", ID1);
//			printf("ID2=%s", ID2);
//			printf("LinkNum=%s", LinkNum);
			LinkNumInt = atoi(LinkNum);
//			printf("As int=%d\n", LinkNumInt);

			if( strcmp(ID1, "0_1\n")==0 && strcmp(ID2, "0_1\n")==0){	//we have two unknots
				if(LinkNumInt==0){
					MaybeUnlinkedUnknots++;
				}
				else{
					LinkedUnknots++;
				}
			}
			else if((strcmp(ID1, "0_1\n")==0 && strcmp(ID2, "0_1\n")!=0) || (strcmp(ID1, "0_1\n")!=0 && strcmp(ID2, "0_1\n")==0)){	//we have one unknot and one knot
				if(LinkNumInt==0){
					MaybeUnlinkedOneKnot++;
				}
				else{
					LinkedOneKnot++;
				}
			}
			else if( strcmp(ID1, "0_1\n")!=0 && strcmp(ID2, "0_1\n")!=0){	//we have two knots
				if(LinkNumInt==0){
					MaybeUnlinkedTwoKnots++;
				}
				else{
					LinkedTwoKnots++;
				}
			}
			


		}



		fclose(fpID);
		fclose(fpLinkNum);
	}
	
	printf("\nUnknownCounter=%d\n\n", UnknownCounter);

	printf("LinkedUnknots=%d\n", LinkedUnknots);
	printf("LinkedOneKnot=%d\n", LinkedOneKnot);
	printf("LinkedTwoKnots=%d\n", LinkedTwoKnots);
	printf("MaybeUnlinkedUnknots=%d\n", MaybeUnlinkedUnknots);
	printf("MaybeUnlinkedOneKnot=%d\n", MaybeUnlinkedOneKnot);
	printf("MaybeUnlinkedTwoKnots=%d\n", MaybeUnlinkedTwoKnots);

	printf("\nmain ending\n\n");

	return 0;

}






