//Takes a Link file (LP format), calculates linking number for each link, and tally's up totals.

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>



/*
Uses: linkingnumber2x1.c
Uses: linkingnumber3x1.c
Uses: linkingnumber4x1.c
Uses: linkingnumber2x2.c

These functions use the global variables:
ordertemplate[side][vM][vL]	(odd numbers mean the go to the right, even to the left)
ordertemplate2[side][vM][vL]
Vedges[M][vL]	(+1 means up, -1 means down)
Vedges2[M][vL]

*/


#define L 2
#define M 2
#define span 31
#define totfiles 100

#define MaxLinkNum 101	//MaxLinkNum and bigger will all be grouped together


///////////STRUCTURES//////////////

struct Point {
	int x;
	int y;
	int z;
};



/////////////END OF STRUCTURES/////

/////////////////GLOBAL VARIABLES///////////////////

int ordertemplate[2][M+1][L+1];
int ordertemplate2[2][M+1][L+1];
short int Vedges[M][L+1];
short int Vedges2[M][L+1];

////////////////END OF GLOBAL VARIABLES//////////////

///////////////////FUNCTIONS//////////////////////

short int linkingnumber(void);

//////////////////END OF FUNCTIONS/////////////////

int main(void) {
	printf("\nmain starting\n\n");

	unsigned long int LineLength = (L+1)*(M+1)*(span+1)*6+1;
	unsigned long int MaxPolyLength = (L+1)*(M+1)*(span+1);

	int curfile;

	char Filename[100];		//files with the links
//	char Outputname[100];	//file where the unknown links will go

	char* line;		//poly in LP format
	char* line2;		//poly2 in LP format

	line = (char*)malloc(sizeof(char)*(LineLength));
	if(line==NULL){
		fprintf(stderr, "LP line too big. Exitting\n");
		exit(1);
	}
	line2 = (char*)malloc(sizeof(char)*(LineLength));
	if(line2==NULL){
		fprintf(stderr, "LP line2 too big. Exitting\n");
		exit(1);
	}

	char* temp;		//used for reading in the polygons
	temp = (char*)malloc(sizeof(char)*(LineLength));
	if(temp==NULL){
		fprintf(stderr, "temp too big. Exitting\n");
		exit(1);
	}

	FILE* fpLinkfile;
//	FILE* fpOutputfile;





	int i, j, k;

	int UnknownCounter=0;
	int length;
	int length2;
	struct Point poly1[MaxPolyLength];
	struct Point poly2[MaxPolyLength];

	int bigtemplate[span+2][M+1][L+1];	//holds all horizontal (x-direction) edges and their directions
								//filled  with 1 means poly1 right, 2 means poly1 left, 3 means poly2 right, 4 means poly2 left.
								//first index=0 will be empty (x=-1 -> x=0)
								//first index=1 is for x=0 -> x=1
								//first index=span is for x=span-1 -> x=span
								//first index=span+1 will be empty (x=span -> x=span+1)

	short int bigVedges[span+1][M][L+1];	//holds vertical (z-direction) edges and their directions
									//filled with +1 means poly1 up, -1 means poly2 down, +2 means poly2 up, -2 means poly2 down.

	short int LinkingNum=0;
	int curspan;

	long int Tally[MaxLinkNum+1];
	for(i=0; i<=MaxLinkNum; i++){
		Tally[i]=0;
	}



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

		//Output file
/*		sprintf(Outputname, "MC2SAPsHamL%dM%dspan%drun1num%d_LP_unknown_LN.txt", L, M, span, curfile);	//CHANGEABLE
		fpOutputfile = fopen(Outputname, "w");
		if(fpOutputfile != NULL){
			printf("Created '%s'\n", Outputname);
		}
		else{
			printf("Couldn't open '%s', fpOutputfile is pointing at NULL. QUITTING!\n", Outputname);
			exit(1);
		}
*/
		//Read in poly1 and poly2
		while( fgets(line, LineLength, fpLinkfile)!=NULL ){
			fgets(line2, LineLength, fpLinkfile);
			UnknownCounter++;

			//Convert line to poly1
			sscanf(line, "%d %[^\n]\n", &length, temp);
//			printf("length=%d\n", length);
//			printf("temp=%s\n", temp);

			for(i=0; i<length; i++){
				sscanf(temp, "%d %d %d %[^\n]\n", &poly1[i].x, &poly1[i].y, &poly1[i].z, temp);
			}

//			for(i=0; i<length; i++){
//				printf("%d %d %d\n", poly1[i].x, poly1[i].y, poly1[i].z);
//			}


			//Convert line2 to poly2
			sscanf(line2, "%d %[^\n]\n", &length2, temp);
//			printf("length2=%d\n", length2);
//			printf("temp=%s\n", temp);

			for(i=0; i<length2; i++){
				sscanf(temp, "%d %d %d %[^\n]\n", &poly2[i].x, &poly2[i].y, &poly2[i].z, temp);
			}

//			for(i=0; i<length2; i++){
//				printf("%d %d %d\n", poly2[i].x, poly2[i].y, poly2[i].z);
//			}


			//Fill giant arrays which represents edges traversed (including the direction)
			//initialize to zeros
			for(i=0; i<span+2; i++){
				for(j=0; j<M+1; j++){
					for(k=0; k<L+1; k++){
						bigtemplate[i][j][k]=0;
					}
				}
			}
			for(i=0; i<span+1; i++){
				for(j=0; j<M; j++){
					for(k=0; k<L+1; k++){
						bigVedges[i][j][k]=0;
					}
				}
			}
			


			//fill for poly1
			for(i=0; i<length-1; i++){
				if(poly1[i].x < poly1[i+1].x){
					bigtemplate[poly1[i+1].x][poly1[i].z][poly1[i].y] = 1;
				}
				else if(poly1[i].x > poly1[i+1].x){
					bigtemplate[poly1[i].x][poly1[i].z][poly1[i].y] = 2;
				}
				else if(poly1[i].z < poly1[i+1].z){
					bigVedges[poly1[i].x][poly1[i].z][poly1[i].y] = 1;
				}
				else if(poly1[i].z > poly1[i+1].z){
					bigVedges[poly1[i].x][poly1[i+1].z][poly1[i].y] = -1;
				}
			}
			//do last point to first point of poly1 also
			if(poly1[length-1].x < poly1[0].x){
				bigtemplate[poly1[0].x][poly1[0].z][poly1[0].y] = 1;
			}
			else if(poly1[length-1].x > poly1[0].x){
				bigtemplate[poly1[length-1].x][poly1[0].z][poly1[0].y] = 2;
			}
			else if(poly1[length-1].z < poly1[0].z){
				bigVedges[poly1[0].x][poly1[length-1].z][poly1[0].y] = 1;
			}
			else if(poly1[length-1].z > poly1[0].z){
				bigVedges[poly1[0].x][poly1[0].z][poly1[0].y] = -1;
			}



			//fill for poly2
			for(i=0; i<length2-1; i++){
				if(poly2[i].x < poly2[i+1].x){
					bigtemplate[poly2[i+1].x][poly2[i].z][poly2[i].y] = 3;
				}
				else if(poly2[i].x > poly2[i+1].x){
					bigtemplate[poly2[i].x][poly2[i].z][poly2[i].y] = 4;
				}
				else if(poly2[i].z < poly2[i+1].z){
					bigVedges[poly2[i].x][poly2[i].z][poly2[i].y] = 2;
				}
				else if(poly2[i].z > poly2[i+1].z){
					bigVedges[poly2[i].x][poly2[i+1].z][poly2[i].y] = -2;
				}
			}
			//do last point to first point of poly2 also
			if(poly2[length2-1].x < poly2[0].x){
				bigtemplate[poly2[0].x][poly2[0].z][poly2[0].y] = 3;
			}
			else if(poly2[length2-1].x > poly2[0].x){
				bigtemplate[poly2[length2-1].x][poly2[0].z][poly2[0].y] = 4;
			}
			else if(poly2[length2-1].z < poly2[0].z){
				bigVedges[poly2[0].x][poly2[length2-1].z][poly2[0].y] = 2;
			}
			else if(poly2[length2-1].z > poly2[0].z){
				bigVedges[poly2[0].x][poly2[0].z][poly2[0].y] = -2;
			}
		
/*
			//print check
			for(i=0; i<span+2; i++){
				printf("i=%d\n", i);
				for(j=0; j<M+1; j++){
					for(k=0; k<L+1; k++){
						printf("%d ", bigtemplate[i][j][k]);
					}
					printf("\n");
				}
				printf("\n\n");
			}
			for(i=0; i<span+1; i++){
				for(j=0; j<M; j++){
					for(k=0; k<L+1; k++){
						bigVedges[i][j][k]=0;
					}
				}
			}

*/




			LinkingNum=0;
			//Fill ordertemplates and Vedges and calculate linking number, for each span
			for(curspan=0; curspan<=span; curspan++){	//for each span

				//initialize to zeros
				for(i=0; i<=1; i++){
					for(j=0; j<=M; j++){
						for(k=0; k<=L; k++){
							ordertemplate[i][j][k]=0;
							ordertemplate2[i][j][k]=0;
						}
					}
				}
				for(j=0; j<M; j++){
					for(k=0; k<=L; k++){
						Vedges[j][k]=0;
						Vedges2[j][k]=0;
					}
				}


				//fill ordertemplates
				for(i=0; i<=M; i++){
					for(j=0; j<=L; j++){
						//fill left side of ordertemplate/ordertemplate2
						if(bigtemplate[curspan][i][j]<3){
							ordertemplate[0][i][j] = bigtemplate[curspan][i][j];
						}
						else{
							ordertemplate2[0][i][j] = bigtemplate[curspan][i][j];
						}

						//fill right side of ordertemplate/ordertemplate2
						if(bigtemplate[curspan+1][i][j]<3){
							ordertemplate[1][i][j] = bigtemplate[curspan+1][i][j];
						}
						else{
							ordertemplate2[1][i][j] = bigtemplate[curspan+1][i][j];
						}
					}
				}

				//fill Vedges
				for(i=0; i<M; i++){
					for(j=0; j<=L; j++){
						//fill Vedges/Vedges2
						if(abs(bigVedges[curspan][i][j])<2){
							Vedges[i][j] = bigVedges[curspan][i][j];
						}
						else{
							Vedges2[i][j] = bigVedges[curspan][i][j];
						}
					}
				}

/*				//print test
				printf("ordtemp[0]:\n");
				for(j=0; j<=M; j++){
					for(k=0; k<=L; k++){
						printf("%d ", ordertemplate[0][j][k]);
					}
					printf("\n");
				}
				printf("\n");
				printf("ordtemp2[0]:\n");
				for(j=0; j<=M; j++){
					for(k=0; k<=L; k++){
						printf("%d ", ordertemplate2[0][j][k]);
					}
					printf("\n");
				}
				printf("\n");
				printf("ordtemp[1]:\n");
				for(j=0; j<=M; j++){
					for(k=0; k<=L; k++){
						printf("%d ", ordertemplate[1][j][k]);
					}
					printf("\n");
				}
				printf("\n");
				printf("ordtemp2[1]:\n");
				for(j=0; j<=M; j++){
					for(k=0; k<=L; k++){
						printf("%d ", ordertemplate2[1][j][k]);
					}
					printf("\n");
				}
				printf("\n");

				printf("Vedges:\n");
				for(j=0; j<M; j++){
					for(k=0; k<=L; k++){
						printf("%d ", Vedges[j][k]);
					}
					printf("\n");
				}
				printf("Vedges2:\n");
				for(j=0; j<M; j++){
					for(k=0; k<=L; k++){
						printf("%d ", Vedges2[j][k]);
					}
					printf("\n");
				}
*/


				//calculate linking number
				LinkingNum+=linkingnumber();
//				printf("LinkingNum is now %d\n", LinkingNum);

			}
			

			if(abs(LinkingNum)>=MaxLinkNum){
				Tally[MaxLinkNum]++;
			}
			else{
				Tally[abs(LinkingNum)]++;
			}
			//print linking number to Outputfile
//			printf("Linking number for 2SAP %d is %d\n", UnknownCounter, LinkingNum);


//			fprintf(fpOutputfile, "%d\n", LinkingNum);			
			






		}
				
/*
		for(i=0; i<MaxLinkNum; i++){
			if(Tally[i]>0){
				printf("Tally[%d] = %li\n", i, Tally[i]);
			}
		}
*/
		fclose(fpLinkfile);
//		fclose(fpOutputfile);
	}

	printf("Links Read = %d\n", UnknownCounter);

	printf("Final Tally:\n");
	for(i=0; i<MaxLinkNum; i++){
		if(Tally[i]>0){
			printf("Tally[%d] = %li\n", i, Tally[i]);
		}
	}
	printf("\nmain ending\n\n");

	return 0;

}


#if M==1 && L==2
	#include "linkingnumber2x1.c"
#elif M==1 && L==3
	#include "linkingnumber3x1.c"
#elif M==1 && L==4
	#include "linkingnumber4x1.c"
#elif M==2 && L==2
	#include "linkingnumber2x2.c"
#endif



