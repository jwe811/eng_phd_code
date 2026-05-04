/*This function will take a file (or set of files) of two-string pattern polygons, and it will find
two-string patterns that are full BUT cannot occur in a Hamiltonian SAP. The condition for this is
that some first/last sections (which are two-strings) cannot be the first/last sections in HAM polys
(which depends what tube we're looking at). This function will output the full but not ham
two-string patterns/polys to a file named TwoStringsFullButNotHamL%dM%dspan%d.txt
*/

//assumes tube is in the x-direction, y corresponds to L, z corresponds to M.

//as of now, this will only work for the L=3, M=1 tube.

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>

#define	L 2
#define	M 1
#define	totalspan 8	//span of polygons being looked at

#define	vL (L+1)
#define	vM (M+1)
#define	maxsize (vL*vM*(totalspan+1)+3)

//CHANGEABLE THINGS HERE: (ctrl+F changeable to find other changeable things)
//GLOBAL VARIABLES
int totalfilenum = 1; //total number of files that need to be read in.
int ordertemplate[2][M+1][L+1];	/* The first index is for side (left = 0 ; right = 1) */


//function
unsigned long int	num_section(int);




int main(void)
{

	printf("Finding full but not ham 2-string patterns\n");

	//CHANGEABLE
	int num_invalid=14;
	int invalid[num_invalid];
	
	if(L==3 && M==1){
		invalid[0]=4445;
		invalid[1]=38886;
		invalid[2]=6667;
		invalid[3]=54440;
		invalid[4]=1112;
		invalid[5]=15555;
		invalid[6]=11111;
		invalid[7]=32220;
		invalid[8]=12222;
		invalid[9]=39997;
		invalid[10]=13333;
		invalid[11]=47774;
		invalid[12]=10000;
		invalid[13]=24443;
		invalid[14]=22221;
		invalid[15]=56662;
		invalid[16]=21110;
		invalid[17]=48885;
		invalid[18]=19999;
		invalid[19]=41108;
		invalid[20]=28887;
		invalid[21]=49996;
		invalid[22]=26665;
		invalid[23]=34442;
		invalid[24]=36664;
		invalid[25]=51107;
		invalid[26]=45552;
		invalid[27]=59995;
	}
	else if(L==2 && M==1){
		invalid[0]=38;
		invalid[1]=371;
		invalid[2]=112;
		invalid[3]=741;
		invalid[4]=260;
		invalid[5]=593;
		invalid[6]=297;
		invalid[7]=778;
		invalid[8]=334;
		invalid[9]=963;
		invalid[10]=482;
		invalid[11]=815;
		invalid[12]=704;
		invalid[13]=1037;
	}
	else{
		printf("The tube size L=%d M=%d is not supported for this function. Returning\n");
		exit(1);
	}

	int i, j, k;
	int firstsec;
	int lastsec;

	int polycount=0;
	int NotFullCount=0;
	int FullNotHamCount=0;
	int FullHamCount=0;

	int Edges[2];
//	Edges[0] = #edges in first_sec
//	Edges[1] = #edges in last_sec

	int x, y, z;
	int direc;
	char line[7];

	//point at file and read it
	FILE* fp;
	FILE* fp2;
	FILE* fp3;
	FILE* fp4;
	int filenum = 1;

	char filename[50];
	char outputname[50];

	int occupied[M+1][L+1][totalspan+1]; //will indicate whether what vertices are filled. All intermediate vertices (z=1 to totalspan-1) must be occupied
	int filled;
	int foundinvalid;

	int polygon[maxsize];
	int edgecount;

	while(filenum<=totalfilenum){	//while there are files to read in
		//CHANGEABLE
		sprintf(filename, "ImpCreatorTwoStrings_PKshrunkL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);

		fp = fopen(filename, "r");	//open file to read.
		if(fp != NULL){
			printf("Opened '%s'\n", filename);
		}
		else{
			printf("file pointer is pointing to NULL, couldn't open '%s'. Exiting program.\n", filename);
			exit(1);
		}

		//CHANGEABLE
		sprintf(outputname, "TwoStringsFullButNotHamL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);

		fp2 = fopen(outputname, "w");	//open file to write.
		if(fp2 != NULL){
			fprintf(fp2, "UofS\n");	//first line in file is always "UofS"
			printf("Opened '%s' and printed UofS\n", outputname);
		}
		else{
			printf("file pointer is pointing to NULL, couldn't open '%s'. Exiting program.\n", outputname);
			exit(1);
		}

		//CHANGEABLE
		sprintf(outputname, "TwoStringsFullAndHamL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);

		fp3 = fopen(outputname, "w");	//open file to write.
		if(fp3 != NULL){
			fprintf(fp3, "UofS\n");	//first line in file is always "UofS"
			printf("Opened '%s' and printed UofS\n", outputname);
		}
		else{
			printf("file pointer is pointing to NULL, couldn't open '%s'. Exiting program.\n", outputname);
			exit(1);
		}

		//CHANGEABLE
		sprintf(outputname, "TwoStringsNotFullL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);

		fp4 = fopen(outputname, "w");	//open file to write.
		if(fp4 != NULL){
			fprintf(fp4, "UofS\n");	//first line in file is always "UofS"
			printf("Opened '%s' and printed UofS\n", outputname);
		}
		else{
			printf("file pointer is pointing to NULL, couldn't open '%s'. Exiting program.\n", outputname);
			exit(1);
		}


		fgets(line, 7, fp);	//read first line to get it out of the way
//		printf("first line=%s", line);

		/*ALGORITHM:
			for each polygon, walk through polygon and at each step:
				fill ordertemplate[2][vM][vL]
					ordertemplate[0] is first section
					ordertemplate[1] is last section
				fill	occupied[][][]
		*/

		while( fgets(line, 7, fp)!=NULL ){	//while there's a line to read
			sscanf(line, "%d %d %d", &x, &y, &z);	//read in starting point
//			printf("starting point is %d, %d, %d\n", x, y, z);
			if(x==-999){	//end of file
				break;
			}

			polycount++;

			//reset Edges[]
			Edges[0]=0;
			Edges[1]=0;

			//reset occupied[][][]
			for(i=0; i<=M; i++){
				for(j=0; j<=L; j++){
					for(k=0; k<=totalspan; k++){
						occupied[i][j][k]=0;	//initialize to zero.
					}
				}
			}

			//reset ordertemplate[]
			for(i=0; i<=M; i++){
				for(j=0; j<=L; j++){
					ordertemplate[0][i][j]=0;
					ordertemplate[1][i][j]=0;
				}
			}

			//reset polygon and edgecount
			polygon[0]=x;
			polygon[1]=y;
			polygon[2]=z;
			for(i=3; i<=maxsize-1; i++){
				polygon[i]=0;
			}
			edgecount=0;

			direc=0;
			while(1){	//while still in the polygon

				fgets(line, 7, fp);
				sscanf(line, "%d", &direc);	//get direction

				if(direc==-111){
					break;	//break if end of the polygon
				}

				//ordertemplate[0] or ordertemplate[1] changed? (will do Edges[2] later)
				if(x==0 && direc==1){	//edge is in first section
					Edges[0]++;
					ordertemplate[0][z][y] = Edges[0];
				}
				else if(x==1 && direc==2){	//edge is in first section
					Edges[0]++;
					ordertemplate[0][z][y] = Edges[0];
				}
				else if(x==(totalspan-1) && direc==1){	//edge is in last section
					Edges[1]++;
					ordertemplate[1][z][y] = Edges[1];
				}
				else if(x==(totalspan) && direc==2){	//edge is in last section
					Edges[1]++;
					ordertemplate[1][z][y] = Edges[1];
				}

				//Now, update x, y, z, occupied, and polygon/edgecount
				if(direc==1){
					x=x+1;
					occupied[z][y][x]=1;
				}
				else if(direc==2){
					x=x-1;
					occupied[z][y][x]=1;
				}
				else if(direc==3){
					y=y+1;
					occupied[z][y][x]=1;
				}
				else if(direc==4){
					y=y-1;
					occupied[z][y][x]=1;
				}
				else if(direc==5){
					z=z+1;
					occupied[z][y][x]=1;
				}
				else if(direc==6){
					z=z-1;
					occupied[z][y][x]=1;
				}
				polygon[edgecount+3]=direc;
				edgecount++;
			}

			//polygon is done

			filled=1;
			for(i=0; i<=M; i++){
				for(j=0; j<=L; j++){
					for(k=1; k<=totalspan-1; k++){
						if(occupied[i][j][k]==0){	// an empty intermediate vertex
							filled=0;
							break;
						}
					}
					if(filled==0){
						break;
					}
				}
				if(filled==0){
					break;
				}
			}

			if(filled==1){
//				printf("SAP %d is full.\n", polycount);

				//get first and last section numbers
				firstsec = num_section(0);
				lastsec = num_section(1);
//				printf("First sec = %d, last sec = %d\n", firstsec, lastsec);

				//check if either of these sections are invalid.
				foundinvalid=0;
				for(i=0; i<=num_invalid-1; i++){
					if(firstsec==invalid[i] || lastsec==invalid[i]){
						printf("SAP %d is a full pattern that is not HAM. Could print to file here.\n", polycount);
						foundinvalid=1;
						FullNotHamCount++;

						//print to file
						fprintf(fp2, "%d %d %d\n", polygon[0], polygon[1], polygon[2]);
						j=3;
						while(polygon[j]!=0){
							fprintf(fp2, "%d\n", polygon[j]);
							j++;
						}
						fprintf(fp2, "-111\n");				

						break;
					}
				}
				if(foundinvalid==0){
					printf("SAP %d is HAM\n", polycount);
					FullHamCount++;

					//print to file
					fprintf(fp3, "%d %d %d\n", polygon[0], polygon[1], polygon[2]);
					j=3;
					while(polygon[j]!=0){
						fprintf(fp3, "%d\n", polygon[j]);
						j++;
					}
					fprintf(fp3, "-111\n");	
				}
			}
			else{ //don't have to worry about this pattern
				printf("This SAP (%d) is not full. Ignoring\n", polycount);
				NotFullCount++;

				//print to file
				fprintf(fp4, "%d %d %d\n", polygon[0], polygon[1], polygon[2]);
				j=3;
				while(polygon[j]!=0){
					fprintf(fp4, "%d\n", polygon[j]);
					j++;
				}
				fprintf(fp4, "-111\n");	
			}
		}


		filenum++;
	}
	printf("polycount=%d\n\n", polycount);

	printf("NotFullCount=%d\n", NotFullCount);
	printf("FullHamCount=%d\n", FullHamCount);
	printf("FullNotHamCount=%d\n", FullNotHamCount);

	fprintf(fp2, "-999\n");
	fprintf(fp3, "-999\n");
	fprintf(fp4, "-999\n");
	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	printf("first_last_sec_ham_check is FINISHED!\n");


	return 1;
}



#if vM*vL<7
#include "../../sections/Num_section_6V.c"  /*unsigned long int num_section(int side)*/
/* This function takes the arguement (int side) and assigns a unique number */
/* to the section on side = side */
#elif vM*vL<9
#include "../../sections/Num_section_8V.c"
#elif vM*vL<11
#include "../../sections/Num_section_10V.c"
#endif






