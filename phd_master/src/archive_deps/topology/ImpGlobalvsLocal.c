//This program takes knotted polygons (of interest) and creates "ArcCheck" polygons. For every knotted polygon,
//2 archeck polygons are created from the 2 strands in the pattern of the knotted poly.

//ImpGlobalvsLocal --> tube is in x-direction

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
	const int num_files = 1;	//CHANGEABLE (number of files)

	int i, j, k;
	int curfile;

	for(curfile=1; curfile<=num_files; curfile++){
		//file with the knotted polygons
		char filename[100];
		sprintf(filename, "ImpCreatorTwoStringsPrime_PKshrunkL%dM%dspan%dnum%d.txt", L, M, span, curfile);	//CHANGEABLE
		FILE* fp = fopen(filename, "r");
		if(fp != NULL){
			printf("Opened '%s'\n", filename);
		}
		else{
			printf("Couldn't open '%s', fp is pointing at NULL. QUITTING!\n", filename);
			exit(1);
		}

		//file where the two polygons from the two arcs will be written to.
		char filenametwo[100];
		sprintf(filenametwo, "ArcCheckTwoStringsPrime_PKshrunkL%dM%dspan%dnum%d.txt", L, M, span, curfile);	//CHANGEABLE
		FILE* fptwo = fopen(filenametwo, "w");
		if(fptwo != NULL){
			printf("Created '%s'\n", filenametwo);
		}
		else{
			printf("Couldn't open '%s', fptwo is pointing at NULL. QUITTING!\n", filenametwo);
			exit(1);
		}
		fprintf(fptwo, "UofS\n");

		//file where polygons of interest will be printed to.
		char filenameInterest[100];
		sprintf(filenameInterest, "InterestTwoStringsPrime_PKshrunkL%dM%dspan%dnum%d.txt", L, M, span, curfile);	//CHANGEABLE
		FILE* fpInterest = fopen(filenameInterest, "w");
		if(fpInterest != NULL){
			printf("Created '%s'\n", filenameInterest);
		}
		else{
			printf("Couldn't open '%s', fpInterest is pointing at NULL. QUITTING!\n", filenameInterest);
			exit(1);
		}
		fprintf(fpInterest, "UofS\n");

		char line[7];
		int x, y, z;

		//get rid of UofS line
		fgets(line, 7, fp);

		int polygon[4+(L+1)*(M+1)*(span+1)];	//first 3 entries (0,1,2) will be the starting xyz
		int arcone[(L+1)*(M+1)*(span+1)];		//will store arcs (made the size of the largest possible arc)
		int arctwo[(L+1)*(M+1)*(span+1)];		//first arc[0], arc[1], and arc[2] will be the starting x,y,z of the arc respectively.
		int arconelength, arctwolength;
		int zedgesintop, zedgesinbottom;
		int otherx, othery, otherz;
	

		int direction=0;
		int curpoly=0;
		int curindex;
		int edgesinfirst, edgesinlast;
		int curx, cury, curz;
		int patterncount=0;
		int localcount=0;
		int globalcount=0;

		//read in the next polygon
		while( fgets(line, 7, fp)!=NULL ){
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

			//Check if there are just two edges in the first and last section of the polygon
			//REMEMBER X-DIRECTION IS THE DIRECTION OF THE TUBE
			edgesinfirst=0;
			edgesinlast=0;
			curx = polygon[0];
			cury = polygon[1];
			curz = polygon[2];
			for(i=3; i<=(4+(L+1)*(M+1)*(span+1)-1); i++){
				if(polygon[i]==1){
					if(curx==0){
						edgesinfirst++;
					}
					else if(curx==span-1){
						edgesinlast++;
					}
					curx++;
				}
				else if(polygon[i]==2){
					if(curx==1){
						edgesinfirst++;
					}
					else if(curx==span){
						edgesinlast++;
					}
					curx--;
				}
				else if(polygon[i]==3){
					cury++;
				}
				else if(polygon[i]==4){
					cury--;
				}
				else if(polygon[i]==5){
					curz++;
				}
				else if(polygon[i]==6){
					curz--;
				}
			}

//			printf("edgesinfirst=%d, edgesinlast=%d\n", edgesinfirst, edgesinlast);
			if(edgesinfirst==2 && edgesinlast==2){	//we have a pattern of interest. now test if it is a local or global knot.
				patterncount++;

				//First, print out polygon which has a pattern of interest:
				fprintf(fpInterest, "%d %d %d\n", polygon[0], polygon[1], polygon[2]);
				i=3;
				while(polygon[i] != 0){
					fprintf(fpInterest, "%d\n", polygon[i]);
					i++;
				}
				fprintf(fpInterest, "-111\n");
			



				//Algorithm to find check if a pattern is local or global:
					//NEVERMIND THIS, DOES SOMETHING DIFFERENT NOW (creates two polygons, one from each arc and will knotID them)
					//get both arcs.
					//if longer arc has 22 or more edges, consider if it's a "local knot"
					//do this by using the shorter arc.
					//consider the longer arc "locally knotted" if the shorter arc has three or more edges in either the top or bottom plane (highest x or lowest x) (not including the edges in first and last section)
			
				//initialize arcs to zero
				for(i=0; i<=(L+1)*(M+1)*(span+1)-1; i++){
					arcone[i]=0;
					arctwo[i]=0;
				}
				//get both arcs.
				curx=polygon[0];
				cury=polygon[1];
				curz=polygon[2];

				//utilizing fact that the first direction of polygon will always be positive x, so that's where arcone will start
				arcone[0]=curx;
				arcone[1]=cury;
				arcone[2]=curz;

				i=3;
				while(curx!=span){	
					arcone[i]=polygon[i];
					if(polygon[i]==1){
						curx++;
					}
					else if(polygon[i]==2){
						curx--;
					}
					else if(polygon[i]==3){
						cury++;
					}
					else if(polygon[i]==4){
						cury--;
					}
					else if(polygon[i]==5){
						curz++;
					}
					else if(polygon[i]==6){
						curz--;
					}
					i++;
				}
				//at this point, arcone is filled correctly.			


				//create polygon from arcone[]
				fprintf(fptwo, "%d %d %d\n", arcone[0], arcone[1], arcone[2]);
				j=3;
				while(arcone[j] != 0){
					fprintf(fptwo, "%d\n", arcone[j]);
					j++;
				}
				//now reconnect to starting point, going "outside" of the arc
				otherx = curx;
				othery = cury;
				otherz = curz;
				//we are at (otherx, othery, otherz), need to get to (arcone[0], arcone[1], arcone[2])
				for(j=1; j<=L+1; j++){
					fprintf(fptwo, "3\n");
					othery++;
				}
				for(j=1; j<=span; j++){
					fprintf(fptwo, "2\n");
					otherx--;
				}
				if(otherx != arcone[0]){	//make sure we're at x=0
					printf("Error when walking around outside the tube. Quitting\n");
					exit(1);
				}
				while(otherz != arcone[2] || othery != arcone[1]){
					if(otherz < arcone[2]){
						fprintf(fptwo, "5\n");
						otherz++;
					}
					else if(otherz > arcone[2]){
						fprintf(fptwo, "6\n");
						otherz--;
					}
					if(othery < arcone[1]){
						fprintf(fptwo, "3\n");
						othery++;
					}
					else if(othery > arcone[1]){
						fprintf(fptwo, "4\n");
						othery--;
					}
				}
				fprintf(fptwo, "-111\n");


				//now need to find where arctwo starts. will start once the direction "2" is activated.
				while(polygon[i]!=2){
					if(polygon[i]==1){
						curx++;
					}
					else if(polygon[i]==3){
						cury++;
					}
					else if(polygon[i]==4){
						cury--;
					}
					else if(polygon[i]==5){
						curz++;
					}
					else if(polygon[i]==6){
						curz--;
					}
					i++;
				}

				//at this point, arc two will start (with the direction 2)
				arctwo[0]=curx;
				arctwo[1]=cury;
				arctwo[2]=curz;
				j=3;
				while(curx!=0){
					arctwo[j]=polygon[i];
					if(polygon[i]==1){
						curx++;
					}
					else if(polygon[i]==2){
						curx--;
					}
					else if(polygon[i]==3){
						cury++;
					}
					else if(polygon[i]==4){
						cury--;
					}
					else if(polygon[i]==5){
						curz++;
					}
					else if(polygon[i]==6){
						curz--;
					}
					i++;
					j++;
				}





				//create polygon from arctwo[]
				fprintf(fptwo, "%d %d %d\n", arctwo[0], arctwo[1], arctwo[2]);
				j=3;
				while(arctwo[j] != 0){
					fprintf(fptwo, "%d\n", arctwo[j]);
					j++;
				}
				//now reconnect to starting point, going "outside" of the arc
				otherx = curx;
				othery = cury;
				otherz = curz;
				//we are at (otherx, othery, otherz), need to get to (arctwo[0], arctwo[1], arctwo[2])
				for(j=1; j<=L+1; j++){
					fprintf(fptwo, "3\n");
					othery++;
				}
				for(j=1; j<=span; j++){
					fprintf(fptwo, "1\n");
					otherx++;
				}
				if(otherx != arctwo[0]){	//make sure we're at x=span
					printf("Error when walking around outside the tube. Quitting\n");
					exit(1);
				}
				while(otherz != arctwo[2] || othery != arctwo[1]){
					if(otherz < arctwo[2]){
						fprintf(fptwo, "5\n");
						otherz++;
					}
					else if(otherz > arctwo[2]){
						fprintf(fptwo, "6\n");
						otherz--;
					}
					if(othery < arctwo[1]){
						fprintf(fptwo, "3\n");
						othery++;
					}
					else if(othery > arctwo[1]){
						fprintf(fptwo, "4\n");
						othery--;
					}
				}

				fprintf(fptwo, "-111\n");

			}
		}





	

		printf("There are %d polygons with patterns of interest\n", patterncount);



		fprintf(fptwo, "-999");
		fprintf(fpInterest, "-999");

		fclose(fptwo);
		fclose(fp);
		fclose(fpInterest);
		printf("\nFinished file %d.\n", curfile);
	}
	printf("\nmain finished\n\n");
}



































