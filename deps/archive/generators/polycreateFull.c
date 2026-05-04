/* This program will create all polygons in a LxM tube with a given span. The polygon must have the
given span, no less. It may not fill the entire LxM tube, but it will go the full span that's given.
Will produce all polygons in a file of the form: FULLcreatedpolysL_M_span_.txt, where the _ change as the
system changes.	*/

//FULL VERTEX VERSION: every polygon that's created will have filled ALL of the verticies in the LxM tube.


#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>


/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDYIED **/

#define	L 2				/* number of horizontal edges  (x-direction)	*/
#define	M 1				/* number of vertical edges (y-direction)	*/
#define totalspan 3		/* span of generated polygon (z-direction)	*/

#define	vM (M+1)				/* number of vertices in vertical direction   = M + 1      */
#define	vL (L+1)				/* number of vertices in horizontal direction = L + 1      */

/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/

/***************************************************************************/
/****************** Structures used in this program ************************/
/***************************************************************************/

struct point {	//point has 3 positive integer coordinates
	unsigned short int x;
	unsigned short int y;
	unsigned short int z;
};
struct poly {	//polygon will just be a linked list of directions. (starting point will be known from global variable "startpoint")
	unsigned short int direction;	//1-6 are the directions
	struct poly *nextpoly;		//pointer to the next poly (next direction)
};



/***************************************************************************/
/****************** end of structures **************************************/
/***************************************************************************/

/***************************************************************************/
/*************** Global variables used in this program *********************/
/***************************************************************************/

unsigned short int	occupied[vL][vM][totalspan+1];	//0 indicates vertex is empty, 1 indicates occupied, 2 if occupied by startpoint
struct point		startpoint;			//point where the polygon starts
unsigned short int	maxz;				//will contain the furthest the polygon has gone in the z-direction (must be fullspan for it to be valid)
struct poly		*firstpoly;			//always points to the first "direction"
struct poly		*currentpoly;			//points to the current "direction"
unsigned long int	tot_polys=0;			//total number of polygons generated

int startwithposx;	//1 if first step is posx, 0 otherwise

/***************************************************************************/
/******************* end of global variables  ******************************/
/***************************************************************************/

/***************************************************************************/
/******************* functions used in this program  ***********************/
/***************************************************************************/ 

void goposx(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp);	//go in the positive x direction, when called, DO NOT know if this is a valid move
void gonegx(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp);	//go in the negative x direction, when called, DO NOT know if this is a valid move
void goposy(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp);	//go in the positive y direction, when called, DO NOT know if this is a valid move
void gonegy(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp);	//go in the negative y direction, when called, DO NOT know if this is a valid move
void goposz(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp);	//go in the positive z direction, when called, DO NOT know if this is a valid move
void gonegz(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp);	//go in the negative z direction, when called, DO NOT know if this is a valid move
struct poly* newpoly(void);
void printpoly(void);
void freepoly(void);	//frees everything in the linked list of poly. first and current poly point to NULL when done.
void freefirstpoly(struct poly* polypntr); //frees everything in the linked list after this pointer, and including what's pointed at





/***************************************************************************/
/******************* end of functions used in this program *****************/
/***************************************************************************/




main(void)
{
	system("date");  /* prints the date and time */
	clock();	/* Start clock to time program */

	printf("Going to create all polygons in L=%d, M=%d tube with span=%d\n", L, M, totalspan);
	unsigned long int i, j, k;	//looping variables
	FILE *fp;	//file pointer. will point to the file that is to be written to.

	char filename[30];	//this will hold the filename of the text file that we want to write to
	sprintf(filename, "FULLcreatedpolysL%dM%dspan%d.txt", L, M, totalspan);
	fp = fopen(filename, "w");	//create or overwrite the file "filename"

	if(fp != NULL){
		fprintf(fp, "UofS\n");	//first line in file is always "UofS"
		printf("printed UofS in file '%s'\n", filename);
	}
	else{
		printf("file pointer is pointing to NULL\n");
	}

	//every vertex is originally empty
	for(i=0; i<=L; i++){
		for(j=0; j<=M; j++){
			for(k=0; k<=totalspan; k++){
				occupied[i][j][k] = 0;
			}
		}
	}
//	printf("every vertex is empty\n");

			maxz=0;		//span of zero initially
			//set (0,0,0) as the starting point. since all verticies have to be filled,
			//we don't have to go through other starting points
			occupied[0][0][0] = 2;	//(0,0,0) is occupied by startpoint

			//set startpoint
			startpoint.x=0;
			startpoint.y=0;
			startpoint.z=0;
			printf("(%d,%d,%d) is the starting point\n", startpoint.x, startpoint.y, startpoint.z);

			/*Since the polygon must fill all vertices, first explore going in the positive x-direction.
			When the SAW comes back to the start, it will either come from the y-direction or the z-direction.
			To avoid recording duplicates of polygons (same polygon, just traversed in the opposite direction),
			we now only want the polygons which dont leave or comeback from the (1,0,0) point. That is, we only
			want polygons that leave and comeback from the y and z directions. WLOG, we will leave in the y direction
			and comeback from the z-direction. So if the last edge goes from (1,0,0) to (0,0,0), we don't want to record
			that polygon. (ADDED GLOBAL VARIABLE startwithposx, AND ALSO MADE CHANGES TO gonegx();

			Successfully cut the number of polygons in half, since it got rid of duplicates

			*/

			startwithposx=1;	//1 if first step is posx, 0 otherwise (global variable)
			firstpoly=newpoly();
			currentpoly=firstpoly;
			goposx(startpoint.x, startpoint.y, startpoint.z, fp);

			startwithposx=0;	//1 if first step is posx, 0 otherwise (global variable)
			freepoly();
			firstpoly=newpoly();
			currentpoly=firstpoly;
			goposy(startpoint.x, startpoint.y, startpoint.z, fp);


	fprintf(fp, "-999\n");
	printf("TOTAL NUMBER OF POLYGONS IN A L=%d, M=%d tube with span=%d is %d\n", L, M, totalspan, tot_polys);
	freepoly();	//free poly!
	fclose(fp);	//close file

//	convertUofStoPoints(filename);

	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	system("date");  /* prints the date and time */
}				/* end of main */
/***************************************************************************/

void goposx(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp){
	//go in the positive x direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevx<L){	//if there's room
		if(occupied[prevx+1][prevy][prevz]!=1){	//if it's not occupied
//			printf("can go posx from (%d,%d,%d)\n", prevx, prevy, prevz);
//			printpoly();
			currentpoly->direction = 1;	//add direction
			currentpoly->nextpoly = newpoly();
			currentpoly = currentpoly->nextpoly;
//			printpoly();

			//need to record the directions up to this point. record this in "memory"
			struct poly* firstmemory = newpoly();
			struct poly* currentmemory = firstmemory;
			currentpoly = firstpoly;

			while(currentpoly->nextpoly != NULL){
				currentmemory->direction = currentpoly->direction;

				currentpoly = currentpoly->nextpoly;
				currentmemory->nextpoly = newpoly();
				currentmemory = currentmemory->nextpoly;
			}
			//memory now holds a linked list of directions up to this point.

			if(maxz<prevz){	//check if we have a new maxz
				maxz=prevz;
			}

			//need to record maxz up to this point
			unsigned short int memoryz=maxz;
			unsigned short int i,j,k;
			if(occupied[prevx+1][prevy][prevz]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					unsigned short int fullverts=1;	//1 if every vertex is occupied, 0 if not
					for(i=0; i<=L; i++){
						for(j=0; j<=M; j++){
							for(k=0; k<=totalspan; k++){
								if(occupied[i][j][k]==0){
									fullverts=0;	//there is a vertex that isn't occupied
								}
							}
						}
					}
					if(fullverts==1){
//						printf("we have a polygon with span=%d\n",maxz);
						//fprint startingpoint and then all directions to file "filename"
						fprintf(fp, "%d %d %d\n", startpoint.x, startpoint.y, startpoint.z);
						//fprint all of the directions in currentpoly
						currentpoly=firstpoly;
						while(currentpoly->nextpoly != NULL){
							fprintf(fp, "%d\n", currentpoly->direction);
							currentpoly = currentpoly->nextpoly;
						}
						fprintf(fp, "-111\n");
						tot_polys++;
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx+1][prevy][prevz]=1;
				
				//explore going different directions (negx is not allowed now)
				goposx(prevx+1, prevy, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				goposy(prevx+1, prevy, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegy(prevx+1, prevy, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				goposz(prevx+1, prevy, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegz(prevx+1, prevy, prevz, fp);

				occupied[prevx+1][prevy][prevz]=0;	//reset back to unoccupied
			}
			freefirstpoly(firstmemory);
//			printf("finished going posx from (%d,%d,%d)\n", prevx, prevy, prevz);
		}
		else{
//			printf("can't go posx from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx+1, prevy, prevz);
		}
	}
	else{
//		printf("can't go posx from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void gonegx(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp){
	//go in the negative x direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevx>0){	//if there's room
		if((occupied[prevx-1][prevy][prevz]!=1) && (!(startwithposx==0 && (prevx-1)==0 && prevy==0 && prevz==0))){	//if it's not occupied  AND  !(didnt start with posx AND next vertex would be origin)
//			printf("can go negx from (%d,%d,%d)\n", prevx, prevy, prevz);
//			printpoly();
			currentpoly->direction = 2;	//add direction
			currentpoly->nextpoly = newpoly();
			currentpoly = currentpoly->nextpoly;
//			printpoly();

			//need to record the directions up to this point. record this in "memory"
			struct poly* firstmemory = newpoly();
			struct poly* currentmemory = firstmemory;
			currentpoly = firstpoly;

			while(currentpoly->nextpoly != NULL){
				currentmemory->direction = currentpoly->direction;

				currentpoly = currentpoly->nextpoly;
				currentmemory->nextpoly = newpoly();
				currentmemory = currentmemory->nextpoly;
			}
			//memory now holds a linked list of directions up to this point.

			if(maxz<prevz){	//check if we have a new maxz
				maxz=prevz;
			}

			//need to record maxz up to this point
			unsigned short int memoryz=maxz;
			unsigned short int i,j,k;
			if(occupied[prevx-1][prevy][prevz]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					unsigned short int fullverts=1;	//1 if every vertex is occupied, 0 if not
					for(i=0; i<=L; i++){
						for(j=0; j<=M; j++){
							for(k=0; k<=totalspan; k++){
								if(occupied[i][j][k]==0){
									fullverts=0;	//there is a vertex that isn't occupied
								}
							}
						}
					}
					if(fullverts==1){
//						printf("we have a polygon with span=%d\n",maxz);
						//fprint startingpoint and then all directions to file "filename"
						fprintf(fp, "%d %d %d\n", startpoint.x, startpoint.y, startpoint.z);
						//fprint all of the directions in currentpoly
						currentpoly=firstpoly;
						while(currentpoly->nextpoly != NULL){
							fprintf(fp, "%d\n", currentpoly->direction);
							currentpoly = currentpoly->nextpoly;
						}
						fprintf(fp, "-111\n");
						tot_polys++;
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx-1][prevy][prevz]=1;
				
				//explore going different directions (posx is not allowed now)
				gonegx(prevx-1, prevy, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				goposy(prevx-1, prevy, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegy(prevx-1, prevy, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				goposz(prevx-1, prevy, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegz(prevx-1, prevy, prevz, fp);

				occupied[prevx-1][prevy][prevz]=0;	//reset back to unoccupied
			}
			freefirstpoly(firstmemory);
//			printf("finished going negx from (%d,%d,%d)\n", prevx, prevy, prevz);
		}
		else{
//			printf("can't go negx from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx-1, prevy, prevz);
		}
	}
	else{
//		printf("can't go negx from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void goposy(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp){
	//go in the positive y direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevy<M){	//if there's room
		if(occupied[prevx][prevy+1][prevz]!=1){	//if it's not occupied
//			printf("can go posy from (%d,%d,%d)\n", prevx, prevy, prevz);
//			printpoly();
			currentpoly->direction = 3;	//add direction
			currentpoly->nextpoly = newpoly();
			currentpoly = currentpoly->nextpoly;
//			printpoly();

			//need to record the directions up to this point. record this in "memory"
			struct poly* firstmemory = newpoly();
			struct poly* currentmemory = firstmemory;
			currentpoly = firstpoly;

			while(currentpoly->nextpoly != NULL){
				currentmemory->direction = currentpoly->direction;

				currentpoly = currentpoly->nextpoly;
				currentmemory->nextpoly = newpoly();
				currentmemory = currentmemory->nextpoly;
			}
			//memory now holds a linked list of directions up to this point.

			if(maxz<prevz){	//check if we have a new maxz
				maxz=prevz;
			}

			//need to record maxz up to this point
			unsigned short int memoryz=maxz;
			unsigned short int i,j,k;
			if(occupied[prevx][prevy+1][prevz]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					unsigned short int fullverts=1;	//1 if every vertex is occupied, 0 if not
					for(i=0; i<=L; i++){
						for(j=0; j<=M; j++){
							for(k=0; k<=totalspan; k++){
								if(occupied[i][j][k]==0){
									fullverts=0;	//there is a vertex that isn't occupied
								}
							}
						}
					}
					if(fullverts==1){
//						printf("we have a polygon with span=%d\n",maxz);
						//fprint startingpoint and then all directions to file "filename"
						fprintf(fp, "%d %d %d\n", startpoint.x, startpoint.y, startpoint.z);
						//fprint all of the directions in currentpoly
						currentpoly=firstpoly;
						while(currentpoly->nextpoly != NULL){
							fprintf(fp, "%d\n", currentpoly->direction);
							currentpoly = currentpoly->nextpoly;
						}
						fprintf(fp, "-111\n");
						tot_polys++;
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx][prevy+1][prevz]=1;
				
				//explore going different directions (negy is not allowed now)
				goposx(prevx, prevy+1, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegx(prevx, prevy+1, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				goposy(prevx, prevy+1, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				goposz(prevx, prevy+1, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegz(prevx, prevy+1, prevz, fp);

				occupied[prevx][prevy+1][prevz]=0;	//reset back to unoccupied
			}
			freefirstpoly(firstmemory);
//			printf("finished going posy from (%d,%d,%d)\n", prevx, prevy, prevz);
		}
		else{
//			printf("can't go posy from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx, prevy+1, prevz);
		}
	}
	else{
//		printf("can't go posy from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void gonegy(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp){
	//go in the negative y direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevy>0){	//if there's room
		if(occupied[prevx][prevy-1][prevz]!=1){	//if it's not occupied
//			printf("can go negy from (%d,%d,%d)\n", prevx, prevy, prevz);
//			printpoly();
			currentpoly->direction = 4;	//add direction
			currentpoly->nextpoly = newpoly();
			currentpoly = currentpoly->nextpoly;
//			printpoly();

			//need to record the directions up to this point. record this in "memory"
			struct poly* firstmemory = newpoly();
			struct poly* currentmemory = firstmemory;
			currentpoly = firstpoly;

			while(currentpoly->nextpoly != NULL){
				currentmemory->direction = currentpoly->direction;

				currentpoly = currentpoly->nextpoly;
				currentmemory->nextpoly = newpoly();
				currentmemory = currentmemory->nextpoly;
			}
			//memory now holds a linked list of directions up to this point.

			if(maxz<prevz){	//check if we have a new maxz
				maxz=prevz;
			}

			//need to record maxz up to this point
			unsigned short int memoryz=maxz;
			unsigned short int i,j,k;
			if(occupied[prevx][prevy-1][prevz]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					unsigned short int fullverts=1;	//1 if every vertex is occupied, 0 if not
					for(i=0; i<=L; i++){
						for(j=0; j<=M; j++){
							for(k=0; k<=totalspan; k++){
								if(occupied[i][j][k]==0){
									fullverts=0;	//there is a vertex that isn't occupied
								}
							}
						}
					}
					if(fullverts==1){
//						printf("we have a polygon with span=%d\n",maxz);
						//fprint startingpoint and then all directions to file "filename"
						fprintf(fp, "%d %d %d\n", startpoint.x, startpoint.y, startpoint.z);
						//fprint all of the directions in currentpoly
						currentpoly=firstpoly;
						while(currentpoly->nextpoly != NULL){
							fprintf(fp, "%d\n", currentpoly->direction);
							currentpoly = currentpoly->nextpoly;
						}
						fprintf(fp, "-111\n");
						tot_polys++;
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx][prevy-1][prevz]=1;
				
				//explore going different directions (posy is not allowed now)
				goposx(prevx, prevy-1, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegx(prevx, prevy-1, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegy(prevx, prevy-1, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				goposz(prevx, prevy-1, prevz, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegz(prevx, prevy-1, prevz, fp);

				occupied[prevx][prevy-1][prevz]=0;	//reset back to unoccupied
			}
			freefirstpoly(firstmemory);
//			printf("finished going negy from (%d,%d,%d)\n", prevx, prevy, prevz);
		}
		else{
//			printf("can't go negy from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx, prevy-1, prevz);
		}
	}
	else{
//		printf("can't go negy from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void goposz(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp){
	//go in the positive z direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevz<totalspan){	//if there's room
		if(occupied[prevx][prevy][prevz+1]!=1){	//if it's not occupied
//			printf("can go posz from (%d,%d,%d)\n", prevx, prevy, prevz);
//			printpoly();
			currentpoly->direction = 5;	//add direction
			currentpoly->nextpoly = newpoly();
			currentpoly = currentpoly->nextpoly;
//			printpoly();

			//need to record the directions up to this point. record this in "memory"
			struct poly* firstmemory = newpoly();
			struct poly* currentmemory = firstmemory;
			currentpoly = firstpoly;

			while(currentpoly->nextpoly != NULL){
				currentmemory->direction = currentpoly->direction;

				currentpoly = currentpoly->nextpoly;
				currentmemory->nextpoly = newpoly();
				currentmemory = currentmemory->nextpoly;
			}
			//memory now holds a linked list of directions up to this point.

			if(maxz<prevz+1){	//check if we have a new maxz
				maxz=prevz+1;
			}

			//need to record maxz up to this point
			unsigned short int memoryz=maxz;
			unsigned short int i,j,k;
			if(occupied[prevx][prevy][prevz+1]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					unsigned short int fullverts=1;	//1 if every vertex is occupied, 0 if not
					for(i=0; i<=L; i++){
						for(j=0; j<=M; j++){
							for(k=0; k<=totalspan; k++){
								if(occupied[i][j][k]==0){
									fullverts=0;	//there is a vertex that isn't occupied
								}
							}
						}
					}
					if(fullverts==1){
//						printf("we have a polygon with span=%d\n",maxz);
						//fprint startingpoint and then all directions to file "filename"
						fprintf(fp, "%d %d %d\n", startpoint.x, startpoint.y, startpoint.z);
						//fprint all of the directions in currentpoly
						currentpoly=firstpoly;
						while(currentpoly->nextpoly != NULL){
							fprintf(fp, "%d\n", currentpoly->direction);
							currentpoly = currentpoly->nextpoly;
						}
						fprintf(fp, "-111\n");
						tot_polys++;
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx][prevy][prevz+1]=1;
				
				//explore going different directions (negz is not allowed now)
				goposx(prevx, prevy, prevz+1, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegx(prevx, prevy, prevz+1, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				goposy(prevx, prevy, prevz+1, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegy(prevx, prevy, prevz+1, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				goposz(prevx, prevy, prevz+1, fp);

				occupied[prevx][prevy][prevz+1]=0;	//reset back to unoccupied
			}
			freefirstpoly(firstmemory);
//			printf("finished going posz from (%d,%d,%d)\n", prevx, prevy, prevz);
		}
		else{
//			printf("can't go posz from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx, prevy, prevz+1);
		}
	}
	else{
//		printf("can't go posz from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void gonegz(unsigned short int prevx, unsigned short int prevy, unsigned short int prevz, FILE* fp){
	//go in the negative z direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevz>0){	//if there's room
		if(occupied[prevx][prevy][prevz-1]!=1){	//if it's not occupied
//			printf("can go negz from (%d,%d,%d)\n", prevx, prevy, prevz);
//			printpoly();
			currentpoly->direction = 6;	//add direction
			currentpoly->nextpoly = newpoly();
			currentpoly = currentpoly->nextpoly;
//			printpoly();

			//need to record the directions up to this point. record this in "memory"
			struct poly* firstmemory = newpoly();
			struct poly* currentmemory = firstmemory;
			currentpoly = firstpoly;

			while(currentpoly->nextpoly != NULL){
				currentmemory->direction = currentpoly->direction;

				currentpoly = currentpoly->nextpoly;
				currentmemory->nextpoly = newpoly();
				currentmemory = currentmemory->nextpoly;
			}
			//memory now holds a linked list of directions up to this point.

			if(maxz<prevz){	//check if we have a new maxz
				maxz=prevz;
			}

			//need to record maxz up to this point
			unsigned short int memoryz=maxz;
			unsigned short int i,j,k;
			if(occupied[prevx][prevy][prevz-1]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					unsigned short int fullverts=1;	//1 if every vertex is occupied, 0 if not
					for(i=0; i<=L; i++){
						for(j=0; j<=M; j++){
							for(k=0; k<=totalspan; k++){
								if(occupied[i][j][k]==0){
									fullverts=0;	//there is a vertex that isn't occupied
								}
							}
						}
					}
					if(fullverts==1){
//						printf("we have a polygon with span=%d\n",maxz);
						//fprint startingpoint and then all directions to file "filename"
						fprintf(fp, "%d %d %d\n", startpoint.x, startpoint.y, startpoint.z);
						//fprint all of the directions in currentpoly
						currentpoly=firstpoly;
						while(currentpoly->nextpoly != NULL){
							fprintf(fp, "%d\n", currentpoly->direction);
							currentpoly = currentpoly->nextpoly;
						}
						fprintf(fp, "-111\n");
						tot_polys++;
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx][prevy][prevz-1]=1;
				
				//explore going different directions (posz is not allowed now)
				goposx(prevx, prevy, prevz-1, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegx(prevx, prevy, prevz-1, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				goposy(prevx, prevy, prevz-1, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegy(prevx, prevy, prevz-1, fp);

				//copy memory into poly
				freepoly();
				firstpoly=newpoly();
				currentpoly = firstpoly;			int startwithposx=1;	//1 if first step is posx, 0 otherwise
				currentmemory = firstmemory;
				while(currentmemory->nextpoly != NULL){
					currentpoly->direction = currentmemory->direction;
					currentmemory = currentmemory->nextpoly;
					currentpoly->nextpoly = newpoly();
					currentpoly = currentpoly->nextpoly;
				}
				maxz=memoryz;
				gonegz(prevx, prevy, prevz-1, fp);

				occupied[prevx][prevy][prevz-1]=0;	//reset back to unoccupied
			}
			freefirstpoly(firstmemory);
//			printf("finished going negz from (%d,%d,%d)\n", prevx, prevy, prevz);
		}
		else{
//			printf("can't go negz from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx, prevy, prevz-1);
		}
	}
	else{
//		printf("can't go negz from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}

struct poly* newpoly(void)		/* returns a pointer to a poly structure */
/* this is a memory allocation function */
{
	struct poly* nextnewpoly;

	nextnewpoly = ((struct poly *) calloc(1, sizeof(struct poly)));
	if (nextnewpoly == NULL) {
		fprintf(stderr, "unable to allocate memory");
		exit(1);
	}
	return nextnewpoly;
}

void printpoly(void){
	//prints whatever's in poly at the time
	currentpoly=firstpoly;
	printf("printing poly: ");
	while(currentpoly->nextpoly != NULL){
		printf("%d, ", currentpoly->direction);
		currentpoly = currentpoly->nextpoly;
	}
	printf("\n");
}

void freepoly(void){
	//frees everything in the poly linked list
	struct poly* polytofree = firstpoly;
	currentpoly = firstpoly;
	firstpoly=NULL;
	while(currentpoly != NULL){
		currentpoly = currentpoly->nextpoly;
		free(polytofree);
		polytofree=currentpoly;
	}
}

void freefirstpoly(struct poly* polypntr){
	//frees everything in the linked list after this pointer, and including what's pointed at
	struct poly* polytofree = polypntr;
	while(polypntr != NULL){
		polypntr = polypntr->nextpoly;
		free(polytofree);
		polytofree = polypntr;
	}
}






















