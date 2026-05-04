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
#define	totalspan 4		/* span of generated polygon (z-direction)	*/

#define	vM (M+1)			/* number of vertices in vertical direction   = M + 1      */
#define	vL (L+1)			/* number of vertices in horizontal direction = L + 1      */

#define 	maxpolys 900000	/* maximum number of polygons in each .txt file	*/

/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/

/***************************************************************************/
/****************** Structures used in this program ************************/
/***************************************************************************/


struct point {	//point has 3 positive integer coordinates
	int x;
	int y;
	int z;
};


/***************************************************************************/
/****************** end of structures **************************************/
/***************************************************************************/

/***************************************************************************/
/*************** Global variables used in this program *********************/
/***************************************************************************/

int			occupied[vL][vM][totalspan+1];	//0 indicates vertex is empty, 1 indicates occupied, 2 if occupied by startpoint
struct point	startpoint;			//point where the polygon starts
int			maxz;				//will contain the furthest the polygon has gone in the z-direction (must be fullspan for it to be valid)
int			tot_polys=0;			//total number of polygons generated
int 			filenum;				//used to rename the file, when there are too many polygons
char 		filename[36];			//this will hold the filename of the text file that we want to write to

int			startwithposx;			//1 if first step is posx, 0 otherwise

int			thepolygon[vM*vL*(totalspan+1)];	//this global array will hold the polygon (list of directions). maximum number of directions will be vM*vL
										//will be initialized to all zeros once main starts.
int			curious=0;


/***************************************************************************/
/******************* end of global variables  ******************************/
/***************************************************************************/

/***************************************************************************/
/******************* functions used in this program  ***********************/
/***************************************************************************/ 

void goposx(int prevx, int prevy, int prevz, FILE* fp);	//go in the positive x direction, when called, DO NOT know if this is a valid move when called
void gonegx(int prevx, int prevy, int prevz, FILE* fp);	//go in the negative x direction, when called, DO NOT know if this is a valid move when called
void goposy(int prevx, int prevy, int prevz, FILE* fp);	//go in the positive y direction, when called, DO NOT know if this is a valid move when called
void gonegy(int prevx, int prevy, int prevz, FILE* fp);	//go in the negative y direction, when called, DO NOT know if this is a valid move when called
void goposz(int prevx, int prevy, int prevz, FILE* fp);	//go in the positive z direction, when called, DO NOT know if this is a valid move when called
void gonegz(int prevx, int prevy, int prevz, FILE* fp);	//go in the negative z direction, when called, DO NOT know if this is a valid move when called

void printthepolygon();	//prints the set of directions that are current in thepolygon[]


/***************************************************************************/
/******************* end of functions used in this program *****************/
/***************************************************************************/




main(void)
{
	system("date");  /* prints the date and time */
	clock();	/* Start clock to time program */

	printf("Going to create all polygons in L=%d, M=%d tube with span=%d\n", L, M, totalspan);
	int i, j, k;	//looping variables
	FILE *fp;	//file pointer. will point to the file that is to be written to.

	filenum=1;	//initialize to 1

	sprintf(filename, "FULLcreatedpolysL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);
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

	//initialize thepolygon[] to be all zeros
	for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
		thepolygon[i]=0;
	}

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
	goposx(startpoint.x, startpoint.y, startpoint.z, fp);

	startwithposx=0;	//1 if first step is posx, 0 otherwise (global variable)
	goposy(startpoint.x, startpoint.y, startpoint.z, fp);


	fprintf(fp, "-999\n");
	printf("TOTAL NUMBER OF POLYGONS IN A L=%d, M=%d tube with span=%d is %d\n", L, M, totalspan, tot_polys);

	fclose(fp);	//close file

//	convertUofStoPoints(filename);

	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	printf("curious=%d\n",curious);
	system("date");  /* prints the date and time */
}				/* end of main */
/***************************************************************************/

void goposx(int prevx, int prevy, int prevz, FILE* fp){
	int i, j, k;
	//go in the positive x direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevx<L){	//if there's room
		if(occupied[prevx+1][prevy][prevz]!=1){	//if it's not occupied
//			printf("can go posx from (%d,%d,%d)\n", prevx, prevy, prevz);

			//add direction
			i=0;
			while(thepolygon[i]!=0){		//go to end of polygon
				i++;
			}
			thepolygon[i]=1;
			
//			printthepolygon();

			if(maxz<prevz){	//check if we have a new maxz
				maxz=prevz;
			}

			//need to record maxz up to this point
			int memoryz=maxz;
			if(occupied[prevx+1][prevy][prevz]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					int fullverts=1;	//1 if every vertex is occupied, 0 if not
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
						//fprint all of the directions in thepolygon[]
						for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
							fprintf(fp, "%d\n", thepolygon[i]);
						}
						fprintf(fp, "-111\n");
						tot_polys++;
						if(thepolygon[vM*vL*(totalspan+1)-1]==6){
							curious++;
						}



						//check if we need to create another file (cause too many polys)
						if(tot_polys>(filenum*maxpolys)){	//then we need to close up the current file, and create a new one to write to.
							fprintf(fp, "-999\n");
							fclose(fp);
							filenum++;
							sprintf(filename, "FULLcreatedpolysL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);
							fp = fopen(filename, "w");	//create or overwrite the file "filename"

							if(fp != NULL){
								fprintf(fp, "UofS\n");	//first line in file is always "UofS"
								printf("printed UofS in file '%s' (tot_polys=%d)\n", filename, tot_polys);
							}
							else{
								printf("file pointer is pointing to NULL\n");
							}
						}
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx+1][prevy][prevz]=1;
				
				//explore going different directions (negx is not allowed now)
				goposx(prevx+1, prevy, prevz, fp);

				maxz=memoryz;
				goposy(prevx+1, prevy, prevz, fp);

				maxz=memoryz;
				gonegy(prevx+1, prevy, prevz, fp);

				maxz=memoryz;
				goposz(prevx+1, prevy, prevz, fp);

				maxz=memoryz;
				gonegz(prevx+1, prevy, prevz, fp);

				occupied[prevx+1][prevy][prevz]=0;	//reset back to unoccupied
			}
//			printf("finished going posx from (%d,%d,%d)\n", prevx, prevy, prevz);
			//erase the last direction in thepolygon[]
			int i=0;
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				if(thepolygon[i+1]==0){
					thepolygon[i]=0;
				}
			}
			thepolygon[vM*vL*(totalspan+1)-1]=0;


//			printthepolygon();
		}
		else{
//			printf("can't go posx from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx+1, prevy, prevz);
		}
	}
	else{
//		printf("can't go posx from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void gonegx(int prevx, int prevy, int prevz, FILE* fp){
	int i, j, k;
	//go in the negative x direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevx>0){	//if there's room
		if((occupied[prevx-1][prevy][prevz]!=1) && (!(startwithposx==0 && (prevx-1)==0 && prevy==0 && prevz==0))){	//if it's not occupied  AND  !(didnt start with posx AND next vertex would be origin)
//			printf("can go negx from (%d,%d,%d)\n", prevx, prevy, prevz);

			//add direction
			i=0;
			while(thepolygon[i]!=0){		//go to end of polygon
				i++;
			}
			thepolygon[i]=2;
			
//			printthepolygon();

			if(maxz<prevz){	//check if we have a new maxz
				maxz=prevz;
			}

			//need to record maxz up to this point
			int memoryz=maxz;
			if(occupied[prevx-1][prevy][prevz]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					int fullverts=1;	//1 if every vertex is occupied, 0 if not
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
						//fprint all of the directions in thepolygon[]
						for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
							fprintf(fp, "%d\n", thepolygon[i]);
						}
						fprintf(fp, "-111\n");
						tot_polys++;
						if(thepolygon[vM*vL*(totalspan+1)-1]==6){
							curious++;
						}

						//check if we need to create another file (cause too many polys)
						if(tot_polys>(filenum*maxpolys)){	//then we need to close up the current file, and create a new one to write to.
							fprintf(fp, "-999\n");
							fclose(fp);
							filenum++;
							sprintf(filename, "FULLcreatedpolysL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);
							fp = fopen(filename, "w");	//create or overwrite the file "filename"

							if(fp != NULL){
								fprintf(fp, "UofS\n");	//first line in file is always "UofS"
								printf("printed UofS in file '%s' (tot_polys=%d)\n", filename, tot_polys);
							}
							else{
								printf("file pointer is pointing to NULL\n");
							}
						}
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx-1][prevy][prevz]=1;
				
				//explore going different directions (posx is not allowed now)
				gonegx(prevx-1, prevy, prevz, fp);

				maxz=memoryz;
				goposy(prevx-1, prevy, prevz, fp);

				maxz=memoryz;
				gonegy(prevx-1, prevy, prevz, fp);

				maxz=memoryz;
				goposz(prevx-1, prevy, prevz, fp);

				maxz=memoryz;
				gonegz(prevx-1, prevy, prevz, fp);

				occupied[prevx-1][prevy][prevz]=0;	//reset back to unoccupied
			}
//			printf("finished going negx from (%d,%d,%d)\n", prevx, prevy, prevz);
			//erase the last direction in thepolygon[]
			int i=0;
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				if(thepolygon[i+1]==0){
					thepolygon[i]=0;
				}
			}
			thepolygon[vM*vL*(totalspan+1)-1]=0;


//			printthepolygon();
		}
		else{
//			printf("can't go negx from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx-1, prevy, prevz);
		}
	}
	else{
//		printf("can't go negx from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void goposy(int prevx, int prevy, int prevz, FILE* fp){
	int i, j, k;
	//go in the positive y direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevy<M){	//if there's room
		if(occupied[prevx][prevy+1][prevz]!=1){	//if it's not occupied
//			printf("can go posy from (%d,%d,%d)\n", prevx, prevy, prevz);

			//add direction
			i=0;
			while(thepolygon[i]!=0){		//go to end of polygon
				i++;
			}
			thepolygon[i]=3;
			
//			printthepolygon();

			if(maxz<prevz){	//check if we have a new maxz
				maxz=prevz;
			}

			//need to record maxz up to this point
			int memoryz=maxz;
			if(occupied[prevx][prevy+1][prevz]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					int fullverts=1;	//1 if every vertex is occupied, 0 if not
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
						//fprint all of the directions in thepolygon[]
						for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
							fprintf(fp, "%d\n", thepolygon[i]);
						}
						fprintf(fp, "-111\n");
						tot_polys++;
						if(thepolygon[vM*vL*(totalspan+1)-1]==6){
							curious++;
						}


						//check if we need to create another file (cause too many polys)
						if(tot_polys>(filenum*maxpolys)){	//then we need to close up the current file, and create a new one to write to.
							fprintf(fp, "-999\n");
							fclose(fp);
							filenum++;
							sprintf(filename, "FULLcreatedpolysL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);
							fp = fopen(filename, "w");	//create or overwrite the file "filename"

							if(fp != NULL){
								fprintf(fp, "UofS\n");	//first line in file is always "UofS"
								printf("printed UofS in file '%s' (tot_polys=%d)\n", filename, tot_polys);
							}
							else{
								printf("file pointer is pointing to NULL\n");
							}
						}
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx][prevy+1][prevz]=1;
				
				//explore going different directions (negy is not allowed now)
				goposx(prevx, prevy+1, prevz, fp);

				maxz=memoryz;
				gonegx(prevx, prevy+1, prevz, fp);

				maxz=memoryz;
				goposy(prevx, prevy+1, prevz, fp);

				maxz=memoryz;
				goposz(prevx, prevy+1, prevz, fp);

				maxz=memoryz;
				gonegz(prevx, prevy+1, prevz, fp);

				occupied[prevx][prevy+1][prevz]=0;	//reset back to unoccupied
			}
//			printf("finished going posy from (%d,%d,%d)\n", prevx, prevy, prevz);
			//erase the last direction in thepolygon[]
			int i=0;
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				if(thepolygon[i+1]==0){
					thepolygon[i]=0;
				}
			}
			thepolygon[vM*vL*(totalspan+1)-1]=0;


//			printthepolygon();
		}
		else{
//			printf("can't go posy from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx, prevy+1, prevz);
		}
	}
	else{
//		printf("can't go posy from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void gonegy(int prevx, int prevy, int prevz, FILE* fp){
	int i, j, k;
	//go in the negative y direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevy>0){	//if there's room
		if(occupied[prevx][prevy-1][prevz]!=1){	//if it's not occupied
//			printf("can go negy from (%d,%d,%d)\n", prevx, prevy, prevz);

			//add direction
			i=0;
			while(thepolygon[i]!=0){		//go to end of polygon
				i++;
			}
			thepolygon[i]=4;
			
//			printthepolygon();

			if(maxz<prevz){	//check if we have a new maxz
				maxz=prevz;
			}

			//need to record maxz up to this point
			int memoryz=maxz;
			if(occupied[prevx][prevy-1][prevz]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					int fullverts=1;	//1 if every vertex is occupied, 0 if not
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
						//fprint all of the directions in thepolygon[]
						for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
							fprintf(fp, "%d\n", thepolygon[i]);
						}
						fprintf(fp, "-111\n");
						tot_polys++;
						if(thepolygon[vM*vL*(totalspan+1)-1]==6){
							curious++;
						}


						//check if we need to create another file (cause too many polys)
						if(tot_polys>(filenum*maxpolys)){	//then we need to close up the current file, and create a new one to write to.
							fprintf(fp, "-999\n");
							fclose(fp);
							filenum++;
							sprintf(filename, "FULLcreatedpolysL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);
							fp = fopen(filename, "w");	//create or overwrite the file "filename"

							if(fp != NULL){
								fprintf(fp, "UofS\n");	//first line in file is always "UofS"
								printf("printed UofS in file '%s' (tot_polys=%d)\n", filename, tot_polys);
							}
							else{
								printf("file pointer is pointing to NULL\n");
							}
						}
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx][prevy-1][prevz]=1;
				
				//explore going different directions (posy is not allowed now)
				goposx(prevx, prevy-1, prevz, fp);

				maxz=memoryz;
				gonegx(prevx, prevy-1, prevz, fp);

				maxz=memoryz;
				gonegy(prevx, prevy-1, prevz, fp);

				maxz=memoryz;
				goposz(prevx, prevy-1, prevz, fp);

				maxz=memoryz;
				gonegz(prevx, prevy-1, prevz, fp);

				occupied[prevx][prevy-1][prevz]=0;	//reset back to unoccupied
			}
//			printf("finished going negy from (%d,%d,%d)\n", prevx, prevy, prevz);
			//erase the last direction in thepolygon[]
			int i=0;
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				if(thepolygon[i+1]==0){
					thepolygon[i]=0;
				}
			}
			thepolygon[vM*vL*(totalspan+1)-1]=0;


//			printthepolygon();
		}
		else{
//			printf("can't go negy from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx, prevy-1, prevz);
		}
	}
	else{
//		printf("can't go negy from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void goposz(int prevx, int prevy, int prevz, FILE* fp){
	int i, j, k;
	//go in the positive z direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevz<totalspan){	//if there's room
		if(occupied[prevx][prevy][prevz+1]!=1){	//if it's not occupied
//			printf("can go posz from (%d,%d,%d)\n", prevx, prevy, prevz);

			//add direction
			i=0;
			while(thepolygon[i]!=0){		//go to end of polygon
				i++;
			}
			thepolygon[i]=5;
			
//			printthepolygon();

			if(maxz<prevz+1){	//check if we have a new maxz
				maxz=prevz+1;
			}

			//need to record maxz up to this point
			int memoryz=maxz;
			if(occupied[prevx][prevy][prevz+1]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					int fullverts=1;	//1 if every vertex is occupied, 0 if not
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
						//fprint all of the directions in thepolygon[]
						for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
							fprintf(fp, "%d\n", thepolygon[i]);
						}
						fprintf(fp, "-111\n");
						tot_polys++;
						if(thepolygon[vM*vL*(totalspan+1)-1]==6){
							curious++;
						}


						//check if we need to create another file (cause too many polys)
						if(tot_polys>(filenum*maxpolys)){	//then we need to close up the current file, and create a new one to write to.
							fprintf(fp, "-999\n");
							fclose(fp);
							filenum++;
							sprintf(filename, "FULLcreatedpolysL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);
							fp = fopen(filename, "w");	//create or overwrite the file "filename"

							if(fp != NULL){
								fprintf(fp, "UofS\n");	//first line in file is always "UofS"
								printf("printed UofS in file '%s' (tot_polys=%d)\n", filename, tot_polys);
							}
							else{
								printf("file pointer is pointing to NULL\n");
							}
						}
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx][prevy][prevz+1]=1;
				
				//explore going different directions (negz is not allowed now)
				goposx(prevx, prevy, prevz+1, fp);

				maxz=memoryz;
				gonegx(prevx, prevy, prevz+1, fp);

				maxz=memoryz;
				goposy(prevx, prevy, prevz+1, fp);

				maxz=memoryz;
				gonegy(prevx, prevy, prevz+1, fp);

				maxz=memoryz;
				goposz(prevx, prevy, prevz+1, fp);

				occupied[prevx][prevy][prevz+1]=0;	//reset back to unoccupied
			}
//			printf("finished going posz from (%d,%d,%d)\n", prevx, prevy, prevz);
			//erase the last direction in thepolygon[]
			int i=0;
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				if(thepolygon[i+1]==0){
					thepolygon[i]=0;
				}
			}
			thepolygon[vM*vL*(totalspan+1)-1]=0;

//			printthepolygon();
		}
		else{
//			printf("can't go posz from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx, prevy, prevz+1);
		}
	}
	else{
//		printf("can't go posz from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void gonegz(int prevx, int prevy, int prevz, FILE* fp){

	int i, j, k;
	//go in the negative z direction, when called, DO NOT know if this is a valid move
	//test if this is a valid move
	if(prevz>0){	//if there's room
		if(occupied[prevx][prevy][prevz-1]!=1){	//if it's not occupied
//			printf("can go negz from (%d,%d,%d)\n", prevx, prevy, prevz);

			//add direction
			i=0;
			while(thepolygon[i]!=0){		//go to end of polygon
				i++;
			}
			thepolygon[i]=6;
			
//			printthepolygon();

			if(maxz<prevz){	//check if we have a new maxz
				maxz=prevz;
			}

			//need to record maxz up to this point
			int memoryz=maxz;
			if(occupied[prevx][prevy][prevz-1]==2){	//if we are now at the startingpoint
//				printf("we are back at the starting point\n");
				if(maxz==totalspan){
					//need to check if all verticies are occupied
					int fullverts=1;	//1 if every vertex is occupied, 0 if not
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
						//fprint all of the directions in thepolygon[]
						for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
							fprintf(fp, "%d\n", thepolygon[i]);
						}
						fprintf(fp, "-111\n");
						tot_polys++;
						if(thepolygon[vM*vL*(totalspan+1)-1]==6){
							curious++;
						}


						//check if we need to create another file (cause too many polys)
						if(tot_polys>(filenum*maxpolys)){	//then we need to close up the current file, and create a new one to write to.
							fprintf(fp, "-999\n");
							fclose(fp);
							filenum++;
							sprintf(filename, "FULLcreatedpolysL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);
							fp = fopen(filename, "w");	//create or overwrite the file "filename"

							if(fp != NULL){
								fprintf(fp, "UofS\n");	//first line in file is always "UofS"
								printf("printed UofS in file '%s' (tot_polys=%d)\n", filename, tot_polys);
							}
							else{
								printf("file pointer is pointing to NULL\n");
							}
						}
					}
				}
			}
			else{	//not at the starting point
//				printf("not at starting point yet\n");
				occupied[prevx][prevy][prevz-1]=1;
				
				//explore going different directions (posz is not allowed now)
				goposx(prevx, prevy, prevz-1, fp);

				maxz=memoryz;
				gonegx(prevx, prevy, prevz-1, fp);

				maxz=memoryz;
				goposy(prevx, prevy, prevz-1, fp);

				maxz=memoryz;
				gonegy(prevx, prevy, prevz-1, fp);

				maxz=memoryz;
				gonegz(prevx, prevy, prevz-1, fp);

				occupied[prevx][prevy][prevz-1]=0;	//reset back to unoccupied
			}
//			printf("finished going negz from (%d,%d,%d)\n", prevx, prevy, prevz);
			//erase the last direction in thepolygon[]
			int i=0;
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				if(thepolygon[i+1]==0){
					thepolygon[i]=0;
				}
			}
			thepolygon[vM*vL*(totalspan+1)-1]=0;

//			printthepolygon();
		}
		else{
//			printf("can't go negz from (%d,%d,%d) because (%d,%d,%d) is already occupied\n", prevx, prevy, prevz, prevx, prevy, prevz-1);
		}
	}
	else{
//		printf("can't go negz from (%d,%d,%d) because there's no room\n", prevx, prevy, prevz);
	}
}


void printthepolygon(){
	printf("Printing thepolygon: ");
	int i=0;
	for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
		printf("%d, ", thepolygon[i]);
	}
	printf("\n");
}



















