//This program will check to make sure the given poly file contains only valid polygons

//restriction: coords[][][] indices can't be negative. could be a problem if SAP travels into negative coordinates.

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>

int main(void){

//////////////////////////////////
	const int L=4;
	const int M=1;
	const int span=4;
	//	const int run=1;
	const int filenum=4;
//////////////////////////////////

	char filename[100];
	FILE* fp;
	char line[10];

	int coords[span+1][L+1][M+1];  //0 if empty, 1 if filled.

	int i, j, k;
	int x, y, z;
	int sx, sy, sz;	//start x, start y, start z
	int direction;
	int numpolys=0;
	int HamCounter=0;

	int curfile;

	for(curfile=1; curfile<=filenum; curfile++){
	  sprintf(filename, "ImpCreatorPrime_PKshrunkPolysL%dM%dspan%dnum%d_8_19s_shr.txt", L, M, span, curfile);
	  //		sprintf(filename, "MCpolysL%dM%dspan%drun%dnum%d.txt", L, M, span, run, curfile);
		fp = fopen(filename, "r");

		if(fp==NULL){
			printf("file could not be opened. Exitting\n");
			exit(1);
		}
		else{
			printf("opened %s.\n", filename);
		}

		fgets(line, 10, fp);
		while( fgets(line, 10, fp) != NULL ){
			for(i=0; i<=span; i++){
				for(j=0; j<=L; j++){
					for(k=0; k<=M; k++){
						coords[i][j][k]=0;
					}
				}
			}
			sscanf(line, "%d %d %d", &x, &y, &z);
			if(x==-999){
				break;
			}
			numpolys++;
//			coords[x][y][z]=1;
			sx=x;
			sy=y;
			sz=z;

			while(1){
				fgets(line, 10, fp);
				sscanf(line, "%d", &direction);
				if(direction==-111){
					break;
				}
				else if(direction==1){
					x++;
				}
				else if(direction==2){
					x--;
				}
				else if(direction==3){
					y++;
				}
				else if(direction==4){
					y--;
				}
				else if(direction==5){
					z++;
				}
				else if(direction==6){
					z--;
				}
				coords[x][y][z] += 1;
			}
			if(sx==x && sy==y && sz==z){
//				printf("woo\n");
			}
			else{
			  printf("Poly %d didn't end at start point. Exitting\n", numpolys);
			  printf("(x,y,z)=(%d,%d,%d).   (sx,sy,sz)=(%d,%d,%d)\n", x, y, z, sx, sy, sz);
				exit(1);
			}
			

			//Check coords--make sure it isn't >= 2
			int bad = 0;
			int Ham = 1;
			for(i=0; i<=span; i++){
				for(j=0; j<=L; j++){
					for(k=0; k<=M; k++){
						if(coords[i][j][k]>=2){
							bad=1;
						}
						if(coords[i][j][k]==0){
							Ham=0;
						}
//						printf("coords[%d][%d][%d]=%d\n", i, j, k, coords[i][j][k]);
					}
				}
			}

			if(Ham==1){
//				printf("Poly number %d is Hamiltonian.\n", numpolys);
				HamCounter++;
			}

			if(bad==1){
				printf("Poly number %d has a self-intersection! Exitting.\n", numpolys);
				exit(1);
			}
			printf("Poly number %d was all good.\n", numpolys);
		}
	}

	printf("There were %d polys checked. All were good!\n", numpolys);
	printf("There were %d Ham. polys.\n", HamCounter);

	return 0;
}

























