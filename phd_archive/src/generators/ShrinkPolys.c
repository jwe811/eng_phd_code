//This program will take a polygon file and shrink the polygons (via bfacf_shrink_v2.c) . Will then be outputted with suffix "_shr"

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>



	#define L 4
	#define M 1
	#define span 4
//	#define run 1
	#define filenum 4	//total number of files

int bfacf_shrink(int x, int y, int z, int polygon[], int length);

int main(void){
	char filename[100];	//file with polys
	char filename2[100];	//out file that will have shrunk polys
	FILE* fp;
	FILE* fp2;
	char line[10];


	static int Polygon[(L+1)*(M+1)*(span+1)];

	int i;
	int x, y, z;
	int index;
	int direction;
	int numpolys=0;
	int curfile;

	for(curfile=1; curfile<=filenum; curfile++){
		//file that contains polys
	  //		sprintf(filename, "MCpolysHamL%dM%dspan%drun%dnum%d.txt", L, M, span, run, curfile);
		sprintf(filename, "ImpCreatorPrime_PKshrunkPolysL%dM%dspan%dnum%d_8_19s.txt", L, M, span, curfile);
		fp = fopen(filename, "r");
		if(fp==NULL){
			printf("file could not be opened. Exitting\n");
			exit(1);
		}
		else{
			printf("opened %s.\n", filename);
		}

		//file that will be written to
		//		sprintf(filename2, "MCpolysHamL%dM%dspan%drun%dnum%d_shr.txt", L, M, span, run, curfile);
		sprintf(filename2, "ImpCreatorPrime_PKshrunkPolysL%dM%dspan%dnum%d_8_19s_shr.txt", L, M, span, curfile);
		fp2 = fopen(filename2, "w");
		if(fp2==NULL){
			printf("file could not be opened. Exitting\n");
			exit(1);
		}
		else{
			printf("created/overwrote %s.\n", filename2);
		}
		//print UofS to shrunk file
		fprintf(fp2, "UofS\n");

		fgets(line, 10, fp);

		while( fgets(line, 10, fp) != NULL ){
			sscanf(line, "%d %d %d", &x, &y, &z);
//			printf("starting point: (%d,%d,%d)\n", x, y, z);
			if(x==-999){
				break;
			}
			numpolys++;

			//fill polygon[]
			//perform shrink
			//output shrunk poly

			for(i=0; i<=(L+1)*(M+1)*(span+1)-1; i++){
				Polygon[i]=0;
			}
			index=0;
			while(1){
				fgets(line, 10, fp);
				sscanf(line, "%d", &direction);
				if(direction==-111){
					break;
				}
				else{
					Polygon[index]=direction;
					index++;
				}
			}

			//Polygon is filled. index is length
/*			printf("Polygon is: ");
			for(i=0; i<=index-1; i++){
				printf("%d ", Polygon[i]);
			}
			printf("\n");
*/
			bfacf_shrink(x,y,z,Polygon, index);

			//print shrunken Polygon
			fprintf(fp2, "0 0 0\n");
			for(i=0; i<=(L+1)*(M+1)*(span+1)-1; i++){
				if(Polygon[i]==0){
					break;
				}
				else{
					fprintf(fp2, "%d\n", Polygon[i]);
				}
			}
			fprintf(fp2, "-111\n");

/*			int newx=0;
			int newy=0;
			int newz=0;
			int maxx=0;
			int maxy=0;
			int maxz=0;
			int minx=0;
			int miny=0;
			int minz=0;
			for(i=0; i<=(L+1)*(M+1)*(span+1)-1; i++){
				if(Polygon[i]==1){
					newx++;
					if(newx>maxx){
						maxx = newx;
					}
				}
				else if(Polygon[i]==2){
					newx--;
					if(newx<minx){
						minx = newx;
					}
				}
				else if(Polygon[i]==3){
					newy++;
					if(newy>maxy){
						maxy = newy;
					}
				}
				else if(Polygon[i]==4){
					newy--;
					if(newy<miny){
						miny = newy;
					}
				}
				else if(Polygon[i]==5){
					newz++;
					if(newz>maxz){
						maxz = newz;
					}
				}
				else if(Polygon[i]==6){
					newz--;
					if(newz<minz){
						minz = newz;
					}
				}
			}
			printf("maxx=%d, maxy=%d, maxz=%d\n", maxx, maxy, maxz);
			printf("minx=%d, miny=%d, minz=%d\n", minx, miny, minz);
*/
		}
		fprintf(fp2, "-999\n");
	}

	printf("There were %d polys.\n", numpolys);

	return 0;
}





#include "bfacf_shrink_v3.c"




















