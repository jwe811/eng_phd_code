//Takes a file(s) with polygons in UofS format. Outputs how of each length of polygons there are.






#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>


main(void)
{
	printf("\nmain starting\n\n");

	int i,j,k;
	int x,y,z;
	char filename[50];
	FILE* fp;
	char line[7];
	int filenum=1;
	int polycount=0;
	int direc=0;
	int curlength=0;


	//CHANGEABLE
	const int L=4;
	const int M=1;
	const int totalspan=2;
	const int totalfilenum=1;

	const int maxlength=(L+1)*(M+1)*(totalspan+1);
	int counts[maxlength+1];
	//initialize counts[] to zero
	for(i=0; i<=maxlength; i++){
		counts[i]=0;
	}
	

	while(filenum<=totalfilenum){

		//CHANGEABLE
		sprintf(filename, "createdpolysL%dM%dspan%d.txt", L, M, totalspan);
//		sprintf(filename, "ImpCreatorAllPolysL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);

		fp = fopen(filename, "r");	//Open file to read

		if(fp != NULL){
			printf("Opened file '%s'\n", filename);
		}
		else{
			printf("file pointer is pointing to NULL\n");
			exit(1);
		}

		fgets(line, 7, fp);	//read first line (UofS) to get it out of the way
//		printf("first line=%s", line);

		while( fgets(line, 7, fp)!=NULL ){	//while there's a line to read
			sscanf(line, "%d %d %d", &x, &y, &z);	//read in starting point
			if(x==-999){	//end of file
				break;
			}
			polycount++;

			curlength=0;
			direc=0;
			while(1){	//while still in the polygon
				fgets(line, 7, fp);
				sscanf(line, "%d", &direc);	//get direction

				if(direc==-111){
					break;	//break if end of the polygon
				}
				else{
					curlength++;
				}
			}
			counts[curlength]++;
		}
		filenum++;
	}
	



	printf("%d SAPs read in\n", polycount);
	for(i=0; i<=maxlength; i++){
		if(counts[i]>0){
			printf("There are %d SAPs with length %d\n", counts[i], i);
		}
	}





	printf("\nmain finished\n\n");
}

