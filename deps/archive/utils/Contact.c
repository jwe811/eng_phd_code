/*
	This program will use all text files of the form "CreatorPrimePolysL%dM%dspan%dnum%d.txt". Will convert each file into a pts file, then
	using all of the polygons from these txt files, this program will create a contact map
	for such polygons

	written by Jeremy Eng
*/


#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>


/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED	**/
#define L 2
#define M 1
#define totalspan 5
#define numtxtfiles 2	//number of polygon files that are in the folder that need to be read
#define R 1

int namesize = 50;	//max size of UofSfilename

struct point {
	int x;
	int y;
	int z;
};

main(void){

	char filename[namesize];		//string that will hold filename
	FILE* fp;	//file pointer
	int i, j, k, l, m;
	double sum[(L + 1) * (M + 1) * (totalspan + 1)/2 + 1];
	l = 0;
	for (i = 0; i <= (L + 1) * (M + 1) * (totalspan + 1)/2; i = i + 1)
	{
		sum[i] = 0;
	}
	for(m=1; m<=numtxtfiles; m++){	//for each file
		sprintf(filename, "CreatorPrimePolysL%dM%dspan%dnum%d.txt", L, M, totalspan, m);	//assign filename (UofS)

		convertUofStoPoints(filename);	//convert UofS file to pts file

		//pts format was created. that's the file we want to look at.
		sprintf(filename, "CreatorPrimePolysL%dM%dspan%dnum%dpts.txt", L, M, totalspan, m);	//assign filename (pts)

		fp = fopen(filename, "r");		

		if(fp != NULL){
			printf("Opened: '%s'\n", filename);
		}
		else{
			printf("file pointer fp is pointing to NULL\n");
			exit(1);
		}

		//ALGORITHM FOR CONTACT MAP:
		while (!feof(fp))
		{
			int X[(L + 1) * (M + 1) * (totalspan + 1)/2 + 1];
			for (i = 0; i <= (L + 1) * (M + 1) * (totalspan + 1)/2; i = i + 1)
			{
				X[i] = 0;
			}
//			printf("X[] initialized to all zeros: X[0]=%d\n", X[0]);
			struct point dummy[(L+1)*(M+1)*(totalspan + 1)];
			fscanf (fp, "%d", &j);
//			printf("scanned first thing: %d\n", j);
			k = 0;
			if (j != -999)
			{
 				while (j != -111)
 				{
//					printf("starting loop\n");
      				dummy[k].x = j;
//					printf("x-coord: %d\n", j);
//					printf("%d 1 %d\n", l, X[0]);
    					fscanf (fp, "%d", &j);
//					printf("scanned y-coord: %d\n", j);
//					printf("%d 2 %d\n", l, X[0]);
					dummy[k].y = j;
//					printf("%d 3 %d\n", l, X[0]);
					fscanf (fp, "%d", &j);
//					printf("scanned z-coord: %d\n", j);
//					printf("%d 4 %d\n", l, X[0]);
					dummy[k].z = j;
//					printf("%d 5 %d\n", l, X[0]);
					fscanf (fp, "%d", &j);
//					printf("scanned next thinger, it's j=%d\n", j);
//					printf("%d 6 %d\n", l, X[0]);
					k = k + 1;
//					printf("%d 7 %d\n", l, X[0]);
//					printf("ending loop\n");
    				}
//				printf("LOOP IS FINISHED\n");
				//printf("%d 2 %d\n", l, X[0]);
				struct point omega[k - 1];
				for (i = 0; i < k - 1; i = i+1)
				{
					omega[i] = dummy[i];
				}
				/*for (i = 0; i < k - 1; i = i+1)
				{
					printf("%d", omega[i].x);
					printf("%d", omega[i].y);
					printf("%d", omega[i].z);
					printf("\n");
				}*/
				//printf("%d %d\n", l, X[0]);
				for (i = 0; i <= (L + 1) * (M + 1) * (totalspan + 1)/2; i = i + 1)
				{
					X[i] = 0;
				}
				for (i = 0; i < k - 1; i = i + 1) //go through each pair of points.
				{
					for (j = 0; j < k - 1; j = j + 1)
					{
						int distance;
						if (i != j)
						{
							double x = omega[i].x - omega[j].x;
							double y = omega[i].y - omega[j].y;
							double z = omega[i].z - omega[j].z;
							double d = sqrt(x*x + y*y + z*z);
							if (d <= R) //A contact occurs.
							{
								distance = fmin(abs(i - j), (k - 1) - abs(i - j));
								X[distance] = X[distance] + 1;
								//printf("Contact at %d, %d. X[%d] is now %d\n", i, j, distance, X[distance]);
							}
						}
					}
				}
				//printf("%d %d\n", l, X[0]);
				X[1] = 0;
				int Y = 0;
				for (j = 0; j <= (L + 1) * (M + 1) * (totalspan + 1)/2; j = j + 1)
				{
					X[j] = X[j] / 2;
					Y = Y + X[j];
				}
				//printf("Y = %d\n", Y);
				for (j = 0; j <= (L + 1) * (M + 1) * (totalspan + 1)/2; j = j + 1)
				{
					double ratio = (double)(X[j])/(double)(Y);
					sum[j] = sum[j] + ratio;
					if (j == 0)
					{
						//printf("%d %d\n", l, X[j]);
					}
				}
				l = l + 1;
			}
		}
		fclose(fp);
	}
	double average[(L + 1) * (M + 1) * (totalspan + 1)/2 + 1];
	for (j = 0; j <= (L + 1) * (M + 1) * (totalspan + 1)/2; j = j + 1)
	{
		printf("sum = %f, l = %d\n", sum[j], l);
		average[j] = sum[j] / l;
		printf("average[%d] = %f\n", j, average[j]);
	}
}



#include "convertUofStoPoints.c"

/*
average[0] = 0.000000
average[1] = 0.000000
average[2] = 0.000000
average[3] = 0.249590
average[4] = 0.000000
average[5] = 0.170468
average[6] = 0.000000
average[7] = 0.138103
average[8] = 0.000000
average[9] = 0.113701
average[10] = 0.000000
average[11] = 0.094475
average[12] = 0.000000
average[13] = 0.082269
average[14] = 0.000000
average[15] = 0.076793
average[16] = 0.000000
average[17] = 0.074602
average[18] = 0.000000
*/

/*
average[0] = 0.000000
average[1] = 0.000000
average[2] = 0.000000
average[3] = 0.237934
average[4] = 0.000000
average[5] = 0.159057
average[6] = 0.000000
average[7] = 0.126522
average[8] = 0.000000
average[9] = 0.102900
average[10] = 0.000000
average[11] = 0.084271
average[12] = 0.000000
average[13] = 0.072349
average[14] = 0.000000
average[15] = 0.065604
average[16] = 0.000000
average[17] = 0.061809
average[18] = 0.000000
average[19] = 0.059863
average[20] = 0.000000
average[21] = 0.029692

*/
