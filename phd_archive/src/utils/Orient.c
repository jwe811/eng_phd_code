#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>


/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED	**/
#define L 2
#define M 1
#define totalspan 6
#define numtxtfiles 1	//number of polygon files that are in the folder that need to be read
#define R sqrt(2)

int namesize = 50;	//max size of UofSfilename

struct point {
	int x;
	int y;
	int z;
};

int compare(struct point pt1, struct point pt2) //Determines whether two point records are equal.
{
	if ((pt1.x == pt2.x) && (pt1.y == pt2.y) && (pt1.z == pt2.z))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int mod(int dividend, int divisor) //modulus division, but guarantees a non-negative remainder.
{
	int remainder = dividend % divisor;
	if (remainder < 0)
	{
		remainder = remainder + divisor;
	}
	return remainder;
}

int opposite(int direction) //reverses the direction.
{
	if (direction % 2 == 0)
	{
		return direction - 1;
	}
	else
	{
		return direction + 1;
	}
}

main(void){
	char filename[namesize];	//string that will hold filename
	char output[namesize];
	FILE* fp;	//file pointer
	FILE* out;
	int j, k, index, m;
	for(m=1; m<=numtxtfiles; m++){	//for each file
		//sprintf(output, "CreatorPrimePolysL%dM%dspan%dnum%d_oriented.txt", L, M, totalspan, m); //******************************************************************************
		sprintf(output, "CreatorPrimePolysL2M1span6num1_3_1_mirrorY_oriented.txt"); //****************************************************************************************************
		out = fopen(output, "w");
		if(out != NULL)
		{
			printf("Opened: '%s'\n", output);
		}
		else
		{
			printf("file pointer out is pointing to NULL\n");
			exit(1);
		}		
		//sprintf(filename, "CreatorPrimePolysL%dM%dspan%dnum%d.txt", L, M, totalspan, m); //************************************************************************************
		sprintf(filename, "CreatorPrimePolysL2M1span6num1_3_1_mirrorY.txt"); //****************************************************************************************************
		//There are up to 900,001 polygons per text file.
		fp = fopen(filename, "r");
		if(fp != NULL){
			printf("Opened: '%s'\n", filename);
		}
		else{
			printf("file pointer fp is pointing to NULL\n");
			exit(1);
		}
		//First line of each file will be UofS junk.
		char junk[7];
		fgets(junk, 7, fp);
		fprintf(out, "UofS\n");
		while (!feof(fp))
		{
			//char line[7];
			//fgets(line, 7, fp);
			struct point start;
			int polygon[(L+1)*(M+1)*(totalspan+1)];
			int x, y, z;
			fscanf(fp, "%d", &x);
			if (x != -999)
			{
				fscanf(fp, "%d", &y);
				fscanf(fp, "%d", &z);
				start.x = x;
				start.y = y;
				start.z = z;
				//fprintf(out, "%d %d %d\n", L - x, y, 0);
				for (k = 0; k < (L+1) * (M+1) * (totalspan + 1); k = k+1) //read in the polygon.
				{
					fscanf(fp, "%d", &j);
					polygon[k] = j;
				}
				index = 0;
				struct point origin;
				origin.x = 0;
				origin.y = 0;
				origin.z = 0;
				while (compare(start, origin) == 0) //locate (0, 0, 0).
				{
					if (polygon[index] == 1)
					{
						start.x = start.x + 1;
					}
					else if (polygon[index] == 2)
					{
						start.x = start.x - 1;
					}
					else if (polygon[index] == 3)
					{
						start.y = start.y + 1;
					}
					else if (polygon[index] == 4)
					{
						start.y = start.y - 1;
					}
					else if (polygon[index] == 5)
					{
						start.z = start.z + 1;
					}
					else if (polygon[index] == 6)
					{
						start.z = start.z - 1;
					}
					index = index + 1;
				}
				//Now, polygon[index] will be the first direction and polygon[mod(index - 1, (L+1)*(M+1)*(totalspan+1))] will be the last.
				//Four cases: First direction is 5 (+z), last direction is 6 (-z), first direction is 3 (+y), or last direction is 4 (-y).
				if (polygon[index] == 5) //Case 1
				{
					fprintf(out, "%d %d %d\n", 0, 0, 0);
					for (k = 0; k < (L+1) * (M+1) * (totalspan + 1); k = k+1)
					{
						j = polygon[mod(index + k, (L+1)*(M+1)*(totalspan+1))]; //Go forward around the polygon.
						fprintf(out, "%d\n", j);
					}
				}
				else if (polygon[mod(index - 1, (L+1)*(M+1)*(totalspan+1))] == 6) //Case 2
				{
					fprintf(out, "%d %d %d\n", 0, 0, 0);
					for (k = 0; k < (L+1) * (M+1) * (totalspan + 1); k = k+1)
					{
						j = polygon[mod(index - k - 1, (L+1)*(M+1)*(totalspan+1))]; //Go backwards around the polygon.
						j = opposite(j); //Reserse each direction.
						fprintf(out, "%d\n", j);
					}
				}
				else if (polygon[index] == 3) //Case 3
				{
					fprintf(out, "%d %d %d\n", 0, 0, 0);
					for (k = 0; k < (L+1) * (M+1) * (totalspan + 1); k = k+1)
					{
						j = polygon[mod(index + k, (L+1)*(M+1)*(totalspan+1))]; //Go forward around the polygon.
						fprintf(out, "%d\n", j);
					}
				}
				else if (polygon[mod(index - 1, (L+1)*(M+1)*(totalspan+1))] == 4) //Case 4
				{
					fprintf(out, "%d %d %d\n", 0, 0, 0);
					for (k = 0; k < (L+1) * (M+1) * (totalspan + 1); k = k+1)
					{
						j = polygon[mod(index - k - 1, (L+1)*(M+1)*(totalspan+1))]; //Go backwards around the polygon.
						j = opposite(j); //Reserse each direction.
						fprintf(out, "%d\n", j);
					}
				}
			fscanf(fp, "%d", &j); //Deal with the -111 now.
			fprintf(out, "%d\n", j);
			}
		}
		fprintf(out, "%d", -999);
		fclose(fp);
		fclose(out);
	}
}

//#include "convertUofStoPoints.c"
