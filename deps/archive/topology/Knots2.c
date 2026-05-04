//This program takes a set of poly files and the corresponding knot id files, the seperates the polys by their knot type.
//i.e. Creates files with just 0_1, or just 3_1, so on.



#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>


/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED	**/
#define L 2
#define M 1
#define totalspan 4
#define numtxtfiles 11	//number of polygon files that are in the folder that need to be read
#define R (1)

int namesize = 50;	//max size of UofSfilename

/*struct point {
	int x;
	int y;
	int z;
};*/

main(void){

	char filename[namesize];	//string that will hold filename
	char knotfile[namesize];
	char outputunknot[namesize];
	char outputtrefoil[namesize];
	char outputnegativetrefoil[namesize];
	char outputfigureeight[namesize];
	char outputcinquefoil[namesize];
	char outputnegativecinquefoil[namesize];
	char outputthreetwist[namesize];
	char outputnegativethreetwist[namesize];
	char outputpretzel[namesize];
	char outputnegativepretzel[namesize];
	char outputbothtrefoils[namesize];
	FILE* fp;	//file pointer
	FILE* kf;
	FILE* uk;
	FILE* pt;
	FILE* nt;
	FILE* fe;
	FILE* pc;
	FILE* nc;
	FILE* p3;
	FILE* n3;
	FILE* pp;
	FILE* np;
	FILE* bt;
	int j, l, m;
	j = 0; //******************************************************************************************************
	int lUnKnot = 0;
	int lTrefoil = 0;
	int lNegativeTrefoil = 0;
	int lFigureEight = 0;
	int lCinquefoil = 0;
	int lNegativeCinquefoil = 0;
	int lThreeTwist = 0;
	int lNegativeThreeTwist = 0;
	int lPretzel = 0;
	int lNegativePretzel = 0;
	int lBothTrefoils = 0;
	int unKnotFile = 1;
	int trefoilFile = 1;
	int negativeTrefoilFile = 1;
	int figureEightFile = 1;
	int cinquefoilFile = 1;
	int negativeCinquefoilFile = 1;
	int threeTwistFile = 1;
	int negativeThreeTwistFile = 1;
	int pretzelFile = 1;
	int negativePretzelFile = 1;
	int bothTrefoilsFile = 1;
	sprintf(outputunknot, "CreatorAllPolysL%dM%dspan%dnum%d_0_1.txt", L, M, totalspan, unKnotFile); //Assign filenames to each output file.
	sprintf(outputtrefoil, "CreatorAllPolysL%dM%dspan%dnum%d_3_1.txt", L, M, totalspan, trefoilFile);
	sprintf(outputnegativetrefoil, "CreatorAllPolysL%dM%dspan%dnum%d_3_1s.txt", L, M, totalspan, negativeTrefoilFile);
	sprintf(outputfigureeight, "CreatorAllPolysL%dM%dspan%dnum%d_4_1.txt", L, M, totalspan, figureEightFile);
	sprintf(outputcinquefoil, "CreatorAllPolysL%dM%dspan%dnum%d_5_1.txt", L, M, totalspan, cinquefoilFile);
	sprintf(outputnegativecinquefoil, "CreatorAllPolysL%dM%dspan%dnum%d_5_1s.txt", L, M, totalspan, negativeCinquefoilFile);
	sprintf(outputthreetwist, "CreatorAllPolysL%dM%dspan%dnum%d_5_2.txt", L, M, totalspan, threeTwistFile);
	sprintf(outputnegativethreetwist, "CreatorAllPolysL%dM%dspan%dnum%d_5_2s.txt", L, M, totalspan, negativeThreeTwistFile);
	sprintf(outputpretzel, "CreatorAllPolysL%dM%dspan%dnum%d_6_1.txt", L, M, totalspan, pretzelFile);
	sprintf(outputnegativepretzel, "CreatorAllPolysL%dM%dspan%dnum%d_6_1s.txt", L, M, totalspan, negativePretzelFile);
	sprintf(outputbothtrefoils, "CreatorAllPolysL%dM%dspan%dnum%d_6c_2.txt", L, M, totalspan, figureEightFile);

	uk = fopen(outputunknot, "w");
	if(uk != NULL)
	{
		printf("Opened: '%s'\n", outputunknot);
	}
	else
	{
		printf("file pointer uk is pointing to NULL\n");
		exit(1);
	}
	pt = fopen(outputtrefoil, "w");
	if(pt != NULL)
	{
		printf("Opened: '%s'\n", outputtrefoil);
	}
	else
	{
		printf("file pointer pt is pointing to NULL\n");
		exit(1);
	}
	nt = fopen(outputnegativetrefoil, "w");
	if(nt != NULL)
	{
		printf("Opened: '%s'\n", outputnegativetrefoil);
	}
	else
	{
		printf("file pointer nt is pointing to NULL\n");
		exit(1);
	}
	fe = fopen(outputfigureeight, "w");
	if(fe != NULL)
	{
		printf("Opened: '%s'\n", outputfigureeight);
	}
	else
	{
		printf("file pointer fe is pointing to NULL\n");
		exit(1);
	}
	pc = fopen(outputcinquefoil, "w");
	if(pc != NULL)
	{
		printf("Opened: '%s'\n", outputcinquefoil);
	}
	else
	{
		printf("file pointer pc is pointing to NULL\n");
		exit(1);
	}
	nc = fopen(outputnegativecinquefoil, "w");
	if(nc != NULL)
	{
		printf("Opened: '%s'\n", outputnegativecinquefoil);
	}
	else
	{
		printf("file pointer nc is pointing to NULL\n");
		exit(1);
	}
	p3 = fopen(outputthreetwist, "w");
	if(p3 != NULL)
	{
		printf("Opened: '%s'\n", outputthreetwist);
	}
	else
	{
		printf("file pointer p3 is pointing to NULL\n");
		exit(1);
	}
	n3 = fopen(outputnegativethreetwist, "w");
	if(n3 != NULL)
	{
		printf("Opened: '%s'\n", outputnegativethreetwist);
	}
	else
	{
		printf("file pointer n3 is pointing to NULL\n");
		exit(1);
	}
	pp = fopen(outputpretzel, "w");
	if(pp != NULL)
	{
		printf("Opened: '%s'\n", outputpretzel);
	}
	else
	{
		printf("file pointer pp is pointing to NULL\n");
		exit(1);
	}
	np = fopen(outputnegativepretzel, "w");
	if(np != NULL)
	{
		printf("Opened: '%s'\n", outputnegativepretzel);
	}
	else
	{
		printf("file pointer np is pointing to NULL\n");
		exit(1);
	}
	bt = fopen(outputbothtrefoils, "w");
	if(bt != NULL)
	{
		printf("Opened: '%s'\n", outputbothtrefoils);
	}
	else
	{
		printf("file pointer bt is pointing to NULL\n");
		exit(1);
	}
	for(m=1; m<=numtxtfiles; m++){	//for each file
		sprintf(filename, "CreatorAllPolysL%dM%dspan%dnum%d.txt", L, M, totalspan, m);	//assign filename (UofS)
		//There are up to 900,001 polygons per text file.
		sprintf(knotfile, "CreatorAllPolysL%dM%dspan%dnum%d.txt_knot", L, M, totalspan, m);
		fp = fopen(filename, "r");
		kf = fopen(knotfile, "r");		

		if(fp != NULL){
			printf("Opened: '%s'\n", filename);
		}
		else{
			printf("file pointer fp is pointing to NULL\n");
			exit(1);
		}

		if(kf != NULL){
			printf("Opened: '%s'\n", knotfile);
		}
		else{
			printf("file pointer kf is pointing to NULL\n");
			exit(1);
		}
		//First line of each file will be UofS junk.
		char junk[7];
		fgets(junk, 7, fp);
		while ((!feof(fp)) && (!feof(kf)))
		{
			//char line[7];
			//fgets(line, 7, fp);
			int x, y, z;
			fscanf (fp, "%d", &x);
			fscanf (fp, "%d", &y);
			fscanf (fp, "%d", &z);
			//printf("%d\n", j);
			char c;
			int counter = 0;
			char str[5] = "_____";
			do
			{
				c = getc(kf);
				str[counter] = c;
				counter = counter + 1;
			} while ((c != '\n') && (!feof(kf))); //Read from the knot file.
			//printf("%s\n", str);
			if ((str[0] == '0') && (str[1] == '_') && (str[2] == '1')) //unknot.
			{
				if (j != -999)
				{
					if (lUnKnot > 900000) //Start a new file, increment the file count, and reset the polygon count.
					{
						fprintf(uk, "-999");
						fclose(uk);
						lUnKnot = 0;
						unKnotFile = unKnotFile + 1;
						sprintf(outputunknot, "CreatorAllPolysL%dM%dspan%dnum%d_0_1.txt", L, M, totalspan, unKnotFile); //~~~~~
						uk = fopen(outputunknot, "w");
						if(uk != NULL)
						{
							printf("Opened: '%s'\n", outputunknot);
						}
						else
						{
							printf("file pointer uk is pointing to NULL\n");
							exit(1);
						}
					}
					if (lUnKnot == 0)
					{
						fprintf(uk, "UofS\n"); //header for the start of a new file.
					}
					lUnKnot = lUnKnot + 1;
					fprintf(uk, "%d %d %d\n", x, y, z); //~~~~~~~~~~~~~~~~
					fscanf(fp, "%d", &j); //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					while (j != -111)
					{
						if (j != 0)
						{
							fprintf(uk, "%d\n", j);
						}						
						fscanf (fp, "%d", &j);
						//printf("%d\n", i);
						//printf("%d\n", j);
					}
					fprintf(uk, "-111\n");
				}
			} //One polygon has been read in.
			else if ((str[0] == '3') && (str[1] == '_') && (str[2] == '1')) //trefoil
			{
				if (str[3] == 's') //negative orientation
				{
					if (j != -999)
					{
						if (lNegativeTrefoil > 900000) //Start a new file, increment the file count, and reset the polygon count.
						{
							fprintf(nt, "-999");
							fclose(nt);
							lNegativeTrefoil = 0;
							negativeTrefoilFile = negativeTrefoilFile + 1;
							sprintf(outputnegativetrefoil, "CreatorAllPolysL%dM%dspan%dnum%d_3_1s.txt", L, M, totalspan, negativeTrefoilFile); //~~~~~~~~~~~~~~~~~~~~~~
							nt = fopen(outputnegativetrefoil, "w");
							if(nt != NULL)
							{
								printf("Opened: '%s'\n", outputnegativetrefoil);
							}
							else
							{
								printf("file pointer nt is pointing to NULL\n");
								exit(1);
							}
						}
						if (lNegativeTrefoil == 0)
						{
							fprintf(nt, "UofS\n"); //header for the start of a new file.
						}
						lNegativeTrefoil = lNegativeTrefoil + 1;
						fprintf(nt, "%d %d %d\n", x, y, z); //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
						fscanf(fp, "%d", &j); //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
						while (j != -111)
						{
							if (j != 0)
							{
								fprintf(nt, "%d\n", j);
							}
							fscanf (fp, "%d", &j);
						}
						fprintf(nt, "-111\n");
					}
				}
				else //positive orientation
				{
					if (j != -999)
					{
						if (lTrefoil > 900000) //Start a new file, increment the file count, and reset the polygon count.
						{
							fprintf(pt, "-999");
							fclose(pt);
							lTrefoil = 0;
							trefoilFile = trefoilFile + 1; 
							sprintf(outputtrefoil, "CreatorAllPolysL%dM%dspan%dnum%d_3_1.txt", L, M, totalspan, trefoilFile);
							pt = fopen(outputtrefoil, "w");
							if(pt != NULL)
							{
								printf("Opened: '%s'\n", outputtrefoil);
							}
							else
							{
								printf("file pointer pt is pointing to NULL\n");
								exit(1);
							}
						}
						if (lTrefoil == 0)
						{
							fprintf(pt, "UofS\n"); //header for the start of a new file.
						}
						lTrefoil = lTrefoil + 1;
						fprintf(pt, "%d %d %d\n", x, y, z);
						fscanf(fp, "%d", &j);
						while (j != -111)
						{
							if (j != 0)
							{
								fprintf(pt, "%d\n", j);
							}
							fscanf (fp, "%d", &j);
						}
						fprintf(pt, "-111\n");
					}
				}
			}
			else if ((str[0] == '4') && (str[1] == '_') && (str[2] == '1')) //figure eight
			{
				if (j != -999)
				{
					if (lFigureEight > 900000) //Start a new file, increment the file count, and reset the polygon count.
					{
						fprintf(fe, "-999");
						fclose(fe);
						lFigureEight = 0;
						figureEightFile = figureEightFile + 1;
						sprintf(outputfigureeight, "CreatorAllPolysL%dM%dspan%dnum%d_4_1.txt", L, M, totalspan, figureEightFile);
						fe = fopen(outputfigureeight, "w");
						if(fe != NULL)
						{
							printf("Opened: '%s'\n", outputfigureeight);
						}
						else
						{
							printf("file pointer fe is pointing to NULL\n");
							exit(1);
						}
					}
					if (lFigureEight == 0)
					{
						fprintf(fe, "UofS\n"); //header for the start of a new file.
					}
					lFigureEight = lFigureEight + 1;
					fprintf(fe, "%d %d %d\n", x, y, z);
					fscanf(fp, "%d", &j);
					while (j != -111)
					{
						if (j != 0)
						{
							fprintf(fe, "%d\n", j);
						}
						fscanf (fp, "%d", &j);
					}
					fprintf(fe, "-111\n");
				}
			}
			else if ((str[0] == '5') && (str[1] == '_') && (str[2] == '1')) //cinquefoil
			{
				if (str[3] == 's') //negative orientation
				{
					if (j != -999)
					{
						if (lNegativeCinquefoil > 900000) //Start a new file, increment the file count, and reset the polygon count.
						{
							fprintf(nc, "-999");
							fclose(nc);
							lNegativeCinquefoil = 0;
							negativeCinquefoilFile = negativeCinquefoilFile + 1;
							sprintf(outputnegativecinquefoil, "CreatorAllPolysL%dM%dspan%dnum%d_5_1s.txt", L, M, totalspan, negativeCinquefoilFile);
							nc = fopen(outputnegativecinquefoil, "w");
							if(nc != NULL)
							{
								printf("Opened: '%s'\n", outputnegativecinquefoil);
							}
							else
							{
								printf("file pointer nc is pointing to NULL\n");
								exit(1);
							}
						}
						if (lNegativeCinquefoil == 0)
						{
							fprintf(nc, "UofS\n"); //header for the start of a new file.
						}
						lNegativeCinquefoil = lNegativeCinquefoil + 1;
						fprintf(nc, "%d %d %d\n", x, y, z);
						fscanf(fp, "%d", &j);
						while (j != -111)
						{
							if (j != 0)
							{
								fprintf(nc, "%d\n", j);
							}
							fscanf (fp, "%d", &j);
						}
						fprintf(nc, "-111\n");
					}
				}
				else //positive orientation
				{
					if (j != -999)
					{
						if (lCinquefoil > 900000) //Start a new file, increment the file count, and reset the polygon count.
						{
							fprintf(pc, "-999");
							fclose(pc);
							lCinquefoil = 0;
							cinquefoilFile = cinquefoilFile + 1; 
							sprintf(outputcinquefoil, "CreatorAllPolysL%dM%dspan%dnum%d_5_1.txt", L, M, totalspan, cinquefoilFile);
							pc = fopen(outputcinquefoil, "w");
							if(pc != NULL)
							{
								printf("Opened: '%s'\n", outputcinquefoil);
							}
							else
							{
								printf("file pointer pc is pointing to NULL\n");
								exit(1);
							}
						}
						if (lCinquefoil == 0)
						{
							fprintf(pc, "UofS\n"); //header for the start of a new file.
						}
						lCinquefoil = lCinquefoil + 1;
						fprintf(pc, "%d %d %d\n", x, y, z);
						fscanf(fp, "%d", &j);
						while (j != -111)
						{
							if (j != 0)
							{
								fprintf(pc, "%d\n", j);
							}
							fscanf (fp, "%d", &j);
						}
						fprintf(pc, "-111\n");
					}
				}
			}
			else if ((str[0] == '5') && (str[1] == '_') && (str[2] == '2')) //three-twist
			{
				if (str[3] == 's') //negative orientation
				{
					if (j != -999)
					{
						if (lNegativeThreeTwist > 900000) //Start a new file, increment the file count, and reset the polygon count.
						{
							fprintf(n3, "-999");
							fclose(n3);
							lNegativeThreeTwist = 0;
							negativeThreeTwistFile = negativeThreeTwistFile + 1;
							sprintf(outputnegativethreetwist, "CreatorAllPolysL%dM%dspan%dnum%d_5_2s.txt", L, M, totalspan, negativeThreeTwistFile);
							n3 = fopen(outputnegativethreetwist, "w");
							if(n3 != NULL)
							{
								printf("Opened: '%s'\n", outputnegativethreetwist);
							}
							else
							{
								printf("file pointer n3 is pointing to NULL\n");
								exit(1);
							}
						}
						if (lNegativeThreeTwist == 0)
						{
							fprintf(n3, "UofS\n"); //header for the start of a new file.
						}
						lNegativeThreeTwist = lNegativeThreeTwist + 1;
						fprintf(n3, "%d %d %d\n", x, y, z);
						fscanf(fp, "%d", &j);
						while (j != -111)
						{
							if (j != 0)
							{
								fprintf(n3, "%d\n", j);
							}
							fscanf (fp, "%d", &j);
						}
						fprintf(n3, "-111\n");
					}
				}
				else //positive orientation
				{
					if (j != -999)
					{
						if (lThreeTwist > 900000) //Start a new file, increment the file count, and reset the polygon count.
						{
							fprintf(p3, "-999");
							fclose(p3);
							lThreeTwist = 0;
							threeTwistFile = threeTwistFile + 1; 
							sprintf(outputthreetwist, "CreatorAllPolysL%dM%dspan%dnum%d_5_2.txt", L, M, totalspan, threeTwistFile);
							p3 = fopen(outputthreetwist, "w");
							if(p3 != NULL)
							{
								printf("Opened: '%s'\n", outputthreetwist);
							}
							else
							{
								printf("file pointer p3 is pointing to NULL\n");
								exit(1);
							}
						}
						if (lThreeTwist == 0)
						{
							fprintf(p3, "UofS\n"); //header for the start of a new file.
						}
						lThreeTwist = lThreeTwist + 1;
						fprintf(p3, "%d %d %d\n", x, y, z);
						fscanf(fp, "%d", &j);
						while (j != -111)
						{
							if (j != 0)
							{
								fprintf(p3, "%d\n", j);
							}
							fscanf (fp, "%d", &j);
						}
						fprintf(p3, "-111\n");
					}
				}
			}
			else if ((str[0] == '6') && (str[1] == '_') && (str[2] == '1')) //pretzel
			{
				if (str[3] == 's') //negative orientation
				{
					if (j != -999)
					{
						if (lNegativePretzel > 900000) //Start a new file, increment the file count, and reset the polygon count.
						{
							fprintf(np, "-999");
							fclose(np);
							lNegativePretzel = 0;
							negativePretzelFile = negativePretzelFile + 1;
							sprintf(outputnegativepretzel, "CreatorAllPolysL%dM%dspan%dnum%d_6_1s.txt", L, M, totalspan, negativePretzelFile);
							np = fopen(outputnegativepretzel, "w");
							if(np != NULL)
							{
								printf("Opened: '%s'\n", outputnegativepretzel);
							}
							else
							{
								printf("file pointer np is pointing to NULL\n");
								exit(1);
							}
						}
						if (lNegativePretzel == 0)
						{
							fprintf(np, "UofS\n"); //header for the start of a new file.
						}
						lNegativePretzel = lNegativePretzel + 1;
						fprintf(np, "%d %d %d\n", x, y, z);
						fscanf(fp, "%d", &j);
						while (j != -111)
						{
							if (j != 0)
							{
								fprintf(np, "%d\n", j);
							}
							fscanf (fp, "%d", &j);
						}
						fprintf(np, "-111\n");
					}
				}
				else //positive orientation
				{
					if (j != -999)
					{
						if (lPretzel > 900000) //Start a new file, increment the file count, and reset the polygon count.
						{
							fprintf(pp, "-999");
							fclose(pp);
							lPretzel = 0;
							pretzelFile = pretzelFile + 1; 
							sprintf(outputpretzel, "CreatorAllPolysL%dM%dspan%dnum%d_6_1.txt", L, M, totalspan, pretzelFile);
							pp = fopen(outputpretzel, "w");
							if(pp != NULL)
							{
								printf("Opened: '%s'\n", outputpretzel);
							}
							else
							{
								printf("file pointer pp is pointing to NULL\n");
								exit(1);
							}
						}
						if (lPretzel == 0)
						{
							fprintf(pp, "UofS\n"); //header for the start of a new file.
						}
						lPretzel = lPretzel + 1;
						fprintf(pp, "%d %d %d\n", x, y, z);
						fscanf(fp, "%d", &j);
						while (j != -111)
						{
							if (j != 0)
							{
								fprintf(pp, "%d\n", j);
							}
							fscanf (fp, "%d", &j);
						}
						fprintf(pp, "-111\n");
					}
				}
			}
			else if ((str[0] == '6') && (str[1] == 'c') && (str[2] == '_') && (str[3] == '2')) //both trefoils
			{
				if (j != -999)
				{
					if (lBothTrefoils > 900000) //Start a new file, increment the file count, and reset the polygon count.
					{
						fprintf(bt, "-999");
						fclose(bt);
						lBothTrefoils = 0;
						bothTrefoilsFile = bothTrefoilsFile + 1;
						sprintf(outputbothtrefoils, "CreatorAllPolysL%dM%dspan%dnum%d_6c_2.txt", L, M, totalspan, bothTrefoilsFile);
						bt = fopen(outputbothtrefoils, "w");
						if(bt != NULL)
						{
							printf("Opened: '%s'\n", outputbothtrefoils);
						}
						else
						{
							printf("file pointer bt is pointing to NULL\n");
							exit(1);
						}
					}
					if (lBothTrefoils == 0)
					{
						fprintf(bt, "UofS\n"); //header for the start of a new file.
					}
					lBothTrefoils = lBothTrefoils + 1;
					fprintf(bt, "%d %d %d\n", x, y, z);
					fscanf(fp, "%d", &j);
					while (j != -111)
					{
						if (j != 0)
						{
							fprintf(bt, "%d\n", j);
						}
						fscanf (fp, "%d", &j);
					}
					fprintf(bt, "-111\n");
				}
			}
		}
		fclose(fp);
	}
	//Close each of the final files, after adding -999 to the end if non-empty.
	if (lUnKnot > 0)
	{
		fprintf(uk, "-999");
	}
	fclose(uk);
	if (lTrefoil > 0)
	{
		fprintf(pt, "-999");
	}
	fclose(pt);
	if (lNegativeTrefoil > 0)
	{
		fprintf(nt, "-999");
	}
	fclose(nt);
	if (lFigureEight > 0)
	{
		fprintf(fe, "-999");
	}
	fclose(fe);
	if (lCinquefoil > 0)
	{
		fprintf(pc, "-999");
	}
	fclose(pc);
	if (lNegativeCinquefoil > 0)
	{
		fprintf(nc, "-999");
	}
	fclose(nc);
	if (lThreeTwist > 0)
	{
		fprintf(p3, "-999");
	}
	fclose(p3);
	if (lNegativeThreeTwist > 0)
	{
		fprintf(n3, "-999");
	}
	fclose(n3);
	if (lPretzel > 0)
	{
		fprintf(pp, "-999");
	}
	fclose(pp);
	if (lNegativePretzel > 0)
	{
		fprintf(np, "-999");
	}
	fclose(np);
	if (lBothTrefoils > 0)
	{
		fprintf(bt, "-999");
	}
	fclose(bt);
}

//#include "../utils/convertUofStoPoints.c"
