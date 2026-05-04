/*
This program converts a text file which contains polygons in UofS format into
a new polygon file in length-points (LP) format.

PARAMETERS:
argv[1] = Original UofS polygon filename
argv[2] = Output filename
argv[3] = maximum length of a polygon
*/

#include <stdio.h>
#include <stdlib.h>
//#include <string.h>

int main(int argc, char* argv[]){

	if(argc!=4){
		printf("Only %d parameters passed. Need 3. Exitting!\n", argc-1);
		exit(1);
	}

	int i;	//looping/index variable

	char *USfilename = argv[1];	//UofS file
	char *LPfilename = argv[2];	//Output file
	unsigned int maxLength = atoi(argv[3]);	//max length of a polygon

	FILE* fpUS = fopen(USfilename, "r");
	if(fpUS != NULL){
		printf("Opened '%s'\n", USfilename);
	}
	else{
		printf("Couldn't open '%s' (fpUS). Exitting!\n", USfilename);
		exit(1);
	}

	FILE* fpLP = fopen(LPfilename, "w");
	if(fpLP != NULL){
		printf("Created/overwrote: '%s'\n", LPfilename);
	}
	else{
		printf("Couldn't open '%s' (fpLP). Exitting!\n", LPfilename);
		exit(1);
	}

	int *poly;	//polygon that is read (first three entries are starting x,y,z)
	poly = (int*) calloc( maxLength+3, sizeof(int) );

	char line[100];	//line of text from UofS file as it's read
	unsigned int length;	//length of polygon being read

	int x,y,z;	//holds coordinates as polygon is printed

/////////START READING FILE

	fgets(line, 100, fpUS);	//read "UofS" header to get it out of the way

	while( fgets(line, 100, fpUS)!=NULL ){	//while there's a line to read
		//reset poly
		for(i=0; i<maxLength+3; i++){
			poly[i]=0;
			i++;
		}

		//fill poly[]
		sscanf(line, "%d %d %d", &poly[0], &poly[1], &poly[2]);	//scan in starting point
//		printf("Starting point: %d,%d,%d\n", poly[0], poly[1], poly[2]);

		if(poly[0]==-999){	//then end of file is reached
			break;
		}

		i=3;		//current index of poly[] to be filled
		while( fgets(line, 100, fpUS) != NULL ){
			sscanf(line, "%d", &poly[i]);
			if(poly[i]==-111){	//end of polygon
				break;
			}
			else{
				i++;
			}
		}
		length = i-3;

		//poly is now filled appropriately.
/*		printf("filled poly is: ");
		i=3;
		while(poly[i]!=0){
			printf("%d ", poly[i]); 
			i++;
		}
		printf("The length is %d\n",i-4);
*/

		fprintf(fpLP, "%d", length);		//print length to file

		x=poly[0];
		y=poly[1];
		z=poly[2];
		//print first point
		fprintf(fpLP, " %d %d %d", x, y, z);
		i=3;
		while(i<length+2){	//dont want to print last point (as it will be the same as the first point)
			if(poly[i]==1){
				x++;
			}
			else if(poly[i]==2){
				x--;
			}
			else if(poly[i]==3){
				y++;
			}
			else if(poly[i]==4){
				y--;
			}
			else if(poly[i]==5){
				z++;
			}
			else if(poly[i]==6){
				z--;
			}
			fprintf(fpLP, " %d %d %d", x, y, z);
			i++;
		}
		fprintf(fpLP, "\n");

	}

	fclose(fpUS);
	fclose(fpLP);

	return 0;
}




