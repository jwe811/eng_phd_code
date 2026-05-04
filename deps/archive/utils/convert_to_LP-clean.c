/*
This program takes in a name of a UofS polygon file and will convert it to a polygon
file with form length, followed by x y z coordinates.This new polygon file will
be saved under the same name, just add "_LP" to the end of the name (before the .txt).
*/

//NOTE: only works for polygons with length <(maxlen-3)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char* argv[]){

/*
	const int L = atoi(argv[1]);
	const int M = atoi(argv[2]);
	const int span = atoi(argv[3]);
	const int run = atoi(argv[4]);
	const int numfiles = atoi(argv[5]);	//total number of files
*/
	const int maxlen=15000;	//max length of a polygon

	int i;
	char *filename = argv[1];

	FILE* fpUofS = fopen(filename, "r");

	if(fpUofS != NULL){
		printf("Opened '%s'\n", filename);
	}
	else{
		printf("Couldn't open '%s', fpUofS is pointing at NULL. QUITTING!\n", filename);
		exit(1);
	}

	//create/overwrite new .txt file
	//remove .txt from filename

	//find end of file name
	i=0;
	while(filename[i]!='\0'){
		i++;
	}
	//filename has '\0' in position i of the array
//	printf("filename has %d characters in it.\n", i);

	char UofSfilenamepts[i+4];	//1 for array correction, plus 3 for _LP
	strncpy(UofSfilenamepts, filename, i);
	UofSfilenamepts[i-4]='\0';	//erase .txt
//	printf("UofSfilenamepts=%s\n", UofSfilenamepts);


	strcat(UofSfilenamepts, "_LP.txt");	//add pts.txt
//	printf("UofSfilenamepts now is %s\n", UofSfilenamepts);

	FILE* fpPts = fopen(UofSfilenamepts, "w");

	if(fpPts != NULL){
		printf("Created/overwrote: '%s'\n", UofSfilenamepts);
	}
	else{
		printf("file pointer fpPts is pointing to NULL\n");
		exit(1);
	}



//CONVERT UofS FORMAT INTO PTS FORMAT






//	ALGORITHM
//	store poly in UofS format while recording length
//	print length to file
//	use stored poly to print xyz coordinates

	int poly[maxlen];
	for(i=0; i<maxlen; i++){
		poly[i]=0;
	}
	//poly[0] is first x coord
	//poly[1] is first y coord
	//poly[2] is first z coord

	char line[10];
	int x, y, z;
	int len=0;

	fgets(line, 10, fpUofS);	//read first line to get it out of the way
//	printf("first line=%s", line);

	while( fgets(line, 10, fpUofS)!=NULL ){
		//reset poly
		poly[0]=0;
		poly[1]=0;
		poly[2]=0;
		i=3;
		while(poly[i]!=0){
			poly[i]=0;
			i++;
		}

		//fill poly and count length.
		sscanf(line, "%d %d %d", &poly[0], &poly[1], &poly[2]);
//		printf("poly0=%d,poly1=%d,poly2=%d\n", poly[0], poly[1], poly[2]);
		if(poly[0]==-999){
			break;
		}

		len=0;	//length
		i=3;		//current index of poly to be filled
		while(1){

			fgets(line, 10, fpUofS);
			sscanf(line, "%d", &poly[i]);
//			printf("poly[%d]=%d\n", i, poly[i]);
			if(poly[i]==-111){
				break;
			}
			i++;
			len++;
		}

		//poly is filled appropriately.
/*		printf("filled poly is: ");
		i=3;
		while(poly[i]!=0){
			printf("%d ", poly[i]); 
			i++;
		}
		printf("The length is %d\n",len);
*/

		//print length to file
		fprintf(fpPts, "%d", len);

		x=poly[0];
		y=poly[1];
		z=poly[2];
		//print first point
		fprintf(fpPts, " %d %d %d", x, y, z);
		i=3;
		while(i<len+2){	//dont want to print last point (as it will be the same as the first point)
			if(poly[i]==1){
				x++;
			}
			if(poly[i]==2){
				x--;
			}
			if(poly[i]==3){
				y++;
			}
			if(poly[i]==4){
				y--;
			}
			if(poly[i]==5){
				z++;
			}
			if(poly[i]==6){
				z--;
			}
			fprintf(fpPts, " %d %d %d", x, y, z);
			i++;
		}
		fprintf(fpPts, "\n");
	}

	fclose(fpUofS);
	fclose(fpPts);


	return 0;
}




