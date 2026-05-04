void convertUofStoPoints(char* UofSfilename){
/*
This function takes in a name of a UofS polygon file (first parameter) and will convert it to a polygon
file that is just a sequence of (x,y,z) points (Three integers per line). This new polygon file will
be saved under the same name, just add "pts" to the end of the name (before the .txt).
*/
	int i;
	FILE* fpUofS = fopen(UofSfilename, "r");

	if(fpUofS != NULL){
		printf("Opened '%s'\n", UofSfilename);
	}
	else{
		printf("Couldn't open '%s', fpUofS is pointing at NULL. QUITTING!\n", UofSfilename);
		exit(1);
	}

	//create/overwrite new .txt file
	//remove .txt from UofSfilename

	//find end of file name
	i=0;
	while(UofSfilename[i]!='\0'){
		i++;
	}
	//UofSfilename has '\0' in position i of the array
//	printf("filename has %d characters in it.\n", i);

	char UofSfilenamepts[i+4];	//1 for array correction, plus 3 for pts
	strncpy(UofSfilenamepts, UofSfilename, i);
	UofSfilenamepts[i-4]='\0';	//erase .txt
//	printf("UofSfilenamepts=%s\n", UofSfilenamepts);


	strcat(UofSfilenamepts, "pts.txt");	//add pts.txt
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

	char line[7];
	int x, y, z;
	int direction=0;


	fgets(line, 7, fpUofS);	//read first line to get it out of the way
//	printf("first line=%s", line);

	/*ALGORITHM
	read first point
	rewrite first point to fpPts
	read in direction
	if direction != -111
	record appropriate point
	else record -111 and go to next poly.
	*/
	

	while( fgets(line, 7, fpUofS)!=NULL ){
		sscanf(line, "%d %d %d", &x, &y, &z);
		if(x==-999){
			fprintf(fpPts, "-999");
			break;
		}
		fprintf(fpPts, "%d %d %d\n", x, y, z);
		direction=0;
		while(1){
			fgets(line, 7, fpUofS);
			sscanf(line, "%d", &direction);
			if(direction==1){
				x=x+1;
				fprintf(fpPts, "%d %d %d\n", x, y, z);
			}
			else if(direction==2){
				x=x-1;
				fprintf(fpPts, "%d %d %d\n", x, y, z);
			}
			else if(direction==3){
				y=y+1;
				fprintf(fpPts, "%d %d %d\n", x, y, z);
			}
			else if(direction==4){
				y=y-1;
				fprintf(fpPts, "%d %d %d\n", x, y, z);
			}
			else if(direction==5){
				z=z+1;
				fprintf(fpPts, "%d %d %d\n", x, y, z);
			}
			else if(direction==6){
				z=z-1;
				fprintf(fpPts, "%d %d %d\n", x, y, z);
			}
			else{	//direction=-111 (end of polygon)
				fprintf(fpPts, "%d\n", direction);
				break;
			}
		}
	}

	fclose(fpUofS);
	fclose(fpPts);
}












