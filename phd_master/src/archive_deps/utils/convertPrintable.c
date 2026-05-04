void convertUofStoPrintable(char* UofSfilename){
/*
This function takes in a name of a UofS polygon file (with less than 1000 polygons) (parameter) and will convert it to polygon
files that are sequences of (x,y,z) points (Three integers per line). Each polygon will be its own file.
These new polygon files will be saved under the same name, just add "_i" to the end of the name (before
the .txt), where i is the polygon number in the UofS file, .
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
	
	char UofSfilenameNew[i+4];	//1 for array correction, plus 3 for "_iii"
/*
	strncpy(UofSfilenameNew, UofSfilename, i);
	UofSfilenameNew[i-4]='\0';	//erase .txt
//	printf("UofSfilenameNew=%s\n", UofSfilenameNew);

	strcat(UofSfilenameNew, "_1.txt");	//add _1.txt
//	printf("UofSfilenameNew now is %s\n", UofSfilenameNew);
*/
	FILE* fpNew; //will point to new file


//CONVERT UofS into multiple polygons files

	char line[7];
	int x, y, z;
	int direction=0;


	fgets(line, 7, fpUofS);	//read first line to get it out of the way
//	printf("first line=%s", line);

	/*ALGORITHM
	read first point
	rewrite first point to fpNew
	read in direction
	if direction != -111
	record appropriate point
	else newline and go to next poly and create new file.
	*/

	int curPoly=0;
	char toAdd[9];	//will hold _111.txt
	while( fgets(line, 7, fpUofS)!=NULL ){
		curPoly++;
		sscanf(line, "%d %d %d", &x, &y, &z);
		if(x==-999){
			break;
		}

		strncpy(UofSfilenameNew, UofSfilename, i);
		UofSfilenameNew[i-4]='\0';	//erase .txt
//		printf("UofSfilenameNew=%s\n", UofSfilenameNew);
		sprintf(toAdd, "_%d.txt", curPoly);
		strcat(UofSfilenameNew, toAdd);	//add _iii.txt
//		printf("UofSfilenameNew now is %s\n", UofSfilenameNew);

		FILE* fpNew = fopen(UofSfilenameNew, "w");

		if(fpNew != NULL){
			printf("Created/overwrote: '%s'\n", UofSfilenameNew);
		}
		else{
			printf("file pointer fpNew is pointing to NULL\n");
			exit(1);
		}

		fprintf(fpNew, "%d %d %d\n", x, y, z);
		direction=0;
		while(1){
			fgets(line, 7, fpUofS);
			sscanf(line, "%d", &direction);
			if(direction==1){
				x=x+1;
				fprintf(fpNew, "%d %d %d\n", x, y, z);
			}
			else if(direction==2){
				x=x-1;
				fprintf(fpNew, "%d %d %d\n", x, y, z);
			}
			else if(direction==3){
				y=y+1;
				fprintf(fpNew, "%d %d %d\n", x, y, z);
			}
			else if(direction==4){
				y=y-1;
				fprintf(fpNew, "%d %d %d\n", x, y, z);
			}
			else if(direction==5){
				z=z+1;
				fprintf(fpNew, "%d %d %d\n", x, y, z);
			}
			else if(direction==6){
				z=z-1;
				fprintf(fpNew, "%d %d %d\n", x, y, z);
			}
			else{	//direction=-111 (end of polygon)
//				fgets(line, 7, fpUofS);
				break;
			}
		}
		fclose(fpNew);
	}

	fclose(fpUofS);
}












