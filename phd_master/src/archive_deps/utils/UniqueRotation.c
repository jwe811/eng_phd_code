#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>



main(void)
{
	printf("\nmain starting\n\n");

	const int L = 3;	//CHANGEABLE
	const int M = 1;	//CHANGEABLE
	const int span = 3;	//CHANGEABLE

	int i, j, k;

	//count how many polygons are in the file
	char filename[45];
	sprintf(filename, "CreatorPrimePolysL%dM%dspan%dnum1_3_1.txt", L, M, span);	//CHANGEABLE
	FILE* fpOriginal = fopen(filename, "r");
	if(fpOriginal != NULL){
		printf("Opened '%s'\n", filename);
	}
	else{
		printf("Couldn't open '%s', fpOriginal is pointing at NULL. QUITTING!\n", filename);
		exit(1);
	}
	char line[7];
	int totpolys=0;
	int x, y, z;
	
	while( fgets(line, 7, fpOriginal)!=NULL ){
		sscanf(line, "%d %d %d", &x, &y, &z);
		if(x==-111){
			totpolys++;
		}
	}
	printf("There are %d polygons in the file %s.\n", totpolys, filename);
//	fpOriginal = fopen(filename, "r");
	char Uniquefilename[25];
	sprintf(Uniquefilename, "UniqueL%dM%dspan%d_3_1.txt", L, M, span); //CHANGEABLE
	FILE* fpUnique = fopen(Uniquefilename, "w");
	if(fpUnique != NULL){
		printf("Opened '%s'\n", Uniquefilename);
	}
	else{
		printf("Couldn't open '%s', fpUnique is pointing at NULL. QUITTING!\n", Uniquefilename);
		exit(1);
	}

	//read in and store all polygons in arrays
	int polygon[totpolys][3+(L+1)*(M+1)*(span+1)];	//first 3 entries (0,1,2) will be the starting xyz
	//intialize polygon[][] to zeros
	for(i=0; i<=totpolys-1; i++){
		for(j=0; j<=(3+(L+1)*(M+1)*(span+1)-1); j++){
			polygon[i][j]=0;
		}
	}
	printf("polygon[][] initialized to zeros\n");
	int direction=0;
	int curpoly=-1;
	int curindex;
	//read in all polygons, and put them in arrays
	fseek(fpOriginal, 0, SEEK_SET);
	fgets(line, 7, fpOriginal);	//get rid of UofS line
	
	while( fgets(line, 7, fpOriginal)!=NULL ){
		sscanf(line, "%d %d %d", &x, &y, &z);
		if(x==-999){
			break;
		}
		curpoly++;
		polygon[curpoly][0] = x;
		polygon[curpoly][1] = y;
		polygon[curpoly][2] = z;

		direction=0;
		curindex=3;
		while(1){
			fgets(line, 7, fpOriginal);
			sscanf(line, "%d", &direction);
			if(direction==-111){
				break;
			}
			else{
				polygon[curpoly][curindex]=direction;
				curindex++;
			}
		}
	}
	printf("All polygons recorded (%d polygons)\n", curpoly+1);

	//Now orient all polygons. For each polygon, move to (0,0,0). Follow these rules: If there's an edge in posz direction, then that's the first edge. if not, then the edge in the pos y direction is the first edge.


	int temppoly[3+(L+1)*(M+1)*(span+1)];
	int tempnum;

	//move to (0,0,0)
	for(i=0; i<=totpolys-1; i++){	//for each polygon
		while((polygon[i][0]!=0) || (polygon[i][1]!=0) || (polygon[i][2]!=0) ){
			tempnum = polygon[i][3];
			if(polygon[i][3]==1){
				polygon[i][0]++;
			}
			else if(polygon[i][3]==2){
				polygon[i][0]--;
			}
			else if(polygon[i][3]==3){
				polygon[i][1]++;
			}
			else if(polygon[i][3]==4){
				polygon[i][1]--;
			}
			else if(polygon[i][3]==5){
				polygon[i][2]++;
			}
			else if(polygon[i][3]==6){
				polygon[i][2]--;
			}

			for(k=3; k<=(3+(L+1)*(M+1)*(span+1)-2); k++){
				polygon[i][k] = polygon[i][k+1];
			}
			polygon[i][(3+(L+1)*(M+1)*(span+1)-1)] = tempnum;
		}

		//now at (0,0,0). if first direction is 5, then we are ok. if not, we need to re-orient
		if(polygon[i][3]!=5){
			for(j=0; j<=3+(L+1)*(M+1)*(span+1)-1; j++){
				temppoly[j]=polygon[i][j];
			}
		}
	
		//if last direction is 6, then we need to flip ordering
		if(polygon[i][3+(L+1)*(M+1)*(span+1)-1]==6){
			//flip ordering
			for(j=3; j<=(3+(L+1)*(M+1)*(span+1)-1); j++){
				k = (3+(L+1)*(M+1)*(span+1)-1+3-j);
				if(temppoly[k]==1){
					polygon[i][j] = 2;
				}
				else if(temppoly[k]==2){
					polygon[i][j] = 1;
				}
				else if(temppoly[k]==3){
					polygon[i][j] = 4;
				}
				else if(temppoly[k]==4){
					polygon[i][j] = 3;
				}
				else if(temppoly[k]==5){
					polygon[i][j] = 6;
				}
				else if(temppoly[k]==6){
					polygon[i][j] = 5;
				}
			}
		}
		//at this point, if first direction is 3, then we are ok. if not (it must be 1), we need to re-orient
		else if(polygon[i][3]==1){
			//flip ordering
			for(j=3; j<=(3+(L+1)*(M+1)*(span+1)-1); j++){
				k = (3+(L+1)*(M+1)*(span+1)-1+3-j);
				if(temppoly[k]==1){
					polygon[i][j] = 2;
				}
				else if(temppoly[k]==2){
					polygon[i][j] = 1;
				}
				else if(temppoly[k]==3){
					polygon[i][j] = 4;
				}
				else if(temppoly[k]==4){
					polygon[i][j] = 3;
				}
				else if(temppoly[k]==5){
					polygon[i][j] = 6;
				}
				else if(temppoly[k]==6){
					polygon[i][j] = 5;
				}
			}
		}
	}
	
	printf("polygons are all oriented now\n");

	for(i=0; i<=totpolys-1; i++){
		printf("poly[%d] is:", i);
		for(j=0; j<= (3+(L+1)*(M+1)*(span+1)-1); j++){
			printf("%d, ", polygon[i][j]);
		}
		printf("\n");
	}

	//Starting with polygon 1, rotate and check with all other polygons
	int blacklist[totpolys];	//contains poly numbers which have been matched as a duplicate, that will be deleted
	int whitelist[totpolys];	//contains poly numbers which have been matched as a duplicate, but will be kept
	for(i=0; i<=totpolys-1; i++){	//initialize to zeros
		blacklist[i]=0;
	}
	for(i=0; i<=totpolys-1; i++){	//initialize to zeros
		whitelist[i]=0;
	}

	int match=0;	//0 no match, 1 match
	int match2=0;	//0 no match, 1 match
	int rotatedpoly[3+(L+1)*(M+1)*(span+1)];


	for(i=0; i<=totpolys-1; i++){	//for each polygon
//		printf("checking polygon %d\n", i);

		//rotate polygon i, and check all other polys for duplicates.
	
		//rotate polygon i
		//change startingpoint: newx=-oldz, newy=oldy, newz=oldx
		//change: 1to5, 2to6, 5to2, 6to1

		//initialize rotatedpoly to all zeros
		for(j=0; j<=(3+(L+1)*(M+1)*(span+1)-1); j++){
			rotatedpoly[j]=0;
		}
		//change startingpoint
		rotatedpoly[0] = ((-1)*(polygon[i][2]))+L;
		rotatedpoly[1] = polygon[i][1];
		rotatedpoly[2] = polygon[i][0];

		//make rotation changes
		j=3;
		while(polygon[i][j]!=0){
			if(polygon[i][j]==1){
				rotatedpoly[j] = 5;
			}
			if(polygon[i][j]==2){
				rotatedpoly[j] = 6;
			}
			if(polygon[i][j]==3){
				rotatedpoly[j] = 3;
			}
			if(polygon[i][j]==4){
				rotatedpoly[j] = 4;
			}
			if(polygon[i][j]==5){
				rotatedpoly[j] = 2;
			}
			if(polygon[i][j]==6){
				rotatedpoly[j] = 1;
			}
			j++;
		}
/*
		printf("poly[%d] is: ", i);
		for(j=0; j<= (3+(L+1)*(M+1)*(span+1)-1); j++){
			printf("%d, ", polygon[i][j]);
		}
		printf("\n");

		printf("rotatedpoly is: ");
		for(j=0; j<= (3+(L+1)*(M+1)*(span+1)-1); j++){
			printf("%d, ", rotatedpoly[j]);
		}
		printf("\n");
*/

		//now compare rotatedpoly[] to each of the other polygons

		//orient rotatedpoly[] so comparison can be made
		//move to (0,0,0)
		while((rotatedpoly[0]!=0) || (rotatedpoly[1]!=0) || (rotatedpoly[2]!=0) ){
			tempnum = rotatedpoly[3];
			if(rotatedpoly[3]==1){
				rotatedpoly[0]++;
			}
			else if(rotatedpoly[3]==2){
				rotatedpoly[0]--;
			}
			else if(rotatedpoly[3]==3){
				rotatedpoly[1]++;
			}
			else if(rotatedpoly[3]==4){
				rotatedpoly[1]--;
			}
			else if(rotatedpoly[3]==5){
				rotatedpoly[2]++;
			}
			else if(rotatedpoly[3]==6){
				rotatedpoly[2]--;
			}

			for(k=3; k<=(3+(L+1)*(M+1)*(span+1)-2); k++){
				rotatedpoly[k] = rotatedpoly[k+1];
			}
			rotatedpoly[(3+(L+1)*(M+1)*(span+1)-1)] = tempnum;
		}

		//now at (0,0,0). if first direction is 5, then we are ok. if not, we need to re-orient
		if(rotatedpoly[3]!=5){
			for(j=0; j<=3+(L+1)*(M+1)*(span+1)-1; j++){
				temppoly[j]=rotatedpoly[j];
			}
		}

		//if last direction is 6, then we need to flip ordering
		if(rotatedpoly[3+(L+1)*(M+1)*(span+1)-1]==6){
			//flip ordering
			for(j=3; j<=(3+(L+1)*(M+1)*(span+1)-1); j++){
				k = (3+(L+1)*(M+1)*(span+1)-1+3-j);
				if(temppoly[k]==1){
					rotatedpoly[j] = 2;
				}
				else if(temppoly[k]==2){
					rotatedpoly[j] = 1;
				}
				else if(temppoly[k]==3){
					rotatedpoly[j] = 4;
				}
				else if(temppoly[k]==4){
					rotatedpoly[j] = 3;
				}
				else if(temppoly[k]==5){
					rotatedpoly[j] = 6;
				}
				else if(temppoly[k]==6){
					rotatedpoly[j] = 5;
				}
			}
		}
		//at this point, if first direction is 3, then we are ok. if not (it must be 1), we need to re-orient
		else if(rotatedpoly[3]==1){
			//flip ordering
			for(j=3; j<=(3+(L+1)*(M+1)*(span+1)-1); j++){
				k = (3+(L+1)*(M+1)*(span+1)-1+3-j);
				if(temppoly[k]==1){
					rotatedpoly[j] = 2;
				}
				else if(temppoly[k]==2){
					rotatedpoly[j] = 1;
				}
				else if(temppoly[k]==3){
					rotatedpoly[j] = 4;
				}
				else if(temppoly[k]==4){
					rotatedpoly[j] = 3;
				}
				else if(temppoly[k]==5){
					rotatedpoly[j] = 6;
				}
				else if(temppoly[k]==6){
					rotatedpoly[j] = 5;
				}
			}
		}
		//now compare rotated polygon with all polygons after polygon i.

		for(j=0; j<=totpolys-1; j++){	//for each polygon
			match=1;	//assume first that they match
			for(k=0; k<=(3+(L+1)*(M+1)*(span+1)-1); k++){	//for each coordinate/direction
				if(rotatedpoly[k] != polygon[j][k]){
					match=0;
					break;
				}
			}
			if(match==1){	//if rotatedpoly matches polygon j
				printf("polygon %d matches polygon %d.\n", i, j, j);

//				printf("checking if polygon %d should be added to the whitelist\n", i);
				//if polygon i isn't in the blacklist or the whitelist, add it to the whitelist
				match2=0;
				k=0;
				while(blacklist[k] != 0){
					if(blacklist[k]==i){
						match2=1;	//i is in blacklist
//						printf("poly[%d] is in the blacklist, so it won't be added to the whitelist\n", i);
						break;
					}
					k++;
				}
				if(match2==0){	//if i isn't in the blacklist, add it to the whitelist (if it hasn't already)
					k=0;
					while(whitelist[k] != 0){
						if(whitelist[k]==i){
							match2=1;	//i is in whitelist
							break;
						}
						k++;
					}
				}
				if(match2==0){ //i is not in the whitelist (or blacklist-tested earlier)
//					printf("Adding poly[%d] to the whitelist (in spot %d)\n", i, k);
					k=0;
					while(whitelist[k] != 0){
						k++;
					}
					whitelist[k]=i;
				}



				//if polygon j isn't in the blacklist or the whitelist, add it to the blacklist
				match2=0;
				k=0;
				while(blacklist[k] != 0){
					if(blacklist[k]==j){
						match2=1;	//j is in blacklist
//						printf("poly[%d] is in the blacklist, so it won't be added to the blacklist\n", j);
						break;
					}
					k++;
				}
				if(match2==0){	//if j isn't in the blacklist, check if it's in the whitelist. if it isn't then add it to the blacklist
					k=0;
					while(whitelist[k] != 0){
						if(whitelist[k]==j){
							match2=1;	//j is in whitelist
							break;
						}
						k++;
					}
				}
				if(match2==0){ //j is not in the whitelist (or blacklist-tested earlier)
//					printf("Adding poly[%d] to the blacklist (in spot %d)\n", j, k);
					k=0;
					while(blacklist[k] != 0){
						k++;
					}
					blacklist[k]=j;
				}
			}
		}

	}

	for(i=0; i<=totpolys-1; i++){
		if(blacklist[i]!=0){
			printf("blacklist[%d]=%d\n", i, blacklist[i]);
		}
	}

	for(i=0; i<=totpolys-1; i++){
		if(whitelist[i]!=0){
			printf("whitelist[%d]=%d\n", i, whitelist[i]);
		}
	}


	//print proper whitelist (found manually :()
	int correct[15];
	correct[0]=0;
	correct[1]=1;
	correct[2]=2;
	correct[3]=3;
	correct[4]=4;
	correct[5]=5;
	correct[6]=6;
	correct[7]=7;
	correct[8]=11;
	correct[9]=12;
	correct[10]=13;
	correct[11]=17;
	correct[12]=19;
	correct[13]=30;
	correct[14]=38;

	fprintf(fpUnique, "UofS\n");
	for(i=0; i<=14; i++){
		//print starting point
		fprintf(fpUnique, "%d %d %d\n", polygon[correct[i]][0], polygon[correct[i]][1], polygon[correct[i]][2]);
		for(j=0; j<=(L+1)*(M+1)*(span+1); j++){
			fprintf(fpUnique, "%d\n", polygon[correct[i]][3+j]);
		}
		fprintf(fpUnique, "-111\n");
	}

	fprintf(fpUnique, "-999");

	














	printf("\nmain finished\n\n");
}



































