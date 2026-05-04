//one difference for this version is that it will be used during Monte Carlo sampling of long span. So maxiter is reduced.

//this function takes a polygon (and it's starting vertex and length) and tries to shrink it using BFACF moves. returns the size of the shrunk poly.
//CHANGES POLYGON--MAY NO LONGER BE IN TUBE!!

//potential bugs: if polygon grows bigger than "length", it may be larger than polygon[]

int shrinkattempt(int one, int two, int three, int polygon[]);
void addedge(int direc, int k, int polygon[], int length);
void removeedge(int k, int polygon[]);


int bfacf_shrink(int x, int y, int z, int polygon[], int length){
//	printf("bfacf_shrink() is running\n");
	int maxiter=10*length;
	double prob=0.0;	//prob of +2 move

	int i, j;
	int chosenedge;
	int one, two, three;	//INDICES!! by convention, one<three (for deleting purposes, delete three then one)
	int successshrink;
	int uno, dos;
	int ninetydegree;
	int x1, y1, z1;
	int x2, y2, z2;
	int x3, y3, z3;
	int occupied;
	int temp;
	int performedzeromove=0;
	int un, deux, trois;
	int randnum;

	time_t t;
	srand((unsigned) time(&t));


	//ALGORITHM to shrink:	NOTE: must be careful when picking the first or last edge.
	//choose random edge.
	//attempt shrink move
	//if can't shrink, attempt zero move (must check self-avoidance)
	
	for(i=0; i<=maxiter-1; i++){	//at each iteration
		if(length<6){ //We know it's an unknot and can stop.
			break;
		}

		//Choose random edge
		chosenedge = rand() % (length);

//		printf("Edge %d has been chosen (polygon[%d]=%d)\n", chosenedge, chosenedge, polygon[chosenedge]);
		//attempt to shrink
		if(chosenedge==0){	//if first edge is chosen
			one=1;
			two=0;
			three=length-1;			
		}
		else if(chosenedge==length-1){ //if last edge is chosen
			one=0;
			two=length-1;
			three=length-2;
		}
		else{
			one=chosenedge-1;
			two=chosenedge;
			three=chosenedge+1;
		}
		successshrink = shrinkattempt(one, two, three, polygon); //returns 1 if can shrink there, zero otherwise
		if(successshrink==1){	//shrink
//			printf("found a place to shrink\n");
			removeedge(three, polygon);
			length--;
			removeedge(one, polygon);
			length--;
		}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		else{	//try to perform zero move
			performedzeromove=0;
//			printf("Trying to perform zero move\n");
			//must be 90 degree turn
			ninetydegree=0;
			if(chosenedge==length-1){
				uno=chosenedge;
				dos=0;
			}
			else{
				uno=chosenedge;
				dos=chosenedge+1;
			}
			if(polygon[uno]==1 && (polygon[dos]==3 || polygon[dos]==4 || polygon[dos]==5 || polygon[dos]==6)){		//+x and (y or z)
				ninetydegree=1;
			}
			else if(polygon[uno]==2 && (polygon[dos]==3 || polygon[dos]==4 || polygon[dos]==5 || polygon[dos]==6)){		//-x and (y or z)
				ninetydegree=1;
			}
			else if(polygon[uno]==3 && (polygon[dos]==1 || polygon[dos]==2 || polygon[dos]==5 || polygon[dos]==6)){		//+y and (x or z)
				ninetydegree=1;
			}
			else if(polygon[uno]==4 && (polygon[dos]==1 || polygon[dos]==2 || polygon[dos]==5 || polygon[dos]==6)){		//-y and (x or z)
				ninetydegree=1;
			}
			else if(polygon[uno]==5 && (polygon[dos]==3 || polygon[dos]==4 || polygon[dos]==1 || polygon[dos]==2)){		//+z and (y or x)
				ninetydegree=1;
			}
			else if(polygon[uno]==6 && (polygon[dos]==3 || polygon[dos]==4 || polygon[dos]==1 || polygon[dos]==2)){		//-z and (y or x)
				ninetydegree=1;
			}
			else{
//				printf("Edge %d DOES NOT start a 90degree turn.\n", uno);
			}

			if(ninetydegree==1){
//				printf("Edge %d starts a 90degree turn.\n", uno);
				//Check if possible move is empty
				//find "start" vertex location of the edge polygon[uno]
				x1=x;	//will hold this location.
				y1=y;
				z1=z;
//				printf("starting at (%d, %d, %d)\n", x, y, z);
				for(j=0; j<chosenedge; j++){
					if(polygon[j]==1){
						x1++;
					}
					else if(polygon[j]==2){
						x1--;
					}
					else if(polygon[j]==3){
						y1++;
					}
					else if(polygon[j]==4){
						y1--;
					}
					else if(polygon[j]==5){
						z1++;
					}
					else if(polygon[j]==6){
						z1--;
					}
//					printf("Added edge %d (polygon[%d]=%d), we are now at (%d, %d, %d)\n", j, j, polygon[j], x1, y1, z1);
				}	//x1,y1,z1 now is the coordinates of the vertex
//				printf("Before applying edge %d (polygon[%d]=%d), we are at (%d,%d,%d)\n", uno, uno, polygon[uno], x1, y1, z1);
				//move to location which must be empty.
				if(polygon[dos]==1){
					x1++;
				}
				else if(polygon[dos]==2){
					x1--;
				}
				else if(polygon[dos]==3){
					y1++;
				}
				else if(polygon[dos]==4){
					y1--;
				}
				else if(polygon[dos]==5){
					z1++;
				}
				else if(polygon[dos]==6){
					z1--;
				}
//				printf("vertex which must be empty is (%d, %d, %d)\n", x1, y1, z1);
				//Now need to find if it is empty
				x2=x;
				y2=y;
				z2=z;
				occupied=0;

				j=0;
				while(polygon[j]>0){
					if(x1==x2 && (y1==y2 && z1==z2)){
						occupied=1;
//						printf("it is NOT empty. can't perform move\n");
						break;
					}
					if(polygon[j]==1){
						x2++;
					}
					else if(polygon[j]==2){
						x2--;
					}
					else if(polygon[j]==3){
						y2++;
					}
					else if(polygon[j]==4){
						y2--;
					}
					else if(polygon[j]==5){
						z2++;
					}
					else if(polygon[j]==6){
						z2--;
					}
					j++;
				}
				if(occupied==0){	//perform zero move (swap)
//					printf("It is empty. Perform zero move\n");

					temp = polygon[uno];
					polygon[uno]=polygon[dos];
					polygon[dos]=temp;
					performedzeromove=1;
				}
			}			
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(successshrink==0 && performedzeromove==0){
			if(((double)rand()/(double)RAND_MAX) <= prob){	//perform +2 move with probability "prob"
//				printf("See if we can do a +2 move\n");
				//check self-avoidance
				if(chosenedge==0){
					un=length-1;
					deux=chosenedge;
					trois=chosenedge+1;
				}
				else if(chosenedge==(length-1)){
					un=chosenedge-1;
					deux=chosenedge;
					trois=0;
				}
				else{
					un = chosenedge - 1;
					deux = chosenedge;
					trois = chosenedge + 1;
				}

				//find "start" vertex location of the edge polygon[deux]
				x1=x;	//will hold this location.
				y1=y;
				z1=z;
//				printf("starting at (%d, %d, %d)\n", x, y, z);
				for(j=0; j<chosenedge; j++){
					if(polygon[j]==1){
						x1++;
					}
					else if(polygon[j]==2){
						x1--;
					}
					else if(polygon[j]==3){
						y1++;
					}
					else if(polygon[j]==4){
						y1--;
					}
					else if(polygon[j]==5){
						z1++;
					}
					else if(polygon[j]==6){
						z1--;
					}
//					printf("Added edge %d (polygon[%d]=%d), we are now at (%d, %d, %d)\n", j, j, polygon[j], x1, y1, z1);
				}	//x1,y1,z1 now is the coordinates of the vertex
//				printf("Before applying edge %d (polygon[%d]=%d), we are at (%d,%d,%d)\n", deux, deux, polygon[deux], x1, y1, z1);
				x2=x1;
				y2=y1;
				z2=z1;
				if(polygon[deux]==1){
					x2++;
				}
				else if(polygon[deux]==2){
					x2--;
				}
				else if(polygon[deux]==3){
					y2++;
				}
				else if(polygon[deux]==4){
					y2--;
				}
				else if(polygon[deux]==5){
					z2++;
				}
				else if(polygon[deux]==6){
					z2--;
				}
//				printf("After applying edge %d (polygon[%d]=%d), we are at (%d,%d,%d)\n", deux, deux, polygon[deux], x2, y2, z2);
//				printf("polygon[%d]=%d\n", un, polygon[un]);
//				printf("polygon[%d]=%d\n", deux, polygon[deux]);
//				printf("polygon[%d]=%d\n", trois, polygon[trois]);

				//check if needed vertices are unoccupied.
				occupied=0;
				if(polygon[deux]==1 || polygon[deux]==2){	//if in the x-direction, choose between attempting +y,-y,+z,-z
					randnum = (rand() % 3) +1; //random number between 1 and 4 inclusive
					if(randnum==1){ //+y. (x1,y1+1,z1) and (x2, y2+1, z2) must be empty
//						printf("trying to add in +y direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1==x3 && (y1+1==y3 && z1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2==x3 && (y2+1==y3 && z2==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform +y move
//							printf("performing add in +y direction\n");
							addedge(4, deux+1, polygon, length);
							length++;
							addedge(3, deux, polygon, length);
							length++;
						}
					}
					else if (randnum==2){
//						printf("trying to add in -y direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1==x3 && (y1-1==y3 && z1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2==x3 && (y2-1==y3 && z2==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform +y move
//							printf("performing add in -y direction\n");
							addedge(3, deux+1, polygon, length);
							length++;
							addedge(4, deux, polygon, length);
							length++;
						}
					}
					else if (randnum==3){
//						printf("trying to add in +z direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1==x3 && (y1==y3 && z1+1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2==x3 && (y2==y3 && z2+1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform +z move
//							printf("performing add in +z direction\n");
							addedge(6, deux+1, polygon, length);
							length++;
							addedge(5, deux, polygon, length);
							length++;
						}
					}
					else if (randnum==4){
//						printf("trying to add in -z direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1==x3 && (y1==y3 && z1-1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2==x3 && (y2==y3 && z2-1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform +z move
//							printf("performing add in -z direction\n");
							addedge(5, deux+1, polygon, length);
							length++;
							addedge(6, deux, polygon, length);
							length++;
						}
					}
					
				}
				else if(polygon[deux]==3 || polygon[deux]==4){	//if in the y-direction, choose between attempting +x,-x,+z,-z
					randnum = (rand() % 3) +1; //random number between 1 and 4 inclusive
					if(randnum==1){ //+x. (x1+1,y1,z1) and (x2+1, y2, z2) must be empty
//						printf("trying to add in +x direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1+1==x3 && (y1==y3 && z1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2+1==x3 && (y2==y3 && z2==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform +x move
//							printf("performing add in +x direction\n");
							addedge(2, deux+1, polygon, length);
							length++;
							addedge(1, deux, polygon, length);
							length++;
						}
					}
					else if (randnum==2){
//						printf("trying to add in -x direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1-1==x3 && (y1==y3 && z1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2-1==x3 && (y2==y3 && z2==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform -x move
//							printf("performing add in -x direction\n");
							addedge(1, deux+1, polygon, length);
							length++;
							addedge(2, deux, polygon, length);
							length++;
						}
					}
					else if (randnum==3){
//						printf("trying to add in +z direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1==x3 && (y1==y3 && z1+1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2==x3 && (y2==y3 && z2+1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform +z move
//							printf("performing add in +z direction\n");
							addedge(6, deux+1, polygon, length);
							length++;
							addedge(5, deux, polygon, length);
							length++;
						}
					}
					else if (randnum==4){
//						printf("trying to add in -z direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1==x3 && (y1==y3 && z1-1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2==x3 && (y2==y3 && z2-1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform -z move
//							printf("performing add in -z direction\n");
							addedge(5, deux+1, polygon, length);
							length++;
							addedge(6, deux, polygon, length);
							length++;
						}
					}
					
				}
				else if(polygon[deux]==5 || polygon[deux]==6){	//if in the z-direction, choose between attempting +x,-x,+y,-y
					randnum = (rand() % 3) +1; //random number between 1 and 4 inclusive
					if(randnum==1){ //+x. (x1+1,y1,z1) and (x2+1, y2, z2) must be empty
//						printf("trying to add in +x direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1+1==x3 && (y1==y3 && z1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2+1==x3 && (y2==y3 && z2==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform +x move
//							printf("performing add in +x direction\n");
							addedge(2, deux+1, polygon, length);
							length++;
							addedge(1, deux, polygon, length);
							length++;
						}
					}
					else if (randnum==2){
//						printf("trying to add in -x direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1-1==x3 && (y1==y3 && z1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2-1==x3 && (y2==y3 && z2==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform -x move
//							printf("performing add in -x direction\n");
							addedge(1, deux+1, polygon, length);
							length++;
							addedge(2, deux, polygon, length);
							length++;
						}
					}
					else if (randnum==3){
//						printf("trying to add in +y direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1==x3 && (y1+1==y3 && z1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2==x3 && (y2+1==y3 && z2==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform +y move
//							printf("performing add in +y direction\n");
							addedge(4, deux+1, polygon, length);
							length++;
							addedge(3, deux, polygon, length);
							length++;
						}
					}
					else if (randnum==4){
//						printf("trying to add in -y direction\n");
						x3=x;
						y3=y;
						z3=z;
						j=0;
						while(polygon[j]>0){
							if(x1==x3 && (y1-1==y3 && z1==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							else if(x2==x3 && (y2-1==y3 && z2==z3)){
								occupied=1;
		//						printf("it is NOT empty. can't perform move\n");
								break;
							}
							if(polygon[j]==1){
								x3++;
							}
							else if(polygon[j]==2){
								x3--;
							}
							else if(polygon[j]==3){
								y3++;
							}
							else if(polygon[j]==4){
								y3--;
							}
							else if(polygon[j]==5){
								z3++;
							}
							else if(polygon[j]==6){
								z3--;
							}
							j++;
						}
						if(occupied==0){	//perform -y move
//							printf("performing add in -y direction\n");
							addedge(3, deux+1, polygon, length);
							length++;
							addedge(4, deux, polygon, length);
							length++;
						}
					}
					
				}












			}
		}
	}















	i=0;
	while(polygon[i]>0){
		i++;
	}
	if(length==i){
		return i;
	}
	else{
		printf("PROBLEM: LENGTHS DON'T MATCH. QUITTING.\n");
		printf("length=%d, i=%d\n", length, i);
		exit(1);
	}
}


void removeedge(int k, int polygon[]){
	//removes edge with index k from polygon[] and makes appropriate adjustments
//	printf("removing edge %d (polygon[%d]=%d)\n", k, k, polygon[k]);
	while(polygon[k] !=0){
		polygon[k]=polygon[k+1];
		k++;
	}
}

void addedge(int direc, int k, int polygon[], int length){
	//adds direction direc into spot k. shifts everything else forward
//	printf("Adding direc %d into spot %d (polygon[%d]=%d before adding)\n", direc, k, k, polygon[k]);
	while((length-1)>=k){
		polygon[length]=polygon[length-1];
		length--;
	}
	polygon[k]=direc;
}

int shrinkattempt(int one, int two, int three, int polygon[]){
	int canshrink=0;
	if(polygon[one]==5 && polygon[two]==3 && polygon[three]==6){	//+z+y-z
		canshrink=1;
	}
	else if(polygon[one]==5 && polygon[two]==4 && polygon[three]==6){	//+z-y-z
		canshrink=1;
	}
	else if(polygon[one]==5 && polygon[two]==1 && polygon[three]==6){	//+z+x-z
		canshrink=1;
	}
	else if(polygon[one]==5 && polygon[two]==2 && polygon[three]==6){	//+z-x-z
		canshrink=1;
	}
/////
	else if(polygon[one]==6 && polygon[two]==3 && polygon[three]==5){	//-z+y+z
		canshrink=1;
	}
	else if(polygon[one]==6 && polygon[two]==4 && polygon[three]==5){	//-z-y+z
		canshrink=1;
	}
	else if(polygon[one]==6 && polygon[two]==1 && polygon[three]==5){	//-z+x+z
		canshrink=1;
	}
	else if(polygon[one]==6 && polygon[two]==2 && polygon[three]==5){	//-z-x+z
		canshrink=1;
	}
/////
	else if(polygon[one]==1 && polygon[two]==3 && polygon[three]==2){	//+x+y-x
		canshrink=1;
	}
	else if(polygon[one]==1 && polygon[two]==4 && polygon[three]==2){	//+x-y-x
		canshrink=1;
	}
	else if(polygon[one]==1 && polygon[two]==5 && polygon[three]==2){	//+x+z-x
		canshrink=1;
	}
	else if(polygon[one]==1 && polygon[two]==6 && polygon[three]==2){	//+x-z-x
		canshrink=1;
	}
/////
	else if(polygon[one]==2 && polygon[two]==3 && polygon[three]==1){	//-x+y+x
		canshrink=1;
	}
	else if(polygon[one]==2 && polygon[two]==4 && polygon[three]==1){	//-x-y+x
		canshrink=1;
	}
	else if(polygon[one]==2 && polygon[two]==5 && polygon[three]==1){	//-x+z+x
		canshrink=1;
	}
	else if(polygon[one]==2 && polygon[two]==6 && polygon[three]==1){	//-x-z+x
		canshrink=1;
	}
/////
	else if(polygon[one]==3 && polygon[two]==1 && polygon[three]==4){	//+y+x-y
		canshrink=1;
	}
	else if(polygon[one]==3 && polygon[two]==2 && polygon[three]==4){	//+y-x-y
		canshrink=1;
	}
	else if(polygon[one]==3 && polygon[two]==5 && polygon[three]==4){	//+y+z-y
		canshrink=1;
	}
	else if(polygon[one]==3 && polygon[two]==6 && polygon[three]==4){	//+y-z-y
		canshrink=1;
	}
/////
	else if(polygon[one]==4 && polygon[two]==1 && polygon[three]==3){	//-y+x+y
		canshrink=1;
	}
	else if(polygon[one]==4 && polygon[two]==2 && polygon[three]==3){	//-y-x+y
		canshrink=1;
	}
	else if(polygon[one]==4 && polygon[two]==5 && polygon[three]==3){	//-y+z+y
		canshrink=1;
	}
	else if(polygon[one]==4 && polygon[two]==6 && polygon[three]==3){	//-y-z+y
		canshrink=1;
	}
	return canshrink;
}
































