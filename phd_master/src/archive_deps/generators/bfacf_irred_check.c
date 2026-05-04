//this function takes a polygon in the tube (L,M) (x and y) and tries to shrink it using only p(0) and p(-2) moves (while staying in the tube). Returns 1 if it can be shrunk. 0 otherwise
//CHANGES POLYGON


int bfacf_irred_check(int x, int y, int z, int polygon[], int length){
//	printf("bfacf_irred_check() is running\n");
	int maxiter=10*length;

	int i, j, k;
	int chosenedge;
	int one, two, three;	//INDICES!! by convention, one<three (for deleting purposes, delete three then one)
	int successshrink;
	int uno, dos;
	int ninetydegree;
	int x1, y1, z1;
	int x2, y2, z2;
	int x3, y3, z3;
	int x4, y4, z4;
	int occupied;
	int temp;
	int performedzeromove;
	int un, deux, trois;
	int options;
	int randnum;

	time_t t;
	srand((unsigned) time(&t));

	//remove -111 at end of polygon
	i=0;
	while(polygon[i]>0){
		i++;
	}

	//check
	if(i!=length){
		printf("problemo\n. quitting\n");
		exit(1);
	}


	//ALGORITHM to shrink:	NOTE: must be careful when picking the first or last edge.
	//choose random edge.
	//attempt shrink move
	//if can't shrink, attempt zero move (must check self-avoidance)
	
	for(i=0; i<=maxiter-1; i++){	//at each iteration

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
			return 1;//can shrink
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

				if(x2>L || y2>M){
					occupied=1;
				}
				if(occupied==0){
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

	}

	return 0; //couldn't shrink

}


removeedge(int k, int polygon[]){
	//removes edge with index k from polygon[] and makes appropriate adjustments
//	printf("removing edge %d (polygon[%d]=%d)\n", k, k, polygon[k]);
	while(polygon[k] !=0){
		polygon[k]=polygon[k+1];
		k++;
	}
}

addedge(int direc, int k, int polygon[], int length){
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
































