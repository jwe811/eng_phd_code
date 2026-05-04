int search(int prevdirec, int curx, int cury, int endhingenumber) {
/*
This function is given a endhinge number (which defines the endhinge in question). It is also given curx and cury, the point
where the walk is at. It must find which direction comes next. In some cases (if prevdirec != 0), there are 2 directions, which it can possibly go.
In this case, it must not return prevdirec. Also, if there is only one direction, and it's the same as prevdirec, then it must return 0.
*/

	int numdirections=0;	//the number of directions that we found.
	int oned=0, twod=0;	//will hold the two possible directions

	currentendhinge[endhingenumber]=firstendhinge[endhingenumber];	//point to beginning of list of edges
	
	//find the edge (or edges) with a point=(curx,cury)
	while(currentendhinge[endhingenumber] != NULL){	//this will loop through all edges in the endhinge.
		if((*currentendhinge[endhingenumber]).edge.f.x == curx && (*currentendhinge[endhingenumber]).edge.f.y == cury){	//if first point in edge matches
			numdirections++;
			if((*currentendhinge[endhingenumber]).edge.f.x - (*currentendhinge[endhingenumber]).edge.l.x == -1	){	//then we went posx
				if(numdirections==1){
					oned=1;
				}
				else if(numdirections==2){
					twod=1;
				}

			}
			else if((*currentendhinge[endhingenumber]).edge.f.x - (*currentendhinge[endhingenumber]).edge.l.x == 1	){	//then we went negx
				if(numdirections==1){
					oned=2;
				}
				else if(numdirections==2){
					twod=2;
				}


			}
			else if((*currentendhinge[endhingenumber]).edge.f.y - (*currentendhinge[endhingenumber]).edge.l.y == -1	){	//then we went posy
				if(numdirections==1){
					oned=3;
				}
				else if(numdirections==2){
					twod=3;
				}


			}
			else if((*currentendhinge[endhingenumber]).edge.f.y - (*currentendhinge[endhingenumber]).edge.l.y == 1	){	//then we went negy
				if(numdirections==1){
					oned=4;
				}
				else if(numdirections==2){
					twod=4;
				}


			}
		}
		else if((*currentendhinge[endhingenumber]).edge.l.x == curx && (*currentendhinge[endhingenumber]).edge.l.y == cury){	//if second point in edge matches
			numdirections++;
			if((*currentendhinge[endhingenumber]).edge.f.x - (*currentendhinge[endhingenumber]).edge.l.x == 1	){	//then we went posx
				if(numdirections==1){
					oned=1;
				}
				else if(numdirections==2){
					twod=1;
				}

			}
			else if((*currentendhinge[endhingenumber]).edge.f.x - (*currentendhinge[endhingenumber]).edge.l.x == -1	){	//then we went negy
				if(numdirections==1){
					oned=2;
				}
				else if(numdirections==2){
					twod=2;
				}


			}
			else if((*currentendhinge[endhingenumber]).edge.f.y - (*currentendhinge[endhingenumber]).edge.l.y == 1	){	//then we went posy
				if(numdirections==1){
					oned=3;
				}
				else if(numdirections==2){
					twod=3;
				}


			}
			else if((*currentendhinge[endhingenumber]).edge.f.y - (*currentendhinge[endhingenumber]).edge.l.y == -1	){	//then we went negy
				if(numdirections==1){
					oned=4;
				}
				else if(numdirections==2){
					twod=4;
				}


			}
		}
		currentendhinge[endhingenumber] = (*currentendhinge[endhingenumber]).nextendhinge;
	}

	//so at this point, we should have a direction in oned, and maybe a direction in twod
	//need to return the appropriate one.

	if(prevdirec==0){
//		printf("no prevdirec (should only be oned) oned=%d, twod=%d\n", oned, twod);
		return oned;
	}
	else{
//		printf("we have a prevdirec.\n");
		if(oned==1 && prevdirec!=2){
			return oned;
		}
		else if(oned==2 && prevdirec!=1){
			return oned;
		}
		else if(oned==3 && prevdirec!=4){
			return oned;
		}
		else if(oned==4 && prevdirec!=3){
			return oned;
		}
		else{
			return twod;
		}
	}

	return 0;

}









