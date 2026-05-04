void recursivegenerate(int curspan, int printablewalk[], struct hinge_span* tspantoadd, int firstsec, int x, int y){

/*this is the recursive function that will loop through all intermediate 2-spans and create all possible polygons.
Will also have to add the right endhinges at this point and print to file here (i think)
curspan is the current span of the polygon. printablewalk is the walk up to this point (before adding the next 2-span)
tspantoadd is the 2span that comes next in the polygon. the walk in this 2-span must be added and then all of the following
2-spans must be looped through (using recursion)
*/

	int memorywalk[pwalksize];	//this will remember what printablewalk is after tspantoadd is added.
	int i, j;
	int nextsection;
	int curx, cury;
	int prevdirec, direcsadded, newdirec;
	int temp, temp2;
	int counter;
	
	//add tspantoadd's walk to printablewalk[]
//	printf("adding a 2-span to printablewalk (sections %d and %d make up to 2-span)\n", firstsec, tspantoadd->outorder);
	addtowalk(printablewalk, tspantoadd->walk);
	curspan++;
//	printf("2-span added to printablewalk. curspan is now %d\n", curspan);

	//copy printablewalk into memorywalk
	for(i=0; i<=pwalksize-1; i++){
		memorywalk[i] = printablewalk[i];
	}

	if(curspan<totalspan){	//need to add more intermediate 2-spans
//		printf("curspan<totalspan (%d<%d). adding another 2-span (section %d and ?)\n", curspan, totalspan, tspantoadd->outorder);
		nextsection = tspantoadd->outorder;
		struct hinge_span* nexttspan;
		nexttspan = first_hinge_span[nextsection];	//point nexttspan at first
		while(nexttspan->nexthinge != NULL){		//loop through all of nextsection's outorders.
			nexttspan = nexttspan->nexthinge;
//			printf("%dth section=%d\n", curspan+1, nexttspan->outorder);
			recursivegenerate(curspan, printablewalk, nexttspan, nextsection, x, y);
			//reset printablewalk to memorywalk
			for(i=0; i<=pwalksize-1; i++){
				printablewalk[i]=memorywalk[i];
			}
		}
	}
	else{
		//at the proper span. now need to add right endhinges, and then print to file.
//		printf("at proper span (%d=%d), now add right endhinges (if possible)\n", curspan, totalspan);

		/*we have printablewalk[]. need right endhinge information basically.
		the last 2-span in the polygon is tspantoadd. from this, we know the number of right connecting endhinges and 
		the numbers to identify which of those endhinges connect to the right.
		*/
		if(tspantoadd->num_con_endhinges_right>0){	//then we have a valid polygon.
//			printf("there are %d right endhinges that work with this polygon.\n", tspantoadd->num_con_endhinges_right);
			//make appropriate changes to printablewalk to make it a valid polygon walk!
			int curendhinge;
			for(curendhinge=1; curendhinge <= tspantoadd->num_con_endhinges_right; curendhinge++){	//for each endhinge that connects
//				printf("found a proper polygon! printablewalk[] thus far is:\n");

/*				i=0;
				while(printablewalk[i]!=0){
					printf("%d, ", printablewalk[i]);
					i++;
				}
				printf("\n");
*/

				i=0;
				while(printablewalk[i] != 0){
					curx=x;
					cury=y;
					while(printablewalk[i] != -222){
						if(printablewalk[i]==1 || printablewalk[i]==-1){
							curx++;
						}
						else if(printablewalk[i]==2 || printablewalk[i]==-2){
							curx--;
						}
						else if(printablewalk[i]==3 || printablewalk[i]==-3){
							cury++;
						}
						else if(printablewalk[i]==4 || printablewalk[i]==-4){
							cury--;
						}
						i++;
						if(printablewalk[i]==0){
							break;
						}
					}
					//at this point, printablewalk[i] is -222 or printablewalk is done
					if(printablewalk[i]==-222){	//if walk has left to the right, need to insert part of the endhinge
//						printf("need to insert for right endhinge\n");
						prevdirec=0;	//will be zero if this is the first insertion
						direcsadded=0;	//initialized to zero appropriately
						newdirec=1;	//just so it will run atleast once

						while(newdirec!=0){
							//search for edge in endhinge that has curx and cury as a point in it
							//newdirec is the direction that is to be added. newdirec will be zero if there are
//							printf("calling search for the next direction (prevdirec=%d, curx=%d, cury=%d, endhingenumber=%d)\n", prevdirec, curx, cury, tspantoadd->con_endhinges_right[curendhinge-1]);
							newdirec = search(prevdirec, curx, cury, tspantoadd->con_endhinges_right[curendhinge-1]);
//							printf("the search came up with direction %d\n", newdirec);
							//make changes to curx and cury.
							if(newdirec==1){
								curx++;
							}
							else if(newdirec==2){
								curx--;
							}
							else if(newdirec==3){
								cury++;
							}
							else if(newdirec==4){
								cury--;
							}
							//add newdirec into printablewalk if it isn't 0.
							if(newdirec!=0){
								j=i;
								temp=printablewalk[j];
								printablewalk[j]=newdirec;
								j++;
								while(temp != 0){
									temp2=printablewalk[j];
									printablewalk[j] = temp;
									temp=temp2;
									j++;
								}
								//newdirec is now inserted into printablewalk in spot i
								direcsadded++;
								prevdirec=newdirec;
							}
							i++;
						}
						i--;
						//printablewalk[i] should be -222 at this point
//						printf("!!!!!!!!!!!printablewalk[%d]=%d,\n", i, printablewalk[i]);
						//if printablewalk[i+1] is zero, then we are at the end. (actually cant be at the end, since -222 means we left to the right
						//elseif we are not at the end, remove the -333 and all negative entries after it
						if(printablewalk[i+1]==0){
//							printf("we are at the end (THIS SHOULDN'T BE RUNNING!)\n");
						}
						else{//remove -222 and all negative numbers afterwards (until a positive)
//							printf("removing the -222 and all negative numbers afterwards\n");
							j=i;
							counter=0;	//counts the number of entries that will have to be removed after the -222
							while(printablewalk[j+1]<0){
								counter++;
								j++;
							}
//							printf("counter=%d\n", counter);
							//remove the -222 and the next "counter" entries.
							while(printablewalk[i+counter+1]!=0){
								printablewalk[i]=printablewalk[i+counter+1];
								i++;
							}

							//fill everything afterwards with zeros
							j=i+counter+1;
							while(i<=j){
								printablewalk[i]=0;
								i++;
							}
//							printf("finished removing the -222 and others.\n");

							//since this isnt the end, restart loop with i=0
							i=0;
						}
					}
				}

//				printf("GOT TO HERE. POLYGON IS CREATED\n");


				//CREATE A COPY OF printablewalk[], SEND IT TO bfacf_shrink, SIZE OF SHRUNK POLY
				i=0;
				while(printablewalk[i] >0){
					i++;
				}
//				printf("this poly has %d edges\n", i);
				int walkcopy[pwalksize];
				for(j=0; j<=pwalksize-1; j++){
					walkcopy[j] = printablewalk[j];
				}
//				printf("made a copy of the walk\n");
				i=bfacf_shrink(x, y, 0, walkcopy, i);


				if(i>=minedges){
					grand++;
				}



/*				i=0;
				printf("FINAL: (%d, %d, 0) ", x, y);
				while(printablewalk[i] != 0){
					printf("%d, ", printablewalk[i]);
					i++;
				}
				printf("\n");
*/
//				printf("resetting printablewalk to memorywalk\n");
				//reset printablewalk to memorywalk
				for(i=0; i<=pwalksize-1; i++){
					printablewalk[i]=memorywalk[i];
				}


			}
		}
		else{
//			printf("last 2-span doesnt have any right endhinges\n");
		}
	}
















}




















