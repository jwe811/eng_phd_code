void recursivegenerate(int curspan, int printablewalk[], struct hinge_span* tspantoadd, int firstsec, int x, int y){

//THIS VERSION (PossKnots) WILL ONLY PRINT THE POLYGONS WITH >="minedges" EDGES
//THIS VERSION WILL USE bfacf_shrink.c TO FIRST SHRINK POLYGON AND ONLY PRINT WITH >="minedges" EDGES

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
		if((tspantoadd->num_con_endhinges_right>0) && (test2secs(tspantoadd->outorder)==1)){	//then we have a valid polygon.	ONLY WANT POLYGONS WITH 2 EDGES IN THE LAST SECTION
//			printf("there are %d right endhinges that work with this polygon.\n", tspantoadd->num_con_endhinges_right);
			//make appropriate changes to printablewalk to make it a valid polygon walk!
			int curendhinge;
			//v3 CHANGE HERE
			for(curendhinge=1; curendhinge <= 1; curendhinge++){	//for each endhinge that connects
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
							printf("we are at the end (THIS SHOULDN'T BE RUNNING!) quitting\n");
							exit(1);
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

//				printf("i=%d, minedges=%d\n", i, minedges);
				int slipknot=0;
				if(i<minedges){
					//we have an unknot. Check if it's a slipknot.
					//printablewalk[] has the polygon.

					int polygon[pwalksize+3];	//contains the polygon with [0][1][2] containing starting co-ords
					polygon[0] = x;
					polygon[1] = y;
					polygon[2] = 0;
					for(i=3; i<=pwalksize+3-1; i++){
						if(printablewalk[i-3]>=0){
							polygon[i]=printablewalk[i-3];
						}
					}


					int arcone[(L+1)*(M+1)*(totalspan+1)*2];		//will store arcs (made the size of the largest possible arc)
					int arctwo[(L+1)*(M+1)*(totalspan+1)*2];		//first arc[0], arc[1], and arc[2] will be the starting x,y,z of the arc respectively.
															//also need room for if we connect arc outside of tube
					for(i=0; i<=(L+1)*(M+1)*(totalspan+1)*2-1; i++){
						arcone[i]=0;
						arctwo[i]=0;
					}

					//get both arcs.
					int curx=polygon[0];
					int cury=polygon[1];
					int curz=polygon[2];

					//utilizing fact that the first direction of polygon will always be positive z, so that's where arcone will start
					arcone[0]=curx;
					arcone[1]=cury;
					arcone[2]=curz;

					i=3;
					while(curz!=totalspan){	
						arcone[i]=polygon[i];
						if(polygon[i]==1){
							curx++;
						}
						else if(polygon[i]==2){
							curx--;
						}
						else if(polygon[i]==3){
							cury++;
						}
						else if(polygon[i]==4){
							cury--;
						}
						else if(polygon[i]==5){
							curz++;
						}
						else if(polygon[i]==6){
							curz--;
						}
						i++;
					}
/*					printf("arc before adding outside edges: ");
					j=0;
					while(arcone[j] != 0){
						printf("%d ", arcone[j]);
						j++;
					}
					printf("\n");
*/					int save=i;	//save this index for polygon when we do arctwo.

					//at this point, arcone is filled correctly.
					//Now create polygon from arcone[]
					//Reconnect to starting point, going "outside" of the arc
					int otherx = curx;
					int othery = cury;
					int otherz = curz;
					//we are at (otherx, othery, otherz), need to get to (arcone[0], arcone[1], arcone[2])
					for(j=1; j<=L+1; j++){
						arcone[i]=1;
						i++;
						otherx++;
					}
					for(j=1; j<=totalspan; j++){
						arcone[i]=6;
						i++;
						otherz--;
					}
					if(otherz != arcone[2]){	//make sure we're at z=0
						printf("Error when walking around outside the tube. Quitting\n");
						exit(1);
					}
					while(otherx != arcone[0] || othery != arcone[1]){
						if(otherx < arcone[0]){
							arcone[i]=1;
							i++;
							otherx++;
						}
						else if(otherx > arcone[0]){
							arcone[i]=2;
							i++;
							otherx--;
						}
						if(othery < arcone[1]){
							arcone[i]=3;
							i++;
							othery++;
						}
						else if(othery > arcone[1]){
							arcone[i]=4;
							i++;
							othery--;
						}
					}

					//check poly length.
//					printf("arc one all connected up has a length of %d\n", i-1);
					if(i>(L+1)*(M+1)*(totalspan+1)*2){
						printf("problem, arcone is too big. exitting\n");
						exit(1);
					}

					//use bfacf_shrink to see if this poly created from arcone is knotted. if >=24 edges, then we have a slipknot.
					//create copy of polygon to send to bfacf
					for(j=0; j<=pwalksize-1; j++){
						walkcopy[j]=0;
					}

					j=0;
					while(arcone[j+3] != 0){
						walkcopy[j]=arcone[j+3];
						j++;
					}


/*					printf("walkcopy: ");
					j=0;
					while(walkcopy[j]!=0){
						printf("%d ", walkcopy[j]);
						j++;
					}
					printf("\n");
*/
					if(bfacf_shrink(x, y, 0, walkcopy, j) >= 24){	//we have a slipknot
						slipknot=1;
						printf("Found a slipknot (arcone)\n");
					}


					if(slipknot==0){	//then we must repeat for arctwo
						//now need to find where arctwo starts. will start once the direction "6" is activated.
						i=save; //where we left off after finishing arcone
						while(polygon[i]!=6){
							if(polygon[i]==1){
								curx++;
							}
							else if(polygon[i]==2){
								curx--;
							}
							else if(polygon[i]==3){
								cury++;
							}
							else if(polygon[i]==4){
								cury--;
							}
							else if(polygon[i]==5){
								curz++;
							}
							i++;
						}


						//at this point, arc two will start (with the direction 6)
						if(polygon[i] != 6){
							printf("Problem. Exitting.\n");
							exit(1);
						}
						arctwo[0]=curx;
						arctwo[1]=cury;
						arctwo[2]=curz;
						j=3;
						while(curz!=0){
							arctwo[j]=polygon[i];
							if(polygon[i]==1){
								curx++;
							}
							else if(polygon[i]==2){
								curx--;
							}
							else if(polygon[i]==3){
								cury++;
							}
							else if(polygon[i]==4){
								cury--;
							}
							else if(polygon[i]==5){
								curz++;
							}
							else if(polygon[i]==6){
								curz--;
							}
							i++;
							j++;
						}


						//create polygon from arctwo[]
						//now reconnect to starting point, going "outside" of the arc
						otherx = curx;
						othery = cury;
						otherz = curz;
						//we are at (otherx, othery, otherz), need to get to (arctwo[0], arctwo[1], arctwo[2])
						for(i=1; i<=L+1; i++){
							arctwo[j]=1;
							j++;
							otherx++;
						}
						for(i=1; i<=totalspan; i++){
							arctwo[j]=5;
							j++;
							otherz++;
						}
						if(otherz != arctwo[2]){	//make sure we're at z=span
							printf("Error when walking around outside the tube. Quitting\n");
							exit(1);
						}
						while(otherx != arctwo[0] || othery != arctwo[1]){
							if(otherx < arctwo[0]){
								arctwo[j]=1;
								j++;
								otherx++;
							}
							else if(otherx > arctwo[0]){
								arctwo[j]=2;
								j++;
								otherx--;
							}
							if(othery < arctwo[1]){
								arctwo[j]=3;
								j++;
								othery++;
							}
							else if(othery > arctwo[1]){
								arctwo[j]=4;
								j++;
								othery--;
							}
						}

						//arctwo is nofilled appropriately.
						//check poly length.
//						printf("arc two all connected up has a length of %d\n", j-1);
						if(j>(L+1)*(M+1)*(totalspan+1)*2){
							printf("problem, arctwo is too big. exitting\n");
							exit(1);
						}

						//use bfacf_shrink to see if this poly created from arctwo is knotted. if >=24 edges, then we have a slipknot.
						//create copy of polygon to send to bfacf
						for(j=0; j<=pwalksize-1; j++){
							walkcopy[j]=0;
						}

						j=0;
						while(arctwo[j+3] != 0){
							walkcopy[j]=arctwo[j+3];
							j++;
						}


/*						printf("walkcopy: ");
						j=0;
						while(walkcopy[j]!=0){
							printf("%d ", walkcopy[j]);
							j++;
						}
						printf("\n");
*/
						if(bfacf_shrink(arctwo[0], arctwo[1], arctwo[2], walkcopy, j) >= 24){	//we have a slipknot
							slipknot=1;
							printf("Found a slipknot (arctwo)\n");
						}






					}




				}
				if(slipknot==1){
					if(num_uptomill==999999){
						num_millions++;
						num_uptomill=0;
					}
					else{
						num_uptomill++;
					}

					//print printablewalk to file!
					fprintf(fp, "%d %d 0\n", x, y);	//starting point
					i=0;
					while(printablewalk[i] != 0){
						fprintf(fp, "%d\n", printablewalk[i]);
						i++;
					}
					grand++;
					filetotal++;

					//check if we need to create another file (cause too many polys)
					if(filetotal>maxpolys){	//then we need to close up the current file, and create a new one to write to.
						fprintf(fp, "-999\n");
						fclose(fp);
						filenum++;
						sprintf(filename, "SlipknotPolysL%dM%dspan%dnum%d.txt", passedL, passedM, totalspan, filenum);
						fp = fopen(filename, "w");	//create or overwrite the file "filename
						//					printf("THIS IS RUNNING\n");
						if(fp != NULL){
							fprintf(fp, "UofS\n");	//first line in file is always "UofS"
							printf("printed UofS in file '%s'\n", filename);
						}
						else{
							printf("file pointer is pointing to NULL\n");
						}
						printf("finished filling. filenum=%d, filetotal=%d\n", filenum, filetotal);
						filetotal=0;

					}

				}
				else{
					smallcounter++;
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




















