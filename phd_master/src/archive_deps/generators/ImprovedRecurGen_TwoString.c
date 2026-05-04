void recursivegenerate(int curspan, int printablewalk[], struct hinge_span* tspantoadd, int firstsec, int y, int z){

/*this is the recursive function that will loop through all intermediate 2-spans and create all possible polygons.
Will also have to add the right endhinges at this point and print to file here (i think)
curspan is the current span of the polygon. printablewalk is the walk up to this point (before adding the next 2-span)
tspantoadd is the 2span that comes next in the polygon. the walk in this 2-span must be added and then all of the following
2-spans must be looped through (using recursion)
*/

	int memorywalk[pwalksize];	//this will remember what printablewalk is after tspantoadd is added.
	int i, j;
	int nextsection;
	int lastsection;
	struct endhinge* endhingetoadd;
	int index, temp1, temp2;
	
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
		if( test2secs(nextsection)==0 ){	//can't be a 2-string
			struct hinge_span* nexttspan;
			nexttspan = first_hinge_span[nextsection];	//point nexttspan at first
			while(nexttspan->nexthinge != NULL){		//loop through all of nextsection's outorders.
				nexttspan = nexttspan->nexthinge;
	//			printf("%dth section=%d\n", curspan+1, nexttspan->outorder);
				recursivegenerate(curspan, printablewalk, nexttspan, nextsection, y, z);
				//reset printablewalk to memorywalk
				for(i=0; i<=pwalksize-1; i++){
					printablewalk[i]=memorywalk[i];
				}
			}
		}
	}
	else{
		//close off polygon with ONLY ONE possible endhinges.
		int foundrighthinge=0;
		//Need final section number.
		lastsection = tspantoadd->outorder;
		if(test2secs(lastsection)==1){
			endhingetoadd = firstendhinge[lastsection];
			while(endhingetoadd->nextendhinge!=NULL && foundrighthinge==0){
				endhingetoadd = endhingetoadd->nextendhinge;
				if(endhingetoadd->side == 1){
					foundrighthinge=1;
					//add in right endhinge
	/*				printf("Adding right endhinge to poly. Printablewalk is: ");
					i=0;
					while(printablewalk[i]!=0){
						printf("%d ", printablewalk[i]);
						i++;
					}
					printf("\n");
					printf("endhinge is: ");
					i=0;
					while(endhingetoadd->endwalk[i]!=0){
						printf("%d ", endhingetoadd->endwalk[i]);
						i++;
					}
					printf("\n");
	*/






					i=0;
					j=0;
					while(printablewalk[i]!=0){
						while(printablewalk[i]!=-222){
							i++;
							if(printablewalk[i]==0){
								break;
							}
						}
						if(printablewalk[i]==-222){
							//insert right endhinge walk here
	//						printf("inserting into spot %d\n", i);
							while((*endhingetoadd).endwalk[j]!=-333){
								if((*endhingetoadd).endwalk[j] == 3 ||(*endhingetoadd).endwalk[j] == 4 ||(*endhingetoadd).endwalk[j] == 5 ||(*endhingetoadd).endwalk[j] == 6){
									//insert (*endhingetoadd).endwalk[j] into printable walk[i]
									index=i;
	//								printf("RUNNING. i=%d, index=%d, j=%d\n", i, index, j);
									//shift printablewalk one over
									temp1=printablewalk[index];
									while(printablewalk[index]!=0){
										temp2=printablewalk[index+1];
										printablewalk[index+1]=temp1;
										temp1=temp2;
										index++;
									}
									//add nextwalk component
									printablewalk[i]=(*endhingetoadd).endwalk[j];
									j++;
									i++;
								}
								else{
									j++;
								}
							}
	//						printf("Done inserting part of the endhinge. Can now remove the negative parts of printablewalk. i=%d, j=%d\n", i, j);
							while(printablewalk[i]<0){
	//							printf("deleting printablewalk[%d]=%d\n", i, printablewalk[i]);
								index=i;
								while(printablewalk[index]!=0){
									printablewalk[index] = printablewalk[index+1];
									index++;
								}
							}
							j++;
							i++;
						}
					}
	/*
					printf("printable walk is now: ");
					i=0;
					while(printablewalk[i]!=0){
						printf("%d ", printablewalk[i]);
						i++;
					}
					printf("\n");
	*/




		//			printf("GOT TO HERE. POLYGON IS CREATED\n");


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
					i=bfacf_shrink(0, y, z, walkcopy, i);

	//				printf("i=%d, minedges=%d\n", i, minedges);
					if(i>=minedges){
						//we have a possibly knotted poly
						if(num_uptomill==999999){
							num_millions++;
							num_uptomill=0;
						}
						else{
							num_uptomill++;
						}

						//print printablewalk to file!
						fprintf(fp, "0 %d %d\n", y, z);	//starting point
						i=0;
						while(printablewalk[i] != 0){
							fprintf(fp, "%d\n", printablewalk[i]);
							i++;
						}
						fprintf(fp, "-111\n");	//end of poly
						grand++;
						filetotal++;;

						//check if we need to create another file (cause too many polys)
						if(filetotal>=maxpolys){	//then we need to close up the current file, and create a new one to write to.
							fprintf(fp, "-999\n");
							fclose(fp);
							filenum++;
							sprintf(filename, "ImpCreatorTwoStrings_PKshrunkL%dM%dspan%dnum%d.txt", passedL, passedM, totalspan, filenum);
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
		/*
					i=0;
					printf("FINAL: (0, %d, %d) ", y, z);
					while(printablewalk[i] != 0){
						printf("%d, ", printablewalk[i]);
						i++;
					}
					printf("\n\n");
	*/
		//			printf("resetting printablewalk to memorywalk\n");
					//reset printablewalk to memorywalk
					for(i=0; i<=pwalksize-1; i++){
						printablewalk[i]=memorywalk[i];
					}
				}
			}
		}
	}





}




















