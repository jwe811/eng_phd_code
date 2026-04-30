//this version will print probability value to a file.

double probkspan2(double Beta, double LmultR, double kappa, FILE* fpWrite)
{
	printf("probkspan2() (section version) is running!\n");

	//CHANGEABLE THINGS HERE: (ctrl+F changeable to find other changeable things)
	//total number of files that need to be read in.
	int totalfilenum = 1;
	int totalspan = 7;	//span of polygons being looked at

	double connectivity_inv=kappa;	//comment out if connectivity_inv is defined in HMdocFinal.c

	int i, j, k;
	int firstsec;
	int lastsec;
	float expoccur;
	float prob;
	float expoccur_sum=0;
	float prob_sum=0;
	float Enum_Enum;
	int polycount=0;

	int Edges[3];
//	Edges[0] = #edges in first_sec
//	Edges[1] = #edges in last_sec
//	Edges[2] = #edges in k-span (not including last section)

	int x, y, z;
	int direc;
	char line[7];

	//point at file and read it
	FILE* fp;
	int filenum = 1;

	char filename[50];


	while(filenum<=totalfilenum){	//while there are files to read in
		//CHANGEABLE (add knottype to filename later)
		sprintf(filename, "GlobalTwoStrings_PKshrunkL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);

		fp = fopen(filename, "r");	//open file to read.
		if(fp != NULL){
			printf("Opened '%s'\n", filename);
		}
		else{
			printf("file pointer is pointing to NULL, couldn't open '%s'. Exiting program.\n", filename);
			exit(1);
		}

		fgets(line, 7, fp);	//read first line to get it out of the way
//		printf("first line=%s", line);

		/*ALGORITHM:
			for each polygon, walk through polygon and at each step:
				fill	int Edges[3]
					Edges[0] = #edges in first_sec
					Edges[1] = #edges in last_sec
					Edges[2] = #edges in k-span (not including last section)
				fill ordertemplate[2][vM][vL]
					ordertemplate[0] is first section
					ordertemplate[1] is last section
				perform calculation
				add to overall sum
		*/

		while( fgets(line, 7, fp)!=NULL ){	//while there's a line to read
			sscanf(line, "%d %d %d", &x, &y, &z);	//read in starting point
//			printf("starting point is %d, %d, %d\n", x, y, z);
			if(x==-999){	//end of file
				break;
			}

			polycount++;
			//reset Edges[]
			Edges[0]=0;
			Edges[1]=0;
			Edges[2]=0;

			//reset ordertemplate[]
			for(i=0; i<=vM-1; i++){
				for(j=0; j<=vL-1; j++){
					ordertemplate[0][i][j]=0;
					ordertemplate[1][i][j]=0;
				}
			}

			direc=0;
			while(1){	//while still in the polygon

				fgets(line, 7, fp);
				sscanf(line, "%d", &direc);	//get direction

				if(direc==-111){
					break;	//break if end of the polygon
				}

				//is Edges[0] or Edges[1] and/or ordertemplate[0] or ordertemplate[1] changed? (will do Edges[2] later)
				if(z==0 && direc==5){	//edge is in first section
					Edges[0]++;
					ordertemplate[0][y][x] = Edges[0];
				}
				else if(z==1 && direc==6){	//edge is in first section
					Edges[0]++;
					ordertemplate[0][y][x] = Edges[0];
				}
				else if(z==(totalspan-1) && direc==5){	//edge is in last section
					Edges[1]++;
					ordertemplate[1][y][x] = Edges[1];
				}
				else if(z==(totalspan) && direc==6){	//edge is in last section
					Edges[1]++;
					ordertemplate[1][y][x] = Edges[1];
				}

				//Now, update x, y, z and is Edges[2] changed?
				if(direc==1){
					x=x+1;
					if(z!=0 && z!=totalspan){	//if not in an endhinge
						Edges[2]++;
					}
				}
				else if(direc==2){
					x=x-1;
					if(z!=0 && z!=totalspan){	//if not in an endhinge
						Edges[2]++;
					}
				}
				else if(direc==3){
					y=y+1;
					if(z!=0 && z!=totalspan){	//if not in an endhinge
						Edges[2]++;
					}
				}
				else if(direc==4){
					y=y-1;
					if(z!=0 && z!=totalspan){	//if not in an endhinge
						Edges[2]++;
					}
				}
				else if(direc==5){
					z=z+1;
					if(z!=totalspan){	//if not in last section
						Edges[2]++;
					}
				}
				else if(direc==6){
					z=z-1;
					if(z!=(totalspan-1)){	//if not in last section
						Edges[2]++;
					}
				}
			}

			firstsec = num_section(0);
			lastsec = num_section(1);

//			printf("firstsec=%d\n", firstsec);
//			printsection(firstsec);

//			printf("lastsec=%d\n", lastsec);
//			printsection(lastsec);

			Enum_Enum = (L_Evector[0][firstsec]) * (R_Evector[0][lastsec]) / LmultR;
//			printf("Enum_Enum=%f=%e\n", Enum_Enum, Enum_Enum);


			expoccur = Enum_Enum * pow(connectivity_inv, Edges[2]) * exp(fval*(totalspan-1))/ Beta;
			prob = Enum_Enum * pow(connectivity_inv, Edges[2]) * exp(fval*(totalspan-1));

			expoccur_sum += expoccur;
			prob_sum += prob;

//			printf("expoccur=%f=%e\n", expoccur, expoccur);
//			printf("expoccur_sum=%f=%e\n", expoccur_sum, expoccur_sum);

//			printf("prob=%f=%e\n", prob, prob);
//			printf("prob_sum=%f=%e\n\n", prob_sum, prob_sum);
		}


		filenum++;
	}
//	printf("polycount=%d\n", polycount);
//	printf("sum of all of the expected number of occurrences per edge of this file is: %f=%e\n", expoccur_sum, expoccur_sum);
//	printf("sum of all of the probabilities of this file is: %f=%e\n", prob_sum, prob_sum);

	//changeable
	fprintf(fpWrite, "%e ", prob_sum);

	printf("probkspan2() is FINISHED!\n");
}










