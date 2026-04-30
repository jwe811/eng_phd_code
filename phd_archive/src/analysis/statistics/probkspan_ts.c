void probkspan(double Beta, double LmultR, double kappa)
{
	printf("probkspan() (2-span version) is running!\n");

	//CHANGEABLE THINGS HERE: (ctrl+F changeable to find other changeable things)
	//total number of files that need to be read in.
	int totalfilenum = 1;
	int totalspan = 6;	//span of polygons being looked at

	double connectivity_inv=kappa;

	int i, j, k;
	int firstsec=0;
	int secondsec=0;
	int secondlastsec=0;
	int lastsec=0;
	float expoccur=0;
	float prob=0;
	float expoccur_sum=0;
	float prob_sum=0;
	float Enum_Enum=0;
	int polycount=0;

	int Edges[7];
//	Edges[0] = #edges in first section
//	Edges[1] = #edges in second section
//	Edges[2] = #edges in 2nd last section
//	Edges[3] = #edges in last section
//	Edges[4] = #edges in first section and hinge of first tspan
//	Edges[5] = #edges in first section and hinge of last tspan
//	Edges[6] = #edges in k-span minus #edges in first section and hinge of last tspan (Edges[5])

	int firsttspannum=0;	//first tspan's number
	int lasttspannum=0;		//last tspan's number
	int template[4][vM][vL];	//this template will hold sections 1, 2, 2nd last section, last section.
						//can then copy to ordertemplate when num_section is needed.
	int x, y, z;
	int direc;
	char line[7];

	//point at file and read it
	FILE* fp;
	int filenum = 1;

	char filename[50];


	while(filenum<=totalfilenum){	//while there are files to read in
		//CHANGEABLE (add knottype to filename later)
		sprintf(filename, "CreatorPrimePolysL%dM%dspan%dnum%d_3_1_all.txt", L, M, totalspan, filenum);

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
				fill	int Edges[7]
					Edges[0] = #edges in first section
					Edges[1] = #edges in second section
					Edges[2] = #edges in 2nd last section
					Edges[3] = #edges in last section
					Edges[4] = #edges in first section and hinge of first tspan
					Edges[5] = #edges in first section and hinge of last tspan
					Edges[6] = #edges in k-span (not including last section)
				fill template[2][vM][vL]
					template[0] is first section
					template[1] is second section
					template[2] is 2nd-last section
					template[3] is last section
				find first tspan and second tspan.
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
			Edges[3]=0;
			Edges[4]=0;
			Edges[5]=0;
			Edges[6]=0;

			//reset template[]
			for(i=0; i<=vM-1; i++){
				for(j=0; j<=vL-1; j++){
					template[0][i][j]=0;
					template[1][i][j]=0;
					template[2][i][j]=0;
					template[3][i][j]=0;
				}
			}

			direc=0;
			while(1){	//while still in the polygon

				fgets(line, 7, fp);
				sscanf(line, "%d", &direc);	//get direction

				if(direc==-111){
					break;	//break if end of the polygon
				}

				//is Edges[0, 1, 2, 3, 4, 5] and/or template[0, 1, 2, 3] changed? (will do Edges[6] later)
				if(((z==0) && direc==5 )|| (z==1 && direc==6)){	//edge is in first section
					Edges[0]++;
					Edges[4]++;
					template[0][y][x] = Edges[0];
				}
				else if((z==1 && direc==5) || (z==2 && direc==6)){	//edge is in second section
					Edges[1]++;
					template[1][y][x] = Edges[1];
				}
				else if((z==1) && (direc!=5 && direc!=6)){	//edge is in hinge of the first tspan
					Edges[4]++;
				}
				else if((z==(totalspan-2) && direc==5) || (z==(totalspan-1) && direc==6)){	//edge is in 2nd last section
					Edges[2]++;
					Edges[5]++;
					template[2][y][x] = Edges[2];
				}
				else if((z==(totalspan-1) && direc==5) || (z==(totalspan) && direc==6)){	//edge is in last section
					Edges[3]++;
					template[3][y][x] = Edges[3];
				}
				else if(z==totalspan-1 && (direc!=5 && direc!=6)){
					Edges[5]++;
				}

				//Now, update x, y, z and is Edges[6] changed?
				if(direc==1){
					x=x+1;
					if(z!=0 && z!=totalspan){	//if not in an endhinge
						Edges[6]++;
					}
				}
				else if(direc==2){
					x=x-1;
					if(z!=0 && z!=totalspan){	//if not in an endhinge
						Edges[6]++;
					}
				}
				else if(direc==3){
					y=y+1;
					if(z!=0 && z!=totalspan){	//if not in an endhinge
						Edges[6]++;
					}
				}
				else if(direc==4){
					y=y-1;
					if(z!=0 && z!=totalspan){	//if not in an endhinge
						Edges[6]++;
					}
				}
				else if(direc==5){
					z=z+1;
					if(z!=totalspan){	//if not in last section
						Edges[6]++;
					}
				}
				else if(direc==6){
					z=z-1;
					if(z!=(totalspan-1)){	//if not in last section
						Edges[6]++;
					}
				}
			}

			//Get first tspan number
			//copy template into ordertemplate
			for(i=0; i<=vM-1; i++){
				for(j=0; j<=vL-1; j++){
					ordertemplate[0][i][j] = template[0][i][j];
					ordertemplate[1][i][j] = template[1][i][j];
				}
			}
			firstsec = num_section(0);
			secondsec = num_section(1);
//			printf("firstsec=%d, secondsec=%d, edges=%d\n", firstsec, secondsec, Edges[4]);

			//find tspan number that has firstsec, secondsec, and Edges[4] edges.
			firsttspannum=0;
			for(i=1; i<=num_outsections[firstsec]; i++){
//				printf("%dth outsection of %d is %d (Edges=%d).\n", i, firstsec, tspans_outsection[firstsec][i], tspans_edges[firstsec][i]);
				if((tspans_outsection[firstsec][i] == secondsec) && (Edges[4]==tspans_edges[firstsec][i])){	//if the 2-span has the correct sections numbers and the same edges
//					printf("found. i=%d\n", i);
					firsttspannum = tspans_nrr[firstsec][i];
					break;
				}
			}
			if(firsttspannum==0){
				printf("Didn't find a firsttspan (polycount=%d). firstsec=%d, secondsec=%d, edges=%d. Quitting.\n", polycount, firstsec, secondsec, Edges[4]);
				printsection(firstsec);
				printsection(secondsec);
				exit(1);
			}

			//Get second tspan number
			//copy template into ordertemplate
			for(i=0; i<=vM-1; i++){
				for(j=0; j<=vL-1; j++){
					ordertemplate[0][i][j] = template[2][i][j];
					ordertemplate[1][i][j] = template[3][i][j];
				}
			}
			secondlastsec = num_section(0);
			lastsec = num_section(1);
//			printf("secondlastsec=%d, lastsec=%d, edges=%d\n", secondlastsec, lastsec, Edges[5]);

			//find tspan number that has secondlastsec, lastsec, and Edges[5] edges.
			lasttspannum=0;
			for(i=1; i<=num_outsections[secondlastsec]; i++){
//				printf("%dth outsection of %d is %d (Edges=%d).\n", i, secondlastsec, tspans_outsection[secondlastsec][i], tspans_edges[secondlastsec][i]);
				if((tspans_outsection[secondlastsec][i] == lastsec) && (Edges[5]==tspans_edges[secondlastsec][i])){	//if the 2-span has the correct sections numbers and the same edges
//					printf("found. i=%d\n", i);
					lasttspannum = tspans_nrr[secondlastsec][i];
					break;
				}
			}
			if(lasttspannum==0){
				printf("Didn't find a lasttspan (polycount=%d). secondlastsec=%d, lastsec=%d, edges=%d. Quitting.\n", polycount, secondlastsec, lastsec, Edges[5]);
				printsection(secondlastsec);
				printsection(lastsec);
				exit(1);
			}

//			printf("firsttspannum=%d\n", firsttspannum);
//			printf("lasttspannum=%d\n", lasttspannum);

			Enum_Enum = (L_Evector[0][firsttspannum]) * (R_Evector[0][lasttspannum]) / LmultR;
//			printf("Enum_Enum=%f=%e\n", Enum_Enum, Enum_Enum);


			expoccur = Enum_Enum * pow(connectivity_inv, Edges[6]-Edges[5]) * exp(fval*(totalspan-2))/ Beta;
			prob = Enum_Enum * pow(connectivity_inv, Edges[6]-Edges[5]) * exp(fval*(totalspan-2));

			expoccur_sum += expoccur;
			prob_sum += prob;

//			printf("expoccur=%f=%e\n", expoccur, expoccur);
//			printf("expoccur_sum=%f=%e\n", expoccur_sum, expoccur_sum);

//			printf("prob=%f=%e\n", prob, prob);
//			printf("prob_sum=%f=%e\n\n", prob_sum, prob_sum);
		}


		filenum++;
	}

	printf("polycount=%d\n", polycount);
	printf("sum of all of the expected number of occurrences per edge of this file is: %f=%e\n", expoccur_sum, expoccur_sum);
	printf("sum of all of the probabilities of this file is: %f=%e\n", prob_sum, prob_sum);
	printf("probkspan() is FINISHED!\n");
}










