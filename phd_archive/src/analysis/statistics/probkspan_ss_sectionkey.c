//This version uses converts section number into sectionkey number

void probkspan(double Beta, double LmultR, double kappa)
{
	printf("probkspan() (section version) is running!\n");

	//CHANGEABLE THINGS HERE: (ctrl+F changeable to find other changeable things)
	//total number of files that need to be read in.
	int totalfilenum = 1;
	int totalspan = 6;	//span of polygons being looked at

//	double connectivity_inv=kappa;	//comment out if connectivity_inv is defined in HMdocFinal.c

	int i, j;
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
		sprintf(filename, "CreatorAll_InterestPKshrunkPolysL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);

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


			//must reorder ordertemplate
			fillreordertemplate(Edges[0], Edges[1]);

			//set firstsec and lastsec based on the reordertemplate
			firstsec = num_section_norder(0);
			lastsec = num_section_norder(1);

//			printf("OGfirstsec=%d\n", firstsec);
			//change firstsec to it's sectionkey number
			i=1;
			int found=0;
			while(sectionkey[i]!=0){
				if(sectionkey[i]==firstsec){
					found=1;
					firstsec=i;
					break;
				}
				i++;
			}
			if(found==0){	//problem
				printf("didn't find the firstsec section. exitting\n");
				exit(1);
			}
		//	printf("NEWfirstsec=%d\n", firstsec);

//			printf("OGlastsec=%d\n", lastsec);
			//change lastsec to it's sectionkey number
			i=1;
			found=0;
			while(sectionkey[i]!=0){
				if(sectionkey[i]==lastsec){
					found=1;
					lastsec=i;
					break;
				}
				i++;
			}
			if(found==0){	//problem
				printf("didn't find the lastsec section. exitting\n");
				exit(1);
			}
		//	printf("NEWfirstsec=%d\n", firstsec);
			



/*
			printf("firstsec=%d\n", firstsec);
			printsection(firstsec);

			printf("lastsec=%d\n", lastsec);
			printsection(lastsec);

			printf("left[%d] = %e\n", firstsec, L_Evector[0][firstsec]);
			printf("right[%d] = %e\n", lastsec, R_Evector[0][lastsec]);
			printf("LmultR = %e\n", LmultR);
			printf("Edges[2] = %d\n", Edges[2]);
			printf("fval=%f\n", fval);
//			printf("Beta=%f\n", Beta);
			printf("x_0=%f\n", connectivity_inv);

*/
			Enum_Enum = (L_Evector[0][firstsec]) * (R_Evector[0][lastsec]) / LmultR;
/*			printf("Enum_Enum=%f=%e\n", Enum_Enum, Enum_Enum);
			printf("x_0^edges=%e\n", pow(connectivity_inv, Edges[2]));
			printf("e^f(h-1)=%e\n", exp(fval*(totalspan-1)));
			printf("prev * prev=%e\n", pow(connectivity_inv, Edges[2]) * exp(fval*(totalspan-1)));  
*/

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
	printf("polycount=%d\n", polycount);
	printf("sum of all of the expected number of occurrences per edge of this file is: %f=%e\n", expoccur_sum, expoccur_sum);
	printf("sum of all of the probabilities of this file is: %f=%e\n", prob_sum, prob_sum);
	printf("probkspan() is FINISHED!\n");
}










