void groupings(double Beta,double LmultR,double kappa) {
	/*
	This functions looks at all 2-spans, and looks at their expected number of occurrences. It groups 2-spans
	that have the same the expected number of occurrences together. It will output the different groupings
	to the user.
	*/
	//only works in TS mode

	printf("STARTING groupings() (exactly equal)\n");
	printf("There are %d valid 2-spans.\n", max_tspans);

	double exp_occur[max_tspans+1];	//will hold the expected number of occurrences for each 2-span
	int i, j, k, m;

	for(i=1; i<=max_tspans; i++){
		exp_occur[i] = expoccur(Beta, LmultR, kappa, i);
	}

	//group tspans with the same expected number of occurrences together.
	int matches[max_tspans];	//will hold the tspan numbers that already have matches

	for(i=0; i<=max_tspans-1; i++){
		matches[i]=0;				//initialize to zeros
	}
	
	int matchfoundi=0;	//0=match not found, 1=match found.
	int matchfoundj=0;
	int matchesindex=0;

	int firstindextotal=5000;
	int secondindextotal=100;
	int groups[firstindextotal][secondindextotal];		//first index is the "group number", 2nd index contains the tspan number that's in that group
	for(i=0; i<=firstindextotal-1; i++){
		for(j=0; j<=secondindextotal-1; j++){
			groups[i][j]=0;	//initialize groups to zeros
		}
	}
	printf("groups[][] initialized to all zeros\n");
	int numgroups=0;			//number of groups that exist
	int curgroup=0;			//group that's of interest
	int groupi, groupj;

	for(i=1; i<=max_tspans; i++){
		printf("i=%d\n", i);
		for(j=i+1; j<=max_tspans; j++){
//			printf("j=%d\n", j);
			if(exp_occur[i]==exp_occur[j]){	//tspan i matches tspan j
//				printf("tspan %d matches %d\n", i, j);
				//check if matches[] contains i yet; if not, then add it to matches
				//also need to create groups of tspans that have matching expected occurrences
				matchfoundi=0;
				for(k=0; k<=matchesindex-1; k++){
					if(matches[k]==i){
//						printf("tspan %d already is in a group\n", i);
						matchfoundi=1;	//i is in matches[] (i==matches[matchesindex])
						break;
					}
				}

				if(matchfoundi==0){	//need to add i to matches[]
//					printf("adding tspan %d to matches[]\n", i);
					matches[matchesindex]=i;
					matchesindex++;
				}
				else{	//i already has a match, so it's in a group.
//					printf("finding which group tspan %d is in. numgroups=%d.\n", i, numgroups);
					for(curgroup=0; curgroup<=numgroups-1; curgroup++){		//find which group i is in
						for(m=0; m<=secondindextotal-1; m++){
							if(groups[curgroup][m]==i){	//i is in group: curgroup
								groupi=curgroup;		//i is in groupi.
//								printf("tspan %d is in group %d\n", i, groupi);
								goto groupifound;
							}
						}
					}
					groupifound:;
				}
				
				//check if matches contains j yet; if not, then add it to matches
				matchfoundj=0;
				for(k=0; k<=matchesindex-1; k++){
					if(matches[k]==j){
//						printf("tspan %d already is in a group\n", j);
						matchfoundj=1;	//j is in matches[]
						break;
					}
				}
				if(matchfoundj==0){	//need to add j to matches[]
//					printf("adding tspan %d to matches[]\n", j);
					matches[matchesindex]=j;
					matchesindex++;
				}
				else{	//j already has a match, so it's in a group.
//					printf("finding which group tspan %d is in\n", j);
					for(curgroup=0; curgroup<=numgroups-1; curgroup++){		//find which group j is in
						for(m=0; m<=secondindextotal-1; m++){
							if(groups[curgroup][m]==j){	//j is in group: curgroup
								groupj=curgroup;		//j is in groupj.
//								printf("tspan %d is in group %d\n", j, groupj);
								goto groupjfound;
							}
						}
					}
					groupjfound:;
				}
				//4 CASES:
				//1. i is in a group, j isn't
				//2. j is in a group, i isn't
				//3. neither i or j are in a group
				//4. i and j are already in a group (the same group)
//				printf("matchfoundi=%d, matchfoundj=%d\n", matchfoundi, matchfoundj);
				if(matchfoundi==1 && matchfoundj==0){	//1. need to add j to group groupi
//					printf("need to add tspan %d to group %d\n", j, groupi);
					for(m=0; m<=secondindextotal-1; m++){	//find proper 2nd index in groupi to add j
						if(groups[groupi][m]==0){
							groups[groupi][m]=j;
//							printf("added tspan %d to group %d\n", j, groupi);
							break;
						}
					}
				}
				else if(matchfoundi==0 && matchfoundj==1){	//2. need to add i to group groupj
//					printf("need to add tspan %d to group %d\n", i, groupj);
					for(m=0; m<=secondindextotal-1; m++){	//find proper 2nd index in groupj to add i
						if(groups[groupj][m]==0){
							groups[groupj][m]=i;
//							printf("added tspan %d to group %d\n", i, groupj);
							break;
						}
					}
				}
				else if(matchfoundi==0 && matchfoundj==0){	//3. neither i or j are in a group
//					printf("neither i or j are in a group\n");
					groups[numgroups][0]=i;
					groups[numgroups][1]=j;
//					printf("added tpans %d and %d to a new group: %d\n", i, j, numgroups);
					numgroups++;

				}
				else{	//4. i and j are already in a group (the same group)
					//don't need to do anything?
//					printf("tspans %d and %d are already in a group, don't do anything\n", i , j);
				}

			}
		}
	}

	int inmatches=0;	//0 if i is not in matches[], 1 if i is in matches[]
	int counter=0;
	for(i=1; i<=max_tspans; i++){
		inmatches=0;
		for(k=0; k<=matchesindex-1; k++){
			if(i==matches[k]){
				inmatches=1;
				break;
			}
		}
		if(inmatches==0){
//			printf("tspan %d has no matches (exp_occur[%d])=%e\n", i, i, exp_occur[i]);
			counter++;
		}
	}
	printf("%d tspans have matches\n", matchesindex);
	printf("%d tspans have no matches\n", counter);
	printf("There are %d groups\n\n", numgroups);

	int minm=8;
	int maxm=0;
	//print groups
	for(curgroup=0; curgroup<=numgroups-1; curgroup++){
//		printf("group %d: ", curgroup);
		for(m=0; m<=secondindextotal-1; m++){
			if(groups[curgroup][m]!=0){
//				printf("%d, ", groups[curgroup][m]);
			}
			else{
				break;
			}
		}
//		printf("\n");
		if(m<minm){
			minm=m;
		}
		if(m>maxm){
			maxm=m;
		}
	}
	printf("minm=%d\n", minm);
	printf("maxm=%d\n", maxm);

}






