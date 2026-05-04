int test2secs(int sectionnum) {
/* This function takes in a section number and returns 1 if it consists of only 2 edges... 0 otherwise
*/
	int i, j;
	int i1, j1, i2, j2;

	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			ordertemplate[0][i][j]=0;	//initialize ordertemplate to be empty
		}
	}

	//CASE 1: 2edges in section
	for(i1=0; i1<=M; i1++){
		for(j1=0; j1<=L; j1++){
			ordertemplate[0][i1][j1]=1;	//loop through all possible locations for the first edge
			for(i2=0; i2<=M; i2++){
				for(j2=0; j2<=L; j2++){		//loop through all possible locations for the second edge
					if((i1!=i2) || (j1!=j2)){	//that isn't in the location where the first edge is.
						ordertemplate[0][i2][j2]=2;
						if(num_section(0)==sectionnum){
							return 1;
						}
						ordertemplate[0][i2][j2]=0;
					}
				}
			}
			ordertemplate[0][i1][j1]=0;
		}
	}


	
	return 0;
}
