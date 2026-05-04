double prob_check(double Beta,double LmultR, double kappa) {
//checks if the expected num occurrences of all 2-spans add up to one, without dividing by beta
	double totsum=0;
	int i;

	for(i=1; i<=max_sections; i++){
		if(num_outsections[i]>0){	//for each valid section
			totsum += R_Evector[0][i] * L_Evector[0][i] / LmultR;
			printf("Added Section %d's prob, now totsum= %f\n", i, totsum);
		}
	}
//	printf("LmultR=%f\n", LmultR);
	return(totsum);
}
