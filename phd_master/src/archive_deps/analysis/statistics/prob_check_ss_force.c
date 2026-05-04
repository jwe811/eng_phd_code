double prob_check(double Beta,double LmultR, double kappa, double fval) {
//checks if the expected num occurrences of all 2-spans add up to one, without dividing by beta
	
#if !defined(connectivity_inv)	/* if connectivity_inv is not defined declare it here */

	double			connectivity_inv=kappa;

#endif

	double totsum=0;
	int i, j;

	for(i=1; i<=max_sections; i++){
		if(num_outsections[i]>0){	//for each valid section
			totsum += R_Evector[0][i] * L_Evector[0][i] / LmultR;
			printf("Added Section %d's prob, now totsum= %f\n", i, totsum);
		}
	}
//	printf("LmultR=%f\n", LmultR);
	return(totsum);
}
