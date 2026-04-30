double transprobcheck(double Beta,double LmultR,double kappa) {
	//checks if the transition probabilites for each section add up to one, when summed over all of the connecting sections
	#if !defined(connectivity_inv)	/* if connectivity_inv is not defined declare it here */
		double			connectivity_inv=kappa;
	#endif

	double totsum;
	int i, j, k;
	int secondsec;
	for(i=1; i<=max_sections; i++){
		if(num_outsections[i]>0){
			for(j=1; j<=num_outsections[i]; j++){
				//we have a tspan.
				secondsec = tspans_outsection[i][j];
				if(tspans_nrr[i][j]==1){
					printf("Looking at tspan #%d, made of section %d to %d\n", tspans_nrr[i][j], i, secondsec);
				}
				totsum = 0;
				//look at all tspans that this tspan connects to
				for(k=1; k<=num_outsections[secondsec]; k++){
					totsum += pow(connectivity_inv, tspans_edges[i][j]) * exp(fval) * R_Evector[0][tspans_nrr[secondsec][k]] / R_Evector[0][tspans_nrr[i][j]];
//					if(tspans_nrr[i][j]==1){
//						printf("tspan #%d connected to tspan #%d, (sections: %d to %d to %d). ", tspans_nrr[i][j], tspans_nrr[secondsec][k], i, secondsec, tspans_outsection[secondsec][k]);
//						printf("Added %f to totsum. totsum is now %f\n", pow(connectivity_inv, tspans_edges[i][j]) * exp(fval) * R_Evector[0][tspans_nrr[secondsec][k]] / R_Evector[0][tspans_nrr[i][j]], totsum);
//					}
				}
				printf("tspan %d connects to %d other tspans. Probabilities add up to %f\n", tspans_nrr[i][j], num_outsections[secondsec], totsum);
			}

		}
	}
	return(totsum);
}











