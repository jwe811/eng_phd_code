double sumtoone(double Beta,double LmultR, double kappa) {
//checks if the expected num occurrences of all 2-spans add up to one, without dividing by beta
	
#if !defined(connectivity_inv)	/* if connectivity_inv is not defined declare it here */

	double			connectivity_inv=kappa;

#endif

	double totsum=0;
	int i, j;

	for(i=1; i<=max_sections; i++){
		if(num_outsections[i]>0){
			for(j=1; j<=num_outsections[i]; j++){
				totsum += pow(connectivity_inv, tspans_edges[i][j]) * exp(fval) * R_Evector[0][tspans_outsection[i][j]] * L_Evector[0][i] / LmultR;
				printf("x0=%f, m=%lu, sec_num=%d, out_num=%lu, Ej=%f, Ni=%f, adding:%f, tot=%f\n", connectivity_inv, tspans_edges[i][j], i, tspans_outsection[i][j], R_Evector[0][tspans_outsection[i][j]], L_Evector[0][i], pow(connectivity_inv, tspans_edges[i][j]) * exp(fval)* R_Evector[0][tspans_outsection[i][j]] * L_Evector[0][i] / LmultR, totsum);
				printf("Section %d to %lu, adding %f, tot is now %f\n\n", i, tspans_outsection[i][j], pow(connectivity_inv, tspans_edges[i][j]) *exp(fval) * R_Evector[0][tspans_outsection[i][j]] * L_Evector[0][i] / LmultR, totsum);
			}

		}
	}
	printf("LmultR=%f\n", LmultR);
	return(totsum);
}
