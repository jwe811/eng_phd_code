double transprobcheck(double Beta,double LmultR,double kappa, double fval) {
	//checks if the transition probabilites for each tspan add up to one, when summed over all of the connecting tspans
	double fugacity=connectivity_inv;


	double totsum;
	int i, j;
	for(i=1; i<=max_sections; i++){
		if(num_outsections[i]>0){
			totsum = 0;
			for(j=1; j<=num_outsections[i]; j++){
				totsum += pow(connectivity_inv, tspans_edges[i][j]) * exp(fval) * R_Evector[0][tspans_outsection[i][j]] / R_Evector[0][i];
//				printf("section %d connects to section %d, added %f to totsum, totsum is now %f\n", i, tspans_outsection[i][j], pow(connectivity_inv, tspans_edges[i][j]) * exp(fval) * R_Evector[0][tspans_outsection[i][j]] / R_Evector[0][i], totsum);
//				printf("x0=%f, m=%d, sec_num=%d, out_num=%d, Ej=%f, Ei=%f, adding:%f, tot=%f\n", connectivity_inv, tspans_edges[i][j], i, tspans_outsection[i][j], R_Evector[0][tspans_outsection[i][j]], R_Evector[0][i], pow(connectivity_inv, tspans_edges[i][j]) * R_Evector[0][tspans_outsection[i][j]] / R_Evector[0][i], totsum);
			}
			printf("section num=%d, (connects to %d sections), prob sums to=%f\n", i, num_outsections[i], totsum);
		}
	}
	return(totsum);
}











