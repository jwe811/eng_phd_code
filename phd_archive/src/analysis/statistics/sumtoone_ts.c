double sumtoone(double Beta,double LmultR, double kappa) {
//checks if the expected num occurrences of all 2-spans add up to one, without dividing by beta
	
#if !defined(connectivity_inv)	/* if connectivity_inv is not defined declare it here */

	double			connectivity_inv=kappa;

#endif

	int i;
	unsigned long int total_span;
	unsigned long int section_num;
	unsigned long int tspan_num = 1;
	unsigned long int first_tspan = 1;
	double			Enum_Enum;	/* Warning: enum is a keyword in C */

	double totsum = 0;	//will keep track of the probability sum
	double cursum;

	total_span=2;		//only looking at patterns that are of span 2

	for(tspan_num=1; tspan_num<=max_tspans; tspan_num++){	//for each tspan number
		first_tspan = tspan_num; //since we're only looking at patterns that are two spans right now
//		Enum_Enum = (L_Evector[0][first_tspan]) * (R_Evector[0][tspan_num]) / LmultR;

		//find the first section and which outsection from that first section gives that tspan number
		for(section_num = 1; section_num<=max_sections; section_num++){
			for(i=1; i<=num_outsections[section_num]; i++){
				if(tspans_nrr[section_num][i] == tspan_num){
					goto loopdone;
				}
			}
		}
		loopdone:



		Enum_Enum = (L_Evector[0][first_tspan]) * (R_Evector[0][tspan_num]) / LmultR;
		cursum = Enum_Enum;
		totsum += cursum;

		printf("tspan #%d, L_Evector[0][%d])=%f,  (R_Evector[0][%d])=%f, LmultR = %f\n", tspan_num, first_tspan, L_Evector[0][first_tspan], tspan_num, R_Evector[0][tspan_num], LmultR);
		printf("section_num=%d, outnum=%d\n", section_num, i);
		printf("Adding tspan #%d, added %f to totsum, totsum now = %f (edges=%d)\n", tspans_nrr[section_num][i], cursum, totsum, tspans_edges[section_num][i]);

	}

	return(totsum);
}











