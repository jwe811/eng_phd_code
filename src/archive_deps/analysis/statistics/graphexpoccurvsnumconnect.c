void graphexpoccurvsnumconnect(double Beta,double LmultR,double kappa) {
/*
This function will create a .csv file named graphexpoccurvsnumconnect.csv.  First column will be the number of sections that
can follow the 2-span, the second column will be the expected number of occurrences of that 2-span.
*/

	FILE *fp;	//file pointer. will point to the file that is to be written to.
	fp = fopen("graphexpoccurvsnumconnect.csv", "w");	//create or overwrite the file "graphexpoccurvsnumconnect.csv"

	if(fp != NULL){
		printf("graphexpoccurvsnumconnect.csv is created\n");
	}
	else{
		printf("file pointer is pointing to NULL\n");
	}


	int section, nth_section;
	fprintf(fp, "numconnect expoccur\n");
	for(section=1; section<=max_sections; section++){
		for(nth_section=1; nth_section<=num_outsections[section]; nth_section++){
			fprintf(fp, "%d %f\n", num_outsections[tspans_outsection[section][nth_section]], expoccur(Beta, LmultR, kappa, tspans_nrr[section][nth_section]));
		}
	}


}
