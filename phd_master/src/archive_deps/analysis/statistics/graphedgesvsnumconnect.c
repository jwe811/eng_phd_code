void graphedgesvsnumconnect(double Beta,double LmultR,double kappa) {
/*
This function will create a .csv file named graphedgesvsnumconnect.csv.  First column will be the number of sections that
can follow the 2-span, the second column will be the number of edges in the first section and hinge of the 2-span.
*/

	FILE *fp;	//file pointer. will point to the file that is to be written to.
	fp = fopen("graphedgesvsnumconnect.csv", "w");	//create or overwrite the file "graphedgesvsnumconnect.csv"

	if(fp != NULL){
		printf("graphedgesvsnumconnect.csv is created\n");
	}
	else{
		printf("file pointer is pointing to NULL\n");
	}


	int section, nth_section;
	fprintf(fp, "numconnect edges\n");
	for(section=1; section<=max_sections; section++){
		for(nth_section=1; nth_section<=num_outsections[section]; nth_section++){
			fprintf(fp, "%d %d\n", num_outsections[tspans_outsection[section][nth_section]], tspans_edges[section][nth_section]);
		}
	}


}
