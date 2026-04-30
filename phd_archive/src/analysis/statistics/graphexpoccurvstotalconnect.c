void graphexpoccurvstotalconnect(double Beta,double LmultR,double kappa) {
/*
This function will create a .csv file named graphexpoccurvstotalconnect.csv.  First column will be the total number of sections that
can come before the 2-span plus the total number of sections that can come after the 2-span.  The second column will be the expected
number of occurrences of that 2-span.
*/

	FILE *fp;	//file pointer. will point to the file that is to be written to.
	fp = fopen("graphexpoccurvstotalconnect.csv", "w");	//create or overwrite the file "graphexpoccurvstotalconnect.csv"

	if(fp != NULL){
		printf("graphexpoccurvstotalconnect.csv is created\n");
	}
	else{
		printf("file pointer is pointing to NULL\n");
	}


	int section, nth_section, prevsection, nth_prevsection;
	int numofinsections, numconnect, totalconnect;

	fprintf(fp, "totalconnect expoccur\n");
	for(section=1; section<=max_sections; section++){
		for(nth_section=1; nth_section<=num_outsections[section]; nth_section++){
			//for each 2-span
			numconnect=num_outsections[tspans_outsection[section][nth_section]];

			//look at each section and see if it has an outsection which matches "section"
			numofinsections=0;
			totalconnect=0;
			for(prevsection=1; prevsection<=max_sections; prevsection++){
				for(nth_prevsection=1; nth_prevsection<=num_outsections[prevsection]; nth_prevsection++){
					if(tspans_outsection[prevsection][nth_prevsection]==section){
						numofinsections++;
					}
				}
			}
			totalconnect=numconnect+numofinsections;
			fprintf(fp, "%d %f\n", totalconnect, expoccur(Beta, LmultR, kappa, tspans_nrr[section][nth_section]));
		}
	}


}
