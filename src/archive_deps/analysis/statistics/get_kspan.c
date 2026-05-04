void /*actually need to return two-span #*/
     /*or first section # and last section # */
     /*and total number of edges*/

get_tspan(unsigned short int k)
{

	unsigned short int n = 0;
	unsigned short int valid;
	unsigned long  int tspan_1_num;
	unsigned long  int tspan_k_num;
	unsigned long  int sec_1_num;
	unsigned long  int sec_k_num;


	while (n < k) {

		n = n + 1;
		printf("Enter %i row(s) of %i numbers ",vM,vL);
		printf("for section #%i:\n", n);

		until(valid==1) { /*get first section */
			scan_template(0);	/* enter section #n */
			valid=entsecval(0); /*valid and connecting*/
		}

		until(valid==1) {
			/*get_numedges;*/
			valid=1;
		}

		until(valid==1) { /* get second section */
			scan_template(1);	/* enter section #n */
			valid=entsecval(1)*entseccon(0); /*valid and connecting*/
		}


	}
}

unsigned short int
entsecval(unsigned short int side) /* entered section validity test */
{
	unsigned short int section_num = num_section(side);

	if (1 <= section_num && section_num <= max_sections) {	/* ensure in range */
		if (num_outsections[section_num] != 0) {	/* ensure it connects to
								 * other sections */
		return(1); /*entered section is valid*/
		} else {
			printf("This section does not connect to anything.\n");
		}
	} else {
		printf("The entered section is not valid!\n");
	}
	return(0); /*The entered section is not valid.*/
}


unsigned short int		/* entered section connects test */
entseccon(unsigned short int side)
{

	unsigned long int section_1 = numsection(side);
	unsigned long int section_2 = numsection(side ^ 1);
	unsigned short int result = 0;

	for (nth_tspan = 1; nth_tspan <= num_outsections[section_1]; nth_tspan++) {
		if (section_2 == tspans_outsection[section_1][nth_tspan]) {
			printf("\n%6u:%6u : edges=%3u : two span# =%3u\n", \
			       section_num, tspans_outsection[section_num][nth_tspan], \
			       tspans_edges[section_num][nth_tspan], \
			       tspans_nrr[section_num][nth_tspan]);
			result = 1;
		}
	}
	if (result == 0) {
		printf("This section does not connect to the previous section.\n");
	}
	return (result);
}




	




