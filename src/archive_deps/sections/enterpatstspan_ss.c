//USES section-section MATRIX TO FIND EXPECTED OCCURENCES FOR A 2-SPAN

double enterpats(double Beta, double LmultR, double kappa)
{				/* begin */

#if !defined(connectivity_inv)	/* if connectivity_inv is not defined declare it here */

	double			connectivity_inv=kappa;

#endif

	unsigned long int entering_patterns = 1;
	unsigned long int total_span;
	unsigned short int total_edges = 0;	/* counts the number of edges
						 * as the k-span is entered
						 * by */
	/* the user one section at a time */
	unsigned long int nth_section = 1;
	unsigned short int section_validity;	/* 1 if the an entered
						 * section is valid 0
						 * otherwise */
	unsigned long int section_num;
	int             side = 0;	/* WARNING: same variable name used
					 * in main */
	/* this hopefully does not matter */
	unsigned long int section_2_num;
	unsigned long int nth_tspan;
	unsigned long int connecting_tspans;
	unsigned short int hinge_validity;	/* 1 if subsequent entered
						 * section can be connected
						 * to  */
	/* the previously entered section via a valid hinge     */
	/* 0 otherwise */
	unsigned long int tspan_num = 1;
	unsigned long int first_tspan = 1;
	double			Enum_Enum;	/* Warning: enum is a keyword in C */


	while (entering_patterns != 0) {
		total_span=2;
		newline;
		total_edges = 0;
		nth_section = 1;
		printf("The number of sections to be entered is: %i\n", total_span);
		printf("Please enter %i row(s) of %i numbers for each section:\n", vM, vL);

		/*** ENTER FIRST SECTION ***/

		section_validity = NOTVALID;	/* NOTVALID until VALID */
		while (section_validity == NOTVALID) {	/* get a valid first
							 * section */

			printf("\n\nEnter section #%i:\n\n", nth_section);
			scan_template(0);	/* enter section #1 */

			printf("\nThis is what you entered (section#%i=%i):\n", nth_section, num_section(0));
			echo_template(0);	/* echo to user */

			/* check validity */
			section_num = num_section(0);
			if (1 <= section_num && section_num <= max_sections) {	/* ensure in range */
				if (num_outsections[section_num] != 0) {	/* ensure it connects to
										 * other sections */
					section_validity = VALID;
				} else {
					printf("This section does not connect to anything.\n");
				}
			} else {
				printf("The entered section is not valid!\n");
			}
		}

		nth_section++;	/* increment section number, since the user
				 * has entered a valid section */

		/*** FIRST SECTION ENTERED ***/


		side = 1;
		while (nth_section <= total_span) {	/* Accept user input
							 * until hinge edges
							 * entered = 0 */
			section_validity = NOTVALID;
			connecting_tspans = 0;
			while (section_validity == NOTVALID) {

				printf("\n\nEnter section #%i:\n\n", nth_section);
				scan_template(side);	/* enter section
							 * #nth_section */

				printf("\nThis is what you entered (section#%i=%i):\n", nth_section, section_2_num = num_section(side));
				echo_template(side);	/* echo to user */
				/*
				 * check validity and check compatibility
				 * with previous section
				 */
				if (1 <= section_2_num && section_2_num <= max_sections) {	/* ensure in range */
					newline;
					for (nth_tspan = 1; nth_tspan <= num_outsections[section_num]; nth_tspan++) {
						/*
						 * if compatible display
						 * available hinges/two-spans
						 */
						if (section_2_num == tspans_outsection[section_num][nth_tspan]) {
							printf("\n%6u:%6u : edges=%3u : two span# =%3u\n", \
							       section_num, tspans_outsection[section_num][nth_tspan], \
							       tspans_edges[section_num][nth_tspan], \
							       tspans_nrr[section_num][nth_tspan]);
							connecting_tspans++;
						}
					}
					if (connecting_tspans != 0) {
						section_validity = VALID;
					} else {
						printf("\n\nThis section does not connect to the previous section.\n");
					}
				} else {
					printf("\nThe entered section is not valid!\n");
				}	/* else section_validity == NOTVALID */

				if (section_validity == VALID) {
					hinge_validity = NOTVALID;
					while (hinge_validity == NOTVALID) {
						/*
						 * get # edges in hinge via
						 * tspan or otherwise
						 */
						printf("\nEnter a two span # from the above list to complete %i-span (type 0 to quit):", nth_section);
						scanf("%u", &(tspan_num));
						if (tspan_num == 0) {
							nth_section = total_span;	/* all sections entered */
							break;
						}
						if (nth_section == 2) {
							first_tspan = tspan_num;
						}
						/* check validity */
						for (nth_tspan = 1; nth_tspan <= num_outsections[section_num]; nth_tspan++) {
							if (tspan_num == tspans_nrr[section_num][nth_tspan] && \
							    section_2_num == tspans_outsection[section_num][nth_tspan]) {
								hinge_validity = VALID;
								total_edges += tspans_edges[section_num][nth_tspan];
								break;
							}
						}
						if (hinge_validity == NOTVALID) {
							printf("That two span is not in the list\n");
						}
					}
				}
			}	/* while section_validity == NOTVALID */
			nth_section++;	/* increment section number, since
					 * the user has entered a valid
					 * section */
			side ^= side;	/* switch "sides" */
			//section_num = section_2_num;	/* move on to next section */
		}

	//tspan is section_num and section_2_num=tspans_outsection[section_num][nth_tspan], with edges = total_edges

	printf("entered tspan#%d=%d, section %d to section %d=%d, with edges=%d=%d\n", tspans_nrr[section_num][nth_tspan], tspan_num, section_num, section_2_num, tspans_outsection[section_num][nth_tspan], total_edges, tspans_edges[section_num][nth_tspan]);

	Enum_Enum = (L_Evector[0][section_num]) * (R_Evector[0][section_2_num]) / LmultR;
	printf("L_Evector[%i]*R_Evector[%i] = %f = %e\n\n", section_num, section_2_num, Enum_Enum, Enum_Enum);
	printf("The expected number of times this %i-span pattern occurs\n", --nth_section);
	printf("in an n-step polygon is (as n->infinity):\n\n");
	printf("(%f)n = (%e)n\n", Enum_Enum * pow(connectivity_inv, total_edges) * exp(fval) / Beta, Enum_Enum * pow(connectivity_inv, total_edges) * exp(fval) / Beta);
	printf("\nType 1 to enter another arbitrary span 0 to quit:");
	scanf("%u",&entering_patterns);

	}
	return(Enum_Enum);
}				/* end */

