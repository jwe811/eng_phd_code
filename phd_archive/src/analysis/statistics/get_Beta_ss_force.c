double
get_Beta(double fval) /* This is the section_section version */
{
double fugacity=connectivity_inv;

/*uses globals tspans_outsection,tspans_edges */

int	i;	/* used in for loops */
unsigned long int insection_num;  /* the section number of the in_section */
unsigned long int nth_outsection; /* the nth_section which connects to in_section */
unsigned long int outsection_num; /* section number of the nth_section which connects to in_section */

double	Gprime;	/* Used as intermediate variable in calculation. */
		/* From the given fugacity (parameter passed to this function) */
		/* the value of a given matrix entry is calculated */ 
double  Beta=0;	/* Beta as defined in notes from Soteros */
double  LmultR;	/* left eigenvector multiplied by right eigenvector */

/* The following code basically does these linear algebra operations */
	/* R_Evector[1] = Matrix * R_Evector[0] */
	/* Beta = L_Evector[1] * R_Evector[0]   */
/************************************************/

	LmultR = 0;	/* initialize to zero */
	for (i=1;i<= max_sections;i++){ /* take sum over corresponding components of the eigenvectors */
		LmultR += (L_Evector[0][i])*(R_Evector[0][i]);	
	}
	
	for (i = 1; i <= max_sections; i++) {
		R_Evector[1][i] = 0; /* initiallize to zero */
	}

	for (insection_num=1;insection_num<=max_sections;insection_num++){
	/* loop through all possible sections */
		#if defined(CS)
		outsection_num = 0; /*reset so distinct two-spans get recognized */
		#endif

		for (nth_outsection=1;nth_outsection<=num_outsections[insection_num];nth_outsection++){
		/* for each section loops through the sections it connects to */

			/* the in_section and the section it connects to defines an */
			/* entry in the matrix which is calculated below (Gprime) */
			/* the product of this matrix entry and the eigenvector component */
			/* corresponding to the column of the matrix entry makes a */
			/* contribuation to the resultant vector in the matrix-vector */
			/* product */
			/* so  Gprime*(R_Evector[0][outsection_num]) gets added to */
			/* R_Evector[1][insection_num] */
	#if defined (CS)
	if (tspans_outsection[insection_num][nth_outsection] == outsection_num) { /*not a new two-span*/
		Gprime = 0;
	}
	else
	#endif
			Gprime = tspans_edges[insection_num][nth_outsection]*\
				 pow(fugacity,(tspans_edges[insection_num][nth_outsection])-1)*exp(fval);
			/* Gprime = # of edges * fugacity ^ (#edges-1) */

			outsection_num = tspans_outsection[insection_num][nth_outsection];
			/* section number of the nth_section which connects to in_section */

				R_Evector[1][insection_num] += Gprime*(R_Evector[0][outsection_num]);
				
		}
	}

	for (i = 1; i <= max_sections; i++) {	/* loop through each component of the eigenvector */
		Beta += L_Evector[0][i]*R_Evector[1][i];
		/* equivalent to Beta = L_Evector[0] * Matrix * R_Evector[0] */
		/* so that Beta corresponds to the definition given in Soteros' notes */
		/* once scaled by fugacity/LmultR */
	}

return(Beta*fugacity/LmultR);

}
