double
get_Beta(double fugacity) /* This is the 2span_2span version */
{

int			i; /* used in for loops */
unsigned long int 	insection_num; /* the section number of the in_section */
unsigned long int	nth_outsection; /* the nth_section which connects to in_section */
unsigned long int	outsection_num; /* section number of the nth_section which connects to in_section */
unsigned long int	row_tspan_num;  /* row index for two-span matrix */
unsigned long int	col_tspan_num;  /* column index for two-span matrix */
unsigned long int	right_section;	/* */

/*uses globals tspans_outsection,tspans_edges,tspans_nrr */

double	Gprime;       /* Used as intermediate variable in calculation. */
		      /* From the given fugacity (parameter passed to this function) */
		      /* the value of a given matrix entry is calculated */
double  Beta=0;	/* Beta as defined in notes from Soteros */
double  LmultR;	/* left eigenvector multiplied by right eigenvector */

/* The following code basically does these linear algebra operations */
	/* R_Evector[1] = Matrix * R_Evector[0] */
	/* Beta = L_Evector[1] * R_Evector[0]   */
/************************************************/


	/* multiply left and right eigenvector */
	LmultR = 0;
	for (i=1;i<= max_tspans;i++){
		LmultR += (L_Evector[0][i])*(R_Evector[0][i]);
	}
	
	row_tspan_num=0;

	for (i = 1; i <= max_tspans; i++) {
		R_Evector[1][i] = 0; /* initialize to zero */
	}

	for (insection_num=1;insection_num<=max_sections;insection_num++){
	/* loop through all possible sections */

		for (nth_outsection=1;nth_outsection<=num_outsections[insection_num];nth_outsection++){
		/* for each section loops through the sections it connects to */
		/* this determines the two_span */
		/* further looping is required below to loop through */
		/* the sections that connect to this two-span */
		/* to make a three-span */

			row_tspan_num++; /* increment row number */

			/* the in_section and the section it connects to define a */
			/* two-span */
			/* row_tspan_num and col_tspan_num define the entry in the */
			/* two_span matrix which is calculated below (Gprime) */
			/* actually then information is contained in the first two_span */
			/* defined by in_section and the section it connects to */
			/* the product of this matrix entry and the eigenvector component */
			/* corresponding to the column of the matrix entry makes a */
			/* contribuation to the resultant vector in the matrix-vector */
			/* product */
			/* so  Gprime*(R_Evector[0][col_tspan_num]) gets added to */
			/* R_Evector[1][row_tspan_num] */

			Gprime = tspans_edges[insection_num][nth_outsection]*\
				 pow(fugacity,tspans_edges[insection_num][nth_outsection]-1);
			/* Gprime = # of edges * fugacity ^ (#edges-1) */

			outsection_num = tspans_outsection[insection_num][nth_outsection];
			/* section number of the nth_section which connects to in_section */

			for (right_section=1;right_section<=num_outsections[outsection_num];right_section++){
			/* loop through sections that connect to the two-span's last section */
			/* this defines a second two-span */
				col_tspan_num = tspans_nrr[outsection_num][right_section];
				/*tspans_nrr[outsection][right_section] give the number*/
				/* of the two-span which connects to the first two-span */
				/* and hence the column number in the two-span matrix */
				
				R_Evector[1][row_tspan_num] += Gprime*(R_Evector[0][col_tspan_num]);
			}
		}
	}

	for (i = 1; i <= max_tspans; i++) { /* loop through each component of the eigenvector */
		Beta += L_Evector[0][i]*R_Evector[1][i];
		/* equivalent to Beta = L_Evector[0] * Matrix * R_Evector[0] */
		/* so that Beta corresponds to the definition given in Soteros' notes */
		/* once scaled by fugacity/LmultR */

	}

return(Beta*fugacity/LmultR);

}
