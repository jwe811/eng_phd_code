/***************************************************************************/
#define MAXIT_pw_method 1000

double
max_eval_LRvec(double fugacity)
{
/*uses globals tspans_outsection,tspans_edges,tspans_nrr,num_outsections*/


int	i,k;				/*used in for loops*/
unsigned long int 	insection_num;  /* the section number of the in_section */
unsigned long int	nth_outsection; /* the nth_section which connects to in_section */
unsigned long int	outsection_num; /* section number of the nth_section which connects to in_section */
unsigned long int	row_tspan_num;  /* row index for two-span matrix */
unsigned long int	col_tspan_num;  /* column index for two-span matrix */
unsigned long int 	right_tspan;
double	L_lambda_0; /*left eigenvalue zero  */
double  L_lambda_1; /*left eigenvalue one   */
double  R_lambda_0; /*right eigenvalue zero */
double  R_lambda_1; /*right eigenvalue one  */
double	pwfugacity; /*fugacity raised to an appropriate power */


/*uses globals tspans_outsection,tspans_edges,tspans_nrr*/

printf("\nUsing Power Method with fugacity = %f\n\n",fugacity);
printf("Iteration     Left eigenvalue     Difference     Right eigenvalue     Difference\n\n");

for (i=1;i<=max_tspans;i++) { L_Evector[0][i]=1; R_Evector[0][i]=1; }
/*loop through components of eigenvectors and initialize to 1 this is the initial guess in power method*/

for (k=1;k<=MAXIT_pw_method;k++){/*iterate the matrix multiplications a maximum of MAXIT times */

/**/	/* L_Evector[1] = L_Evector[0] * Matrix */
	/* R_Evector[1] = Matrix * R_Evector[0] */
	
	L_lambda_0 = 0.0;
	R_lambda_0 = 0.0;
	row_tspan_num=0;
	for (i=1;i<=max_tspans;i++){ L_Evector[1][i] = 0; R_Evector[1][i] = 0; }

	for (insection_num=1;insection_num<=max_sections;insection_num++){
		/* loop through all possible sections */
		if (num_outsections[insection_num] > 0) {
			// printf("Section %lu has %lu outsections\n", insection_num, num_outsections[insection_num]);
		}
		for (nth_outsection=1;nth_outsection<=num_outsections[insection_num];nth_outsection++){
		/* for each section loops through the sections it connects to */
			row_tspan_num++; /*each step through above insection_num and nth_outsection loops*/
					  /* defines a new row in the two-span matrix*/

			pwfugacity = pow(fugacity,tspans_edges[insection_num][nth_outsection]);
			/*raises fugacity to number of edges */

       /*for each two span that connects to this two span do*/
			outsection_num = tspans_outsection[insection_num][nth_outsection];
			/* section number of the nth_section which connects to in_section */
/**HERE**/
			for (right_tspan=1;right_tspan<=num_outsections[outsection_num];right_tspan++){
			/*loop through sections which connect on the right*/
				col_tspan_num = tspans_nrr[outsection_num][right_tspan];
				/*column number = tspan # of two span made of outsection_num and right_tspan*/
				/*the following picks out the contribution to the resultant vectors in the */
				/*matrix multiplication*/
				L_Evector[1][col_tspan_num] += pwfugacity*(L_Evector[0][row_tspan_num])*exp(fval); /*left  eigenvector*/
				R_Evector[1][row_tspan_num] += pwfugacity*(R_Evector[0][col_tspan_num])*exp(fval); /*right eigenvector*/
			}
		}
	}

	for (i=1;i<=max_tspans;i++){/*check for largest component in resultant vectors*/
		if (L_Evector[1][i] > L_lambda_0) { L_lambda_0 =  L_Evector[1][i]; }
		if (R_Evector[1][i] > R_lambda_0) { R_lambda_0 =  R_Evector[1][i]; }
	/*and set eigenvalue (ie. eventually converges to eigenvalue) to largest component*/
	}
	if (L_lambda_0 == 0.0) L_lambda_0 = 1.0;
	if (R_lambda_0 == 0.0) R_lambda_0 = 1.0;
	for (i=1;i<=max_tspans;i++) { L_Evector[1][i] /= L_lambda_0; R_Evector[1][i] /= R_lambda_0; }
	/* "normalize" the eigenvectors */


/**/	/* L_Evector[0] = L_Evector[1] * Matrix */
	/* R_Evector[0] = Matrix * R_Evector[1] */

	L_lambda_1 = 0.0;
	R_lambda_1 = 0.0;
	row_tspan_num=0;
	for (i=1;i<=max_tspans;i++){ L_Evector[0][i] = 0; R_Evector[0][i] = 0; }

	for (insection_num=1;insection_num<=max_sections;insection_num++){
		/* loop through all possible sections */
		for (nth_outsection=1;nth_outsection<=num_outsections[insection_num];nth_outsection++){
		/* for each section loops through the sections it connects to */
			row_tspan_num++;/*each step through above insection_num and nth_outsection loops*/
					  /* defines a new row in the two-span matrix*/
			pwfugacity = pow(fugacity,tspans_edges[insection_num][nth_outsection]);
			/*raises fugacity to number of edges */

/*for each two span that connects to this two span do*/
			outsection_num = tspans_outsection[insection_num][nth_outsection];

			for (right_tspan=1;right_tspan<=num_outsections[outsection_num];right_tspan++){
			/*loop through sections which connect on the right*/
				col_tspan_num = tspans_nrr[outsection_num][right_tspan];
				/*column number = tspan # of two span made of outsection_num and right_tspan*/
				/*the following picks out the contribution to the resultant vectors in the */
				/*matrix multiplication*/
				L_Evector[0][col_tspan_num] += pwfugacity*(L_Evector[1][row_tspan_num])*exp(fval); /*left  eigenvector*/
				R_Evector[0][row_tspan_num] += pwfugacity*(R_Evector[1][col_tspan_num])*exp(fval); /*right eigenvector*/
			}
		}
	}


	for (i=1;i<=max_tspans;i++){/*check for largest component in resultant vectors*/
		if (L_Evector[0][i] > L_lambda_1) { L_lambda_1 =  L_Evector[0][i]; }
		if (R_Evector[0][i] > R_lambda_1) { R_lambda_1 =  R_Evector[0][i]; }
	/*and set eigenvalue (ie. eventually converges to eigenvalue) to largest component*/
	}
	if (L_lambda_1 == 0.0) L_lambda_1 = 1.0;
	if (R_lambda_1 == 0.0) R_lambda_1 = 1.0;
	for (i=1;i<=max_tspans;i++) { L_Evector[0][i] /= L_lambda_1; R_Evector[0][i] /= R_lambda_1; }
	/* "normalize" the eigenvectors */

/**/
/*check for convergence*/

	if ((-0.0000001 < (L_lambda_1-L_lambda_0)) && ((L_lambda_1-L_lambda_0) < 0.0000001)\
	   &&(-0.0000001 < (R_lambda_1-R_lambda_0)) && ((R_lambda_1-R_lambda_0) < 0.0000001)){
	printf("\n%6i (DONE)",2*k);
	printf("%13f %17f  ",L_lambda_1,L_lambda_1-L_lambda_0);
	printf("%16f %17f\n",R_lambda_1,R_lambda_1-R_lambda_0);
		return(R_lambda_1 - 1); /*R_Evector[0] = right eigenvector corresponding to R_lambda_1*/
	}
	
	/*pre-convergence output to user, to give idea of how the eigenvalue converges*/
	printf("%6i",2*k);
	printf("%20f %17f  ",L_lambda_1,L_lambda_1-L_lambda_0);
	printf("%16f %17f\n",R_lambda_1,R_lambda_1-R_lambda_0);
		
} /*k loop*/

fprintf(stderr, "Maximum number of iterations exceeded in power method\n");
/*error message if trouble is encountered in this function*/
exit(1);
		return(R_lambda_1 - 1); /*R_Evector[0] = right eigenvector corresponding to R_lambda_1*/

} /*end function max_eigenval*/

#undef MAXIT_pw_method
/***************************************************************************/
