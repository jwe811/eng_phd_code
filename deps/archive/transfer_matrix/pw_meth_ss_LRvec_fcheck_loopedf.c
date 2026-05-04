/***************************************************************************/
#define MAXIT_pw_method 500

double
max_eval_LRvec(double fugacity)
{
/*uses globals tspans_outsection,tspans_edges,tspans_nrr,num_outsections*/


int	i,j,k; /*dummy variables */
unsigned long int insection_num;  /*the section number of the in_section*/
unsigned long int nth_outsection; /*the nth_section which connects to in_section */
unsigned long int outsection_num; /*section number of the nth section which connects to in_section*/
double	L_lambda_0; /*left eigenvalue zero  */
double  L_lambda_1; /*left eigenvalue one   */
double  R_lambda_0; /*right eigenvalue zero */
double  R_lambda_1; /*right eigenvalue one  */
double	pwfugacity; /*fugacity raised to an appropriate power */

//printf("\nUsing Power Method with fugacity = %f\n\n",fugacity);
//printf("Iteration     Left eigenvalue     Difference     Right eigenvalue     Difference\n\n");

for (i=1;i<=max_sections;i++) { L_Evector[0][i]=1; R_Evector[0][i]=1;}
/*loop through components of eigenvectors and initialize to 1 this is the initial guess in power method*/

for (k=1;k<=MAXIT_pw_method;k++){ /*iterate the matrix multiplications a maximum of MAXIT times */
/**/	/* L_Evector[1] = L_Evector[0] * Matrix */
	/* R_Evector[1] = Matrix * R_Evector[0] */
	
	L_lambda_0 = 0.1; /*set to just a little bigger than zero to avoid dividing by zero */
	R_lambda_0 = 0.1; /*set to just a little bigger than zero to avoid dividing by zero */
	for (i = 1; i <= max_sections; i++) { /*initialize resultant vectors to zero */
		L_Evector[1][i] = 0;
		R_Evector[1][i] = 0;
	}

	for (insection_num=1;insection_num<=max_sections;insection_num++){
		/* loop through all possible sections */
		#if defined(CS)
		outsection_num = 0; /*reset so distinct two-spans get recognized */
		#endif

		for (nth_outsection=1;nth_outsection<=num_outsections[insection_num];nth_outsection++){
		/* for each section loops through the sections it connects to */

			#if defined (CS)
			if (tspans_outsection[insection_num][nth_outsection] == outsection_num) {
				pwfugacity = 0;
			}
			else
			#endif

//			printf("insection=%d, outsection=%d, edges=%d\n", insection_num, nth_outsection, tspans_edges[insection_num][nth_outsection]);

			pwfugacity = pow(fugacity, tspans_edges[insection_num][nth_outsection]);
				/*raises fugacity to number of edges */

			outsection_num = tspans_outsection[insection_num][nth_outsection];
			/* section number of the nth_section which connects to in_section */

			L_Evector[1][outsection_num] += (L_Evector[0][insection_num])*pwfugacity*exp(fval);
			/* contribution to resultant vector in L_Evector[0] * Matrix */
			R_Evector[1][insection_num] += (R_Evector[0][outsection_num])*pwfugacity*exp(fval);
			/* contribution to resultant vector in Matrix * R_Evector[0] */

		}
	}

	for (i=1;i<=max_sections;i++){/*check for largest component in resultant vectors*/
		if (L_Evector[1][i] > L_lambda_0) { L_lambda_0 =  L_Evector[1][i]; }
		if (R_Evector[1][i] > R_lambda_0) { R_lambda_0 =  R_Evector[1][i]; }
		/*and set eigenvalue (ie. eventually converges to eigenvalue) to largest component*/
	}

	for (i=1;i<=max_sections;i++) { L_Evector[1][i] /= L_lambda_0; R_Evector[1][i] /= R_lambda_0; }
	/* "normalize" the eigenvectors */

/**/	/* L_Evector[0] = L_Evector[1] * Matrix */
	/* R_Evector[0] = Matrix * R_Evector[1] */

	L_lambda_1 = 0.1; /*set to just a little bigger than zero to avoid dividing by zero */
	R_lambda_1 = 0.1; /*set to just a little bigger than zero to avoid dividing by zero */

	for (i=1;i<=max_sections;i++){ L_Evector[0][i] = 0; R_Evector[0][i] = 0; }
	/*initialize resultant vectors to zero */

	for (insection_num=1;insection_num<=max_sections;insection_num++){
	/* loop through all possible sections */
		#if defined(CS)
		outsection_num = 0; /*reset so distinct two-spans get recognized */
		#endif

		for (nth_outsection=1;nth_outsection<=num_outsections[insection_num];nth_outsection++){
		/* for each section loops through the sections it connects to */

			#if defined (CS)
			if (tspans_outsection[insection_num][nth_outsection] == outsection_num) {
				pwfugacity = 0;
			}
			else
			#endif

			pwfugacity = pow(fugacity, tspans_edges[insection_num][nth_outsection]);
			/* matrix entry corresponding to row = insection , column  = outsection */
			outsection_num = tspans_outsection[insection_num][nth_outsection];
			/* section number of the nth_section which connects to in_section */

			L_Evector[0][outsection_num] += pwfugacity*(L_Evector[1][insection_num])*exp(fval);
			R_Evector[0][insection_num] += pwfugacity*(R_Evector[1][outsection_num])*exp(fval);
		}
	}


	for (i=1;i<=max_sections;i++){/*check for largest component in resultant vectors*/
		if (L_Evector[0][i] > L_lambda_1) { L_lambda_1 =  L_Evector[0][i]; }
		if (R_Evector[0][i] > R_lambda_1) { R_lambda_1 =  R_Evector[0][i]; }
		/*and set eigenvalue (ie. eventually converges to eigenvalue) to largest component*/
	}
	for (i=1;i<=max_sections;i++) { L_Evector[0][i] /= L_lambda_1; R_Evector[0][i] /= R_lambda_1; }
	/* "normalize" the eigenvectors */


/**/ /*Check for convergence */
	if ((-0.0000001 < (L_lambda_1-L_lambda_0)) && ((L_lambda_1-L_lambda_0) < 0.0000001)\
	&&(-0.0000001 < (R_lambda_1-R_lambda_0)) && ((R_lambda_1-R_lambda_0) < 0.0000001)){
//		printf("\n%6i (DONE)",2*k);
//		printf("%13f %17f  ",L_lambda_1,L_lambda_1-L_lambda_0);
//		printf("%16f %17f\n",R_lambda_1,R_lambda_1-R_lambda_0);
		return(R_lambda_1 - 1); /*R_Evector[0] = right eigenvector corresponding to R_lambda_1*/
	}
	/*pre-convergence output to user, to give idea of how the eigenvalue converges*/
	
//	printf("%6i",2*k);
//	printf("%20f %17f  ",L_lambda_1,L_lambda_1-L_lambda_0);
//	printf("%16f %17f\n",R_lambda_1,R_lambda_1-R_lambda_0);
		
} /*k loop*/

fprintf(stderr, "Maximum number of iterations exceeded in power method\n");
/*error message if trouble is encountered in this function*/
exit(1);

} /*end function max_eigenval*/

#undef MAXIT_pw_method
/***************************************************************************/
