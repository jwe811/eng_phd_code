/***************************************************************************/
//#define MAXIT_pw_method 1000

///////////////////////USE THIS IF TRANSFER MATRIX IS PERIODIC////////////////////////
//PERFORMS METHOD ON sI+G. RETURNS max eigenvalue of G (which is equal to [max eigenvalue of (sI+G)] - s), then minus 1 cause that's what has been always done.
//sI+G has same eigenvector as G


#define MAXIT_pw_method 10000
#define theshift 100	//shift (s) used to make matrix aperiodic ( (L+1)*(M+1) is odd)
#define accuracy 0.0000000001

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
	double	pwfugacity; /*fugacity raised to an appropriate power */

	///////////////////////////////////////////////////////////////////////////////
	double lambda_min;
	double lambda_max;
	double min_ratio;
	double max_ratio;
	double cur_ratio;

	double themax_R;
	double themax_L;
	double difference;
	///////////////////////////////////////////////////////////////////////////////


	/*uses globals tspans_outsection,tspans_edges,tspans_nrr*/

	printf("\nUsing Modified Power Method with fugacity = %f\n\n",fugacity);

	for(i=1;i<=max_tspans;i++){
		L_Evector[0][i]=1;
		R_Evector[0][i]=1;
	}
	//loop through components of eigenvectors and initialize to 1. this is the initial guess vector


	for (k=1;k<=MAXIT_pw_method;k++){/*iterate the matrix multiplications a maximum of MAXIT times */

		/* L_Evector[1] = L_Evector[0] * Matrix */
		/* R_Evector[1] = Matrix * R_Evector[0] */

		row_tspan_num=0;
		//resultant vectors
		for (i=1;i<=max_tspans;i++){
			L_Evector[1][i] = 0;
			R_Evector[1][i] = 0;
		}

		for (insection_num=1;insection_num<=max_sections;insection_num++){
			/* loop through all possible sections */
			for (nth_outsection=1;nth_outsection<=num_outsections[insection_num];nth_outsection++){
			/* for each section loops through the sections it connects to */
				row_tspan_num++;	/*each step through above insection_num and nth_outsection loops*/
								/* defines a new row in the two-span matrix*/
				pwfugacity = pow(fugacity,tspans_edges[insection_num][nth_outsection]);
				/*raises fugacity to number of edges */

				/*for each two span that connects to this two span do*/
				outsection_num = tspans_outsection[insection_num][nth_outsection];
				/* section number of the nth_section which connects to in_section */
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

		if( ((L+1)*(M+1)) % 2 != 0 ){
//			printf("L*M is odd. Doing different calculation\n");
			//Add I to calculation.
			//q_{k+1} = (sI+G)*q_k = s*q_k + G*q_k
			for (i=1;i<=max_tspans;i++){
				L_Evector[1][i] += theshift*L_Evector[0][i];
				R_Evector[1][i] += theshift*R_Evector[0][i];
			}
		}




	/////////////////////////////////////////////////
		//calculate lambda_min and lambda_max
		themax_L = L_Evector[1][1];
		themax_R = R_Evector[1][1];
		min_ratio = R_Evector[1][1] / R_Evector[0][1];
		max_ratio = R_Evector[1][1] / R_Evector[0][1];
		for (i=2;i<=max_tspans;i++){
			cur_ratio = R_Evector[1][i] / R_Evector[0][i];
			if(cur_ratio > max_ratio){
				max_ratio = cur_ratio;
			}
			if(cur_ratio < min_ratio){
				min_ratio = cur_ratio;
			}
			if(R_Evector[1][i] > themax_R){
				themax_R = R_Evector[1][i];
			}
			if(L_Evector[1][i] > themax_L){
				themax_L = L_Evector[1][i];
			}
		}

		lambda_min = min_ratio;
		lambda_max = max_ratio;
//		printf("Iteration %d. lambda_min=%f, lambda_max=%f\n", k, lambda_min, lambda_max);
		///////////////////////////////////////////////////////

		for (i=1;i<=max_tspans;i++) {
			R_Evector[1][i] /= themax_R;
			L_Evector[1][i] /= themax_L;
		}
		/* "normalize" the eigenvectors */


	/**/	/* L_Evector[0] = L_Evector[1] * Matrix */
		/* R_Evector[0] = Matrix * R_Evector[1] */

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

		if( ((L+1)*(M+1)) % 2 != 0 ){
			//Add I to calculation.
			//q_{k+1} = (sI+G)*q_k = s*q_k + G*q_k
			for (i=1;i<=max_tspans;i++){
				L_Evector[0][i] += theshift*L_Evector[1][i];
				R_Evector[0][i] += theshift*R_Evector[1][i];
			}
		}


	/////////////////////////////////////////////////
		//calculate lambda_min and lambda_max
		themax_L = L_Evector[0][1];
		themax_R = R_Evector[0][1];
		min_ratio = R_Evector[0][1] / R_Evector[1][1];
		max_ratio = R_Evector[0][1] / R_Evector[1][1];
		for (i=2;i<=max_tspans;i++){
			cur_ratio = R_Evector[0][i] / R_Evector[1][i];
			if(cur_ratio > max_ratio){
				max_ratio = cur_ratio;
			}
			if(cur_ratio < min_ratio){
				min_ratio = cur_ratio;
			}
			if(R_Evector[0][i] > themax_R){
				themax_R = R_Evector[0][i];
			}
			if(L_Evector[0][i] > themax_L){
				themax_L = L_Evector[0][i];
			}
		}

		lambda_min = min_ratio;
		lambda_max = max_ratio;
		printf("Iteration %d. lambda_min=%f, lambda_max=%f.  ", k, lambda_min, lambda_max);
	///////////////////////////////////////////////////////

		for (i=1;i<=max_tspans;i++) {
			R_Evector[0][i] /= themax_R;
			L_Evector[0][i] /= themax_L;
		}
		/* "normalize" the eigenvectors */

		//check for convergence
		difference = lambda_max - lambda_min;
		printf("Difference is %f = %e\n", difference, difference);

		if(difference < accuracy){
			if( ((L+1)*(M+1)) % 2 != 0 ){
				return(lambda_max - theshift -1);	//minus s since we added sI (since there was periodicity)
										//minus 1 since that's how it's always been done....
			}
			else{
				return(lambda_max-1);
			}
		}

	} /*k loop*/

	fprintf(stderr, "Maximum number of iterations exceeded in power method\n");
	/*error message if trouble is encountered in this function*/
	exit(1);

} /*end function max_eigenval*/

#undef MAXIT_pw_method
/***************************************************************************/
