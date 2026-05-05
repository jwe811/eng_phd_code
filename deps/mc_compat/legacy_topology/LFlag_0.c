
unsigned short int		/* returns 1 if first section connects to phi on left, 0
				 * otherwise */
LFlag(int (*pointordNum)[3])
{
	/* (*pointordNum)[0] - 1 = number of edges in first section */

	unsigned short int j, k, l;	/*dummy variables */
	unsigned short int validity = 1; /*1 if section connects 0 otherwise */

	if ((*pointordNum)[0] == 3) { /* if there are only two (3-1=2) edges check the orientation */
		for (j = 0; j <= L; j++) { /*loop throuhg each vertex */
					   /*warning: this assumes that the slit is horizontal*/
			if (ordertemplate[0][0][j] == 1) {
				return (1); /*if the edge numbered 1 is reached first the orientation */
					    /* is correct */
			}
			if (ordertemplate[0][0][j] == 2) {
				return (0);	/*if the edge numbered 2 is reached first the orientation */
						/* is not correct */
			}
		}
	}

	/* for more than two edges the following is necessary */

	for (k = 1; k <= ((*pointordNum)[0] / 2) - 1; k++) {
		validity *= noncrossing(0, 1, (*pointordNum)[0] - 1, 2 * k, 2 * k + 1);
		/* the first parameter passed to noncrossing is the side */
		/* the rest are edge numbers in two pairs */
		/* check that the pair 1 , and highest numbered edge */
		/* does not cross and pair 2,3 or 4,5 ,etc... */
		/* if a crossing is encountered noncrossing returns zero and */
		/* validity becomes zero and the *= assignment will not change */
		/* if from being zero */
	}

	for (k = 1; k <= ((*pointordNum)[0] / 2) - 1; k++) {
		for (l = (k + 1); l <= ((*pointordNum)[0] / 2) - 1; l++) {
			validity *= noncrossing(0, 2 * k, 2 * k + 1, 2 * l, 2 * l + 1);
		/* check that the pairs 2,3 and 4,5 and 6,7 ,etc... do not cross each */
		/* other */
		}
	}

	/* now check that the orientation is correct */

	for (j = 0; j <= L; j++) {
		if (ordertemplate[0][0][j] % 2 == 1) { 
			return (validity);
			/*if the edge numbered 1 mod 2  is reached first the orientation */
					    /* is correct */
		}
		if (ordertemplate[0][0][j] != 0 && ordertemplate[0][0][j] % 2 == 0) {
			return (0);
			/*if the edge numbered 0 mod 2  (and not just zero) is reached first the orientation */
			/*is not correct */
		}
	}

}
