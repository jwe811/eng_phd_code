
unsigned short int		/* returns 1 if second section connects to phi on right, 0
				 * otherwise */
RFlag(int (*pointordNum)[3])
{
	/* (*pointordNum)[1] - 1 =  number of edges in ordertemplate[1][1][j] */

	unsigned short int j, k, l;	/*dummy variables*/
	unsigned short int validity = 1; /*1 if section connects 0 otherwise */

	if ((*pointordNum)[1] == 3) { /* if there are only two (3-1) edges check the orientation */
		for (j = 0; j <= L; j++) {/*loop throuhg each vertex */
					   /*warning: this assumes that the slit is horizontal*/
			if (ordertemplate[1][0][j] == 1) { 
				return (1);/*if the edge numbered 1 is reached first the orientation */
					    /* is correct */
			}
			if (ordertemplate[1][0][j] == 2) {
				return (0);	/*if the edge numbered 2 is reached first the orientation */
						/* is not correct */
			}
		}
	}

	/* for more than two edges the following is necessary */

	for (k = 1; k <= (*pointordNum)[1] / 2; k++) {
		for (l = (k + 1); l <= (*pointordNum)[1] / 2; l++) {

			validity *= noncrossing(1, 2 * k - 1, 2 * k, 2 * l - 1, 2 * l);
		/* the first parameter passed to noncrossing is the side */
		/* the rest are edge numbers in two pairs */
		/* check that the pairs 1,2 and 3,4 and 5,6 ,etc... do not cross each */
		/* other */
		}
	}

	/* now check that the orientation is correct */

	for (j = 0; j <= L; j++) {
		if (ordertemplate[1][0][j] % 2 == 1) {
			return (validity);
			/*if the edge numbered 1 mod 2  is reached first the orientation */
					    /* is correct */
		}
		if (ordertemplate[1][0][j] != 0 && ordertemplate[1][0][j] % 2 == (*pointordNum)[1]) {
			return (0);
			/*if the edge numbered 0 mod 2  (and not just zero) is reached first the orientation */
			/*is not correct */
		}
	}
}
