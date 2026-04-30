unsigned short int		/* returns 1 if testhinge connects, 0 otherwise */
LFlag(int (*pointordNum)[3])
{
	if ((*pointordNum)[0] - 1 < vM*vL) {	/* if there are less edges in first section */
						/* than there are number of vertices */
						/* ie. if first section is not full */

/* Strictly speaking this should actually read (*pointordNum) - 1 < vM*vL */
/* But since (*pointordNum) is odd and vM*vL is even it is not necessary */

		return (1);			/* the section connects on the left */
						/* note: first section = left section */
	} else { /*all edges are occupied in this section */
		unsigned short int i, j,b; /*dummy variables*/

		for (i = 0; i <= M; i++) {
			for (j = 0; j < L; j++) { /*the i and j loops serve to loop through*/
						/* all posible pairs of horizontally neighboring edges */
						/* that is i,j and i,j+1  hence the j loop ends at L-1 */

				int diff = (ordertemplate[0][i][j]) - (b=ordertemplate[0][i][j + 1]);
				if (diff == 1) {
					if (b != 0 && (b & 1) == 0) return(1);
				} else if (diff == -1) {
					if (b != 1 && (b & 1) == 1) return(1);
				} else if (diff == (vM*vL-1) || diff == (1-vM*vL)) {
					return(1);
				}
			}
		}

		for (i = 0; i < M; i++) {
			for (j = 0; j <= L; j++) {/*the i and j loops serve to loop through*/
						/* all posible pairs of vertically neighboring edges */
						/* that is i,j and i+1,j  hence the i loop ends at M-1 */

				int diff = (ordertemplate[0][i][j]) - (b=ordertemplate[0][i+1][j]);
				if (diff == 1) {
					if (b != 0 && (b & 1) == 0) return(1);
				} else if (diff == -1) {
					if (b != 1 && (b & 1) == 1) return(1);
				} else if (diff == (vM*vL-1) || diff == (1-vM*vL)) {
					return(1);
				}
	}	}	}
	return (0); /*if after checking all neighboring edges there has been no */
			/*occurance of a close pair such as 2 and 3 */
			/*then the section being considered does not connect to phi */
			/*recall that all edges are filled! see if statement above */
}
