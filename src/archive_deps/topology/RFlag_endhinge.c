unsigned short int		/* returns 1 if testhinge connects, 0 otherwise */
RFlag_endhinge(int (*EndOrdNum)[2])
{
	if ((*EndOrdNum)[0] - 1 < vM*vL) {	/* if there are less edges in first section */
						/* than there are number of vertices */
						/* ie. if first section is not full */

/* Strictly speaking this should actually read (*pointordNum) - 1 < vM*vL */
/* But since (*pointordNum) is odd and vM*vL is even it is not necessary */

		return (1);

	} else { /*all edges are occupied in this section */
		unsigned short int i, j,b; /*dummy variables*/

		for (i = 0; i <= M; i++) {
			for (j = 0; j < L; j++) {/*the i and j loops serve to loop through*/
						/* all posible pairs of horizontally neighboring edges */
						/* that is i,j and i,j+1  hence the j loop ends at L-1 */

				switch ((endtemplate[i][j]) - (b=endtemplate[i][j + 1])) {


				case (1):{
						if ( (b & 1) == 1) /*b is smaller and congruent to 1 mod 2*/
						return(1);
						break;
						/*these conditions are satisfied if 2 neighbors 1 */
						/* or 4 neighbors 3, or 6 neighbors 5, etc... */
						/* note: it doesn't make sense for 0 to neighbor 1 */
						/* but this case is elliminated implicitely since 1 > 0 */ 
					}
				case (-1):{
						if ( (b & 1) == 0) /*b is larger and congruent to 0 mod 2*/
						return(1);
						break;
						/*these conditions are satisfied if 1 neighbors 2 */
						/* or 3 neighbors 4, or 5 neighbors 6, etc... */
						/* note: it doesn't make sense for 1 to neighbor 0 */
						/* but this case is elliminated implicitely since 0 < 1 */
					}
				default:{
						break;   /*the numbers of the two edges do not differ by one */
							 /*no immediate way for the the edges to join together */
					}
				}

			}
		}
		for (i = 0; i < M; i++) {

			for (j = 0; j <= L; j++) {/*the i and j loops serve to loop through*/
						/* all posible pairs of vertically neighboring edges */
						/* that is i,j and i+1,j  hence the i loop ends at M-1 */

				switch ((endtemplate[i][j]) - (b=endtemplate[i+1][j])) {
				case (1):{
						if ( (b & 1) == 1) /*b is smaller and congruent to 1 mod 2*/
						return(1);
						break;
						/*these conditions are satisfied if 2 neighbors 1 */
						/* or 4 neighbors 3, or 6 neighbors 5, etc... */
						/* note: it doesn't make sense for 0 to neighbor 1 */
						/* but this case is elliminated implicitely since 1 > 0 */ 
					}
				case (-1):{
						if ( (b & 1) == 0) /*b is larger and congruent to 0 mod 2*/
						return(1);
						break;
						/*these conditions are satisfied if 1 neighbors 2 */
						/* or 3 neighbors 4, or 5 neighbors 6, etc... */
						/* note: it doesn't make sense for 1 to neighbor 0 */
						/* but this case is elliminated implicitely since 0 < 1 */

					}
				default:{
						break;  /*the numbers of the two edges do not differ by one */
							 /*no immediate way for the the edges to join together */

					}
				}
	}	}	}
	return (0);	/*if after checking all neighboring edges there has been no */
			/*occurance of a close pair such as 1 and 2 */
			/*then the section being considered does not connect to phi */
			/*recall that all edges are filled! see if statement above */
}
