unsigned short int		/* returns 1 if testhinge connects, 0 otherwise */
LFlag_endhinge(int (*EndOrdNum)[2])
{
	if ((*EndOrdNum)[0] - 1 < vM*vL) {	/* if there are less edges in first section */
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

				switch ((endtemplate[i][j]) - (b=endtemplate[i][j + 1])) {
					/* look at the difference between the numbers of the two edges */

				case (1):{	if ( b!= 0 && (b & 1) == 0) /*b is not zero but is congruent to 0 mod 2*/
						return(1); break; }
						/*these conditions are satisfied if 3 neighbors 2 */
						/* or 5 neighbors 4, or 7 neighbors 6, etc... */
						/* note: it doesn't make sense for 1 to neighbor 0 */
						/* so this case is elliminated by  b != 0 in if statement*/

				case (-1):{	if ( b!=1 && (b & 1) == 1) /*b is not one but is congruent to 1 mod 2*/
						return(1); break; }
						/*these conditions are satisfied if 2 neighbors 3 */
						/* or 4 neighbors 5, or 6 neighbors 7, etc... */
						/* note: it doesn't make sense for 0 to neighbor 1 */
						/* so this case is elliminated by b != 1 in if statement*/
				case (vM*vL-1):{ return(1); }
				case (1-vM*vL):{ return(1); }

				default:{	break; } /*the numbers of the two edges do not differ by one */
							 /*no immediate way for the the edges to join together */

				} /*end switch*/
			}
		}

		for (i = 0; i < M; i++) {
			for (j = 0; j <= L; j++) {/*the i and j loops serve to loop through*/
						/* all posible pairs of vertically neighboring edges */
						/* that is i,j and i+1,j  hence the i loop ends at M-1 */

				switch ((endtemplate[i][j]) - (b=endtemplate[i+1][j])) {

				case (1):{	if ( b!= 0 && (b & 1) == 0) /*b is not zero but is congruent to 0 mod 2*/
						return(1); break; }
						/*these conditions are satisfied if 3 neighbors 2 */
						/* or 5 neighbors 4, or 7 neighbors 6, etc... */
						/* note: it doesn't make sense for 1 to neighbor 0 */
						/* so this case is elliminated by b != 0 in if statement*/

				case (-1):{	if ( b!=1 && (b & 1) == 1) /*b is not one but is congruent to 1 mod 2*/
						return(1); break; }
						/*these conditions are satisfied if 2 neighbors 3 */
						/* or 4 neighbors 5, or 6 neighbors 7, etc... */
						/* note: it doesn't make sense for 0 to neighbor 1 */
						/* so this case is elliminated by b != 1 in if statement*/
				case (vM*vL-1):{ return(1); }
				case (1-vM*vL):{ return(1); }

				default:{	break; } /*the numbers of the two edges do not differ by one */
							 /*no immediate way for the the edges to join together */

				} /*end switch*/
	}	}	}
	return (0); /*if after checking all neighboring edges there has been no */
			/*occurance of a close pair such as 2 and 3 */
			/*then the section being considered does not connect to phi */
			/*recall that all edges are filled! see if statement above */
}
