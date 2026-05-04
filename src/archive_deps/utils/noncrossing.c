unsigned short int
noncrossing(int side, int a, int b, int c, int d)
{

	unsigned short int j;	/* dummy variable to loop through each edge */
	unsigned short int crossing_signature = 0; /* a number representing the implied way */
						   /* in which */
	unsigned long int n = 1;	/* bit holder */


	for (j = 0; j <= L; j++) {	/* loop through each edge in the slit */

		if (ordertemplate[side][0][j] == a || ordertemplate[side][0][j] == b) {
		/* if edge is numbered with a number from the first pair then */
		/**/ 
		/* by not adding n to crossing_signature this effectively puts a zero */
		/* in the bit place determined by n */
			n <<= 1; /* now shift to the next higher bit */
		}
		if (ordertemplate[side][0][j] == c || ordertemplate[side][0][j] == d) {
		/* if edge is numbered with number from the second pair then */
		/**/
		/* record a one in the bit location determined by n */
			crossing_signature += n;
			n <<= 1; /* now shift n to the next higher bit */
		}
	}

	/* Basically this function goes along the slit and records a zero when a edge numbered */
	/* from the first pair is encountered, and records a one when a edge numbered from the */
	/* second pair is encountered */
	/**/
	/* example: consider the slit (on side 1) 	 			*/
	/*				  0 1 3 0 2 4 		                */
	/* if noncrossing is called with parameters (1,1,2,3,4) then the result */
	/* is that going from left to right 1 is encountered first then 3, then 2 */
	/* then 4 */
	/* so in binary code the crossing signature is 1010 (this is built up from right to left) */
	/* that is the zero is for the edge numbered one from the first pair which is encountered first */
	/* in the for loop through the vertices. The one is for the edge numbered three from */
	/* the second pair which is encountered second in the for loop etc... */
	/* note that there is an implied crossing of paths so that this is not a valid */
	/* slit. That is the path between 1 and 2 crosses the path between 3 and 4 */
	/* The binary number 1010 is the number 10 in decimal notation thus */
	/* if crossing signature is 10 the slit is not valid */
	/* also if crossing signature is 5 the slit is not valid */
	/* 5 --> 0101 */
	/* We also cannot have a crossing signature of 9 if we are on the left side of the two-span */
	/* and the edge numbered 1 is in the first pair */
	/* Note: when this function is called it is assumed that if there is a edge numbered 1 */
	/* then it will be the second parameter "a" */
	/* The reason this crossing signature is not valid is because */
	/* 9 --> 1 0 0 1  and this would imply that this section cannot connect to phi */
	/* an example of a slit that would produce this crossing signature is */
	/* 4 1 2 3  assuming noncrossing is called with parameters (0,1,2,3,4)*/	

	if (crossing_signature == 10 || crossing_signature == 5 || (side == 0 && a == 1 && crossing_signature == 9) ) {
		return (0); /* the pairs a,b cross the pairs c,d */
		            /* or the slit is not valid because of location of edge numbered 1 */
	} else
		return (1); /* the pairs do not cross each other */
}
