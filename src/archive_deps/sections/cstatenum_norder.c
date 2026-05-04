/***************************************************************************/
unsigned long int
num_section_norder(int side)
{				/* f is the number of arc ends */
	unsigned long int i, j, jj;	/* used as dummy variables in for
					 * loops */
	unsigned short int edge_num, hsum;	/* the number of the edge
						 * being considered in for
						 * loop */
	/* this is zero if there is no edge there */
	unsigned short int seenpair[20] = {0};
	unsigned short int h[20] = {0};
	unsigned short int c[20] = {0};
	unsigned short int n[20] = {0};
	unsigned long int lexicord = 1;
	unsigned short int a, f, q, p, m;

	unsigned short int nchoosem[12][12] = {
		{1},
		{1, 1},
		{1, 2, 1},
		{1, 3, 3, 1},
		{1, 4, 6, 4, 1},
		{1, 5, 10, 10, 5, 1},
		{1, 6, 15, 20, 15, 6, 1},
		{1, 7, 21, 35, 35, 21, 7, 1},
		{1, 8, 28, 56, 70, 56, 28, 8, 1},
		{1, 9, 36, 84, 126, 126, 84, 36, 9, 1},
		{1,10,45,120,210,252,210,120,45,10,1},
		{1,11,55,165,330,462,462,330,165,55,11,1}
	};


	f = 0;
	for (j = 0; j <= L; j++) {
		if (ordertemplate[side][0][j] != 0) {
			f++;
		}
	}
	a = f / 2;
	jj = 0;

	/* Necessary to initialize the array to contain zero elements ? */

	for (j = 0; j <= L; j++) {
		edge_num = ordertemplate[side][0][j];	/* assign number of edge
							 * determined by i and j
							 * to edge_num */
		/* this will be zero if there is no edge */

		if (edge_num != 0) {	/* if there is an edge at this
					 * position */
			jj++;
			if (seenpair[edge_num] == 0) {
				h[jj] = 1;
			}
			seenpair[edge_num] = 1;
			seenpair[paircomp(edge_num, f)] = 1;
		}
	}
	jj = 0;
	hsum = 1;
	for (i = 2; i <= f - 1; i++) {
		if (hsum + h[i] > hsum) {
			jj++;
			n[a - jj] = 2 * jj + 1 - i;
		}
		hsum += h[i];
	}
	jj = 0;
	for (j = 0; j <= L; j++) {
		if (ordertemplate[side][0][j] != 0) {
			jj++;
			c[jj] = j + 1;
		}
	}

	for (m = 1; m <= a - 1; m++) {
		lexicord += nchoosem[2 * m][m] * nchoosem[L + 1][2 * m] / (m + 1);
	}
	for (m = 1; m <= a - 1; m++) {
		lexicord += nchoosem[L + 1][2*a]*nchoosem[2 * m + n[m] - 1][m] * n[m] / (n[m] + m);
	}
	for (q = 1; q <= f; q++) {
		for (p = c[q - 1] + 1; p <= c[q] - 1; p++) {
			lexicord += nchoosem[L + 1 - p][2 * a - q];
		}
	}
	return (lexicord);

}				/* end of function num_section */
/***************************************************************************/

