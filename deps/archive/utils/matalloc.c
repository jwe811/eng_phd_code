/***************************************************************************/
double **
matalloc(int rowlow, int rowhigh, int collow, int colhigh)
{
	int             k;
	double           **x;
	x = (double **) calloc((unsigned) (rowhigh - rowlow + 1), sizeof(double *));
	if (x == NULL) {
		fprintf(stderr, "unable to allocate memory in function matalloc\n");
		exit(1);
	}
	x -= rowlow;
	for (k = rowlow; k <= rowhigh; k++) {
		x[k] = (double *) calloc((unsigned) (colhigh - collow + 1), sizeof(double));
		if (x[k] == NULL) {
			fprintf(stderr, "unable to allocate memory in function matalloc\n");
			exit(1);
		}
		x[k] -= collow;
	}
	return x;
}
/***************************************************************************/
