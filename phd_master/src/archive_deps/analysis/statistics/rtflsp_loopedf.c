/***************************************************************************/
/*From Numerical Recipes in C section 9.2*/
/* Using the false position method, find the root of a function known */
/* to lie between x1 and x2. The root, returned as rtflsp is refined  */
/* until its accuracy is +/- xacc. */

#define MAXIT 500 /* set the maximum allowed number of iterations */

double rtflsp(double(*func)(double,double),double x1,double x2,double xacc, double ff)
{
	int j;
	double fl,fh,xl,xh,swap,dx,del,f,rtf;

	fl=(*func)(x1,ff);
	fh=(*func)(x2,ff);
	if (fl*fh > 0.0) fprintf(stderr, "Root must be bracketed in RTFLSP\n");
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xl=x2;
		xh=x1;
		swap=fl;
		fl=fh;
		fh=swap;
	}
	dx=xh-xl;
	for (j=1;j<=MAXIT;j++) {
		rtf=xl+dx*fl/(fl-fh);
		f=(*func)(rtf,ff);
		if (f < 0.0) {
			del=xl-rtf;
			xl=rtf;
			fl=f;
		} else {
			del=xh-rtf;
			xh=rtf;
			fh=f;
		}
		dx=xh-xl;
		if (fabs(del) < xacc || f == 0.0) return rtf;
	}
	fprintf(stderr, "Maximum number of iterations exceeded in RTFLSP\n");
	exit(1);
}

#undef MAXIT

/***************************************************************************/
