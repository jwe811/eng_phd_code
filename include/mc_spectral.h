#ifndef MC_SPECTRAL_H
#define MC_SPECTRAL_H

typedef struct McTransitionSpectralInput {
	/* The 2SAP sampler solves on transitions, not sections. Each two-span
	   transition is one spectral row; outgoing columns are the transitions that
	   can legally follow it. */
	unsigned long int max_keynum;
	unsigned long int max_tspans;
	unsigned long int *num_outsections;
	unsigned long int **tspans_outsection;
	unsigned long int **tspans_edges;
	unsigned long int **tspans_nrr;
	double *left[2];
	double *right[2];
	double force;
	int L;
	int M;
	int hamiltonian;
} McTransitionSpectralInput;

double mc_spectral_max_eval_2sap(const McTransitionSpectralInput *input, double fugacity);

#endif
