#ifndef TM_SPECTRAL_H
#define TM_SPECTRAL_H

typedef struct TmSpectralProblem {
	/* CSR matrix for the transfer graph. Rows and vector entries are 1-based
	   because the archival vector/matrix helpers reserve index 0. */
	unsigned long int states;
	unsigned long int *row_ptr;
	unsigned long int *out_states;
	unsigned long int *edges;
	unsigned long int max_edge;
	/* Two buffers each for ping-pong power-method iterations. */
	double *left[2];
	double *right[2];
	double convergence_threshold;
	int damping_enabled;
	/* Transpose CSR lets the left-vector multiply run without OpenMP atomics. */
	unsigned long int *transpose_row_ptr;
	unsigned long int *transpose_rows;
	unsigned long int *transpose_edges;
} TmSpectralProblem;

void tm_spectral_prepare_transpose(TmSpectralProblem *problem);
void tm_spectral_free(TmSpectralProblem *problem);
double tm_spectral_max_eval(TmSpectralProblem *problem, double fugacity);
double tm_spectral_beta(const TmSpectralProblem *problem, double kappa);

#endif
