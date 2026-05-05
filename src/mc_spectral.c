#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include "mc_runtime.h"
#include "mc_spectral.h"

#define MC_SPECTRAL_MAXIT_2SAP 1000
#define MC_SPECTRAL_MAXIT_2SAP_HAM 10000
#define MC_SPECTRAL_HAM_SHIFT 100.0
#define MC_SPECTRAL_HAM_ACCURACY 0.0000000001

typedef struct McSpectralCsr {
	unsigned long int rows;
	unsigned long int nnz;
	unsigned long int *row_ptr;
	unsigned long int *cols;
	unsigned long int *edges;
	unsigned long int max_edge;
} McSpectralCsr;

static McSpectralCsr build_transition_csr(const McTransitionSpectralInput *input)
{
	McSpectralCsr csr;
	unsigned long int insection;
	unsigned long int nth;
	unsigned long int row = 0;
	unsigned long int pos = 0;

	memset(&csr, 0, sizeof(csr));
	csr.rows = input->max_tspans;
	csr.row_ptr = mc_xcalloc(csr.rows + 2, sizeof(*csr.row_ptr), "2SAP spectral CSR row pointers");

	for (insection = 1; insection <= input->max_keynum; insection++) {
		for (nth = 1; nth <= input->num_outsections[insection]; nth++) {
			unsigned long int outsection = input->tspans_outsection[insection][nth];
			row++;
			if (row > csr.rows) {
				fprintf(stderr, "Fatal: 2SAP spectral row count exceeded max_tspans\n");
				exit(EXIT_FAILURE);
			}
			csr.row_ptr[row + 1] = csr.row_ptr[row] + input->num_outsections[outsection];
			if (input->tspans_edges[insection][nth] > csr.max_edge) {
				csr.max_edge = input->tspans_edges[insection][nth];
			}
		}
	}
	if (row != csr.rows) {
		fprintf(stderr, "Fatal: 2SAP spectral row count %lu did not match max_tspans %lu\n", row, csr.rows);
		exit(EXIT_FAILURE);
	}

	csr.nnz = csr.row_ptr[csr.rows + 1];
	csr.cols = mc_xcalloc(csr.nnz + 1, sizeof(*csr.cols), "2SAP spectral CSR columns");
	csr.edges = mc_xcalloc(csr.nnz + 1, sizeof(*csr.edges), "2SAP spectral CSR edges");

	row = 0;
	for (insection = 1; insection <= input->max_keynum; insection++) {
		for (nth = 1; nth <= input->num_outsections[insection]; nth++) {
			unsigned long int outsection = input->tspans_outsection[insection][nth];
			unsigned long int right_tspan;
			unsigned long int edge_count = input->tspans_edges[insection][nth];
			row++;
			for (right_tspan = 1; right_tspan <= input->num_outsections[outsection]; right_tspan++) {
				pos++;
				csr.cols[pos] = input->tspans_nrr[outsection][right_tspan];
				csr.edges[pos] = edge_count;
			}
		}
	}

	return csr;
}

static void free_transition_csr(McSpectralCsr *csr)
{
	free(csr->row_ptr);
	free(csr->cols);
	free(csr->edges);
	memset(csr, 0, sizeof(*csr));
}

static double *build_fugacity_powers(double fugacity, double force, unsigned long int max_edge)
{
	unsigned long int i;
	double exp_force = exp(force);
	double *powers = mc_xmalloc((max_edge + 1) * sizeof(*powers), "2SAP spectral fugacity powers");

	for (i = 0; i <= max_edge; i++) {
		powers[i] = pow(fugacity, i) * exp_force;
	}
	return powers;
}

static void multiply_transition_matrix(
	const McSpectralCsr *csr,
	const double *powers,
	const double *left_in,
	double *left_out,
	const double *right_in,
	double *right_out)
{
	unsigned long int row;

	#pragma omp parallel for schedule(static)
	for (row = 1; row <= csr->rows; row++) {
		unsigned long int pos;
		unsigned long int start = csr->row_ptr[row] + 1;
		unsigned long int end = csr->row_ptr[row + 1];
		double right_sum = 0.0;

		for (pos = start; pos <= end; pos++) {
			unsigned long int col = csr->cols[pos];
			double weight = powers[csr->edges[pos]];

			#pragma omp atomic update
			left_out[col] += left_in[row] * weight;
			right_sum += right_in[col] * weight;
		}
		right_out[row] += right_sum;
	}
}

static double max_eval_regular_2sap(const McTransitionSpectralInput *input, const McSpectralCsr *csr, const double *powers, double fugacity)
{
	int k;
	unsigned long int i;
	double L_lambda_0;
	double L_lambda_1;
	double R_lambda_0;
	double R_lambda_1;
	double *left0 = input->left[0];
	double *left1 = input->left[1];
	double *right0 = input->right[0];
	double *right1 = input->right[1];

	printf("\nUsing OpenMP CSR Power Method with fugacity = %f\n\n", fugacity);
	printf("Iteration     Left eigenvalue     Difference     Right eigenvalue     Difference\n\n");

	#pragma omp parallel for schedule(static)
	for (i = 1; i <= input->max_tspans; i++) {
		left0[i] = 1.0;
		right0[i] = 1.0;
	}

	for (k = 1; k <= MC_SPECTRAL_MAXIT_2SAP; k++) {
		L_lambda_0 = 0.1;
		R_lambda_0 = 0.1;
		memset(left1 + 1, 0, input->max_tspans * sizeof(*left1));
		memset(right1 + 1, 0, input->max_tspans * sizeof(*right1));
		multiply_transition_matrix(csr, powers, left0, left1, right0, right1);

		#pragma omp parallel for reduction(max:L_lambda_0, R_lambda_0) schedule(static)
		for (i = 1; i <= input->max_tspans; i++) {
			if (left1[i] > L_lambda_0) {
				L_lambda_0 = left1[i];
			}
			if (right1[i] > R_lambda_0) {
				R_lambda_0 = right1[i];
			}
		}
		#pragma omp parallel for schedule(static)
		for (i = 1; i <= input->max_tspans; i++) {
			left1[i] /= L_lambda_0;
			right1[i] /= R_lambda_0;
		}

		L_lambda_1 = 0.1;
		R_lambda_1 = 0.1;
		memset(left0 + 1, 0, input->max_tspans * sizeof(*left0));
		memset(right0 + 1, 0, input->max_tspans * sizeof(*right0));
		multiply_transition_matrix(csr, powers, left1, left0, right1, right0);

		#pragma omp parallel for reduction(max:L_lambda_1, R_lambda_1) schedule(static)
		for (i = 1; i <= input->max_tspans; i++) {
			if (left0[i] > L_lambda_1) {
				L_lambda_1 = left0[i];
			}
			if (right0[i] > R_lambda_1) {
				R_lambda_1 = right0[i];
			}
		}
		#pragma omp parallel for schedule(static)
		for (i = 1; i <= input->max_tspans; i++) {
			left0[i] /= L_lambda_1;
			right0[i] /= R_lambda_1;
		}

		if (fabs(L_lambda_1 - L_lambda_0) < 0.000000001 && fabs(R_lambda_1 - R_lambda_0) < 0.000000001) {
			printf("\n%6i (DONE)", 2 * k);
			printf("%13f %17f  ", L_lambda_1, L_lambda_1 - L_lambda_0);
			printf("%16f %17f\n", R_lambda_1, R_lambda_1 - R_lambda_0);
			return R_lambda_1 - 1.0;
		}

		printf("%6i", 2 * k);
		printf("%20f %17f  ", L_lambda_1, L_lambda_1 - L_lambda_0);
		printf("%16f %17f\n", R_lambda_1, R_lambda_1 - R_lambda_0);
	}

	fprintf(stderr, "Maximum number of iterations exceeded in power method\n");
	exit(EXIT_FAILURE);
}

static double max_eval_ham_2sap(const McTransitionSpectralInput *input, const McSpectralCsr *csr, const double *powers, double fugacity)
{
	int k;
	unsigned long int i;
	int use_shift = (((input->L + 1) * (input->M + 1)) % 2) != 0;
	double *left0 = input->left[0];
	double *left1 = input->left[1];
	double *right0 = input->right[0];
	double *right1 = input->right[1];

	printf("\nUsing OpenMP CSR Modified Power Method with fugacity = %f\n\n", fugacity);

	#pragma omp parallel for schedule(static)
	for (i = 1; i <= input->max_tspans; i++) {
		left0[i] = 1.0;
		right0[i] = 1.0;
	}

	for (k = 1; k <= MC_SPECTRAL_MAXIT_2SAP_HAM; k++) {
		double lambda_min;
		double lambda_max;
		double min_ratio;
		double max_ratio;
		double themax_R;
		double themax_L;
		double difference;

		memset(left1 + 1, 0, input->max_tspans * sizeof(*left1));
		memset(right1 + 1, 0, input->max_tspans * sizeof(*right1));
		multiply_transition_matrix(csr, powers, left0, left1, right0, right1);
		if (use_shift) {
			#pragma omp parallel for schedule(static)
			for (i = 1; i <= input->max_tspans; i++) {
				left1[i] += MC_SPECTRAL_HAM_SHIFT * left0[i];
				right1[i] += MC_SPECTRAL_HAM_SHIFT * right0[i];
			}
		}

		themax_L = left1[1];
		themax_R = right1[1];
		min_ratio = right1[1] / right0[1];
		max_ratio = right1[1] / right0[1];
		for (i = 2; i <= input->max_tspans; i++) {
			double cur_ratio = right1[i] / right0[i];
			if (cur_ratio > max_ratio) {
				max_ratio = cur_ratio;
			}
			if (cur_ratio < min_ratio) {
				min_ratio = cur_ratio;
			}
			if (right1[i] > themax_R) {
				themax_R = right1[i];
			}
			if (left1[i] > themax_L) {
				themax_L = left1[i];
			}
		}
		lambda_min = min_ratio;
		lambda_max = max_ratio;
		printf("Iteration %d. lambda_min=%f, lambda_max=%f.  ", k, lambda_min, lambda_max);

		#pragma omp parallel for schedule(static)
		for (i = 1; i <= input->max_tspans; i++) {
			right1[i] /= themax_R;
			left1[i] /= themax_L;
		}

		difference = lambda_max - lambda_min;
		printf("Difference is %f = %e\n", difference, difference);

		if (difference < MC_SPECTRAL_HAM_ACCURACY) {
			if (use_shift) {
				return lambda_max - MC_SPECTRAL_HAM_SHIFT - 1.0;
			}
			return lambda_max - 1.0;
		}

		memset(left0 + 1, 0, input->max_tspans * sizeof(*left0));
		memset(right0 + 1, 0, input->max_tspans * sizeof(*right0));
		multiply_transition_matrix(csr, powers, left1, left0, right1, right0);
		if (use_shift) {
			#pragma omp parallel for schedule(static)
			for (i = 1; i <= input->max_tspans; i++) {
				left0[i] += MC_SPECTRAL_HAM_SHIFT * left1[i];
				right0[i] += MC_SPECTRAL_HAM_SHIFT * right1[i];
			}
		}

		themax_L = left0[1];
		themax_R = right0[1];
		min_ratio = right0[1] / right1[1];
		max_ratio = right0[1] / right1[1];
		for (i = 2; i <= input->max_tspans; i++) {
			double cur_ratio = right0[i] / right1[i];
			if (cur_ratio > max_ratio) {
				max_ratio = cur_ratio;
			}
			if (cur_ratio < min_ratio) {
				min_ratio = cur_ratio;
			}
			if (right0[i] > themax_R) {
				themax_R = right0[i];
			}
			if (left0[i] > themax_L) {
				themax_L = left0[i];
			}
		}
		lambda_min = min_ratio;
		lambda_max = max_ratio;
		printf("Iteration %d. lambda_min=%f, lambda_max=%f.  ", k, lambda_min, lambda_max);

		#pragma omp parallel for schedule(static)
		for (i = 1; i <= input->max_tspans; i++) {
			right0[i] /= themax_R;
			left0[i] /= themax_L;
		}

		difference = lambda_max - lambda_min;
		printf("Difference is %f = %e\n", difference, difference);

		if (difference < MC_SPECTRAL_HAM_ACCURACY) {
			if (use_shift) {
				return lambda_max - MC_SPECTRAL_HAM_SHIFT - 1.0;
			}
			return lambda_max - 1.0;
		}
	}

	fprintf(stderr, "Maximum number of iterations exceeded in power method\n");
	exit(EXIT_FAILURE);
}

double mc_spectral_max_eval_2sap(const McTransitionSpectralInput *input, double fugacity)
{
	double result;
	double *powers;
	McSpectralCsr csr = build_transition_csr(input);

	powers = build_fugacity_powers(fugacity, input->force, csr.max_edge);
	if (input->hamiltonian) {
		result = max_eval_ham_2sap(input, &csr, powers, fugacity);
	} else {
		result = max_eval_regular_2sap(input, &csr, powers, fugacity);
	}

	free(powers);
	free_transition_csr(&csr);
	return result;
}

#undef MC_SPECTRAL_MAXIT_2SAP
#undef MC_SPECTRAL_MAXIT_2SAP_HAM
#undef MC_SPECTRAL_HAM_SHIFT
#undef MC_SPECTRAL_HAM_ACCURACY
