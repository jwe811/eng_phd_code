#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include "tm_spectral.h"

static void spectral_fatal(const char *message)
{
	fprintf(stderr, "Fatal: %s\n", message);
	exit(EXIT_FAILURE);
}

static void *spectral_xcalloc(size_t count, size_t size, const char *label)
{
	void *ptr = calloc(count, size);
	if (!ptr) {
		fprintf(stderr, "Fatal: unable to allocate %s (%zu x %zu bytes)\n", label, count, size);
		exit(EXIT_FAILURE);
	}
	return ptr;
}

static void *spectral_xmalloc(size_t size, const char *label)
{
	void *ptr = malloc(size);
	if (!ptr) {
		fprintf(stderr, "Fatal: unable to allocate %s (%zu bytes)\n", label, size);
		exit(EXIT_FAILURE);
	}
	return ptr;
}

static double *build_fugacity_powers(double fugacity, unsigned long int max_edge)
{
	unsigned long int i;
	double *powers = spectral_xmalloc((max_edge + 1) * sizeof(*powers), "fugacity powers");

	for (i = 0; i <= max_edge; i++) {
		powers[i] = pow(fugacity, i);
	}
	return powers;
}

void tm_spectral_prepare_transpose(TmSpectralProblem *problem)
{
	unsigned long int row;
	unsigned long int pos;
	unsigned long int total_transitions = problem->row_ptr[problem->states + 1];
	unsigned long int *write_cursor;

	problem->max_edge = 0;
	problem->transpose_row_ptr = spectral_xcalloc(problem->states + 2, sizeof(*problem->transpose_row_ptr), "transpose CSR row pointer");
	problem->transpose_rows = spectral_xcalloc(total_transitions + 1, sizeof(*problem->transpose_rows), "transpose CSR source rows");
	problem->transpose_edges = spectral_xcalloc(total_transitions + 1, sizeof(*problem->transpose_edges), "transpose CSR edge weights");

	for (row = 1; row <= problem->states; row++) {
		unsigned long int start = problem->row_ptr[row] + 1;
		unsigned long int end = problem->row_ptr[row + 1];
		for (pos = start; pos <= end; pos++) {
			unsigned long int col = problem->out_states[pos];
			if (col < 1 || col > problem->states) {
				spectral_fatal("CSR transition points outside state range");
			}
			problem->transpose_row_ptr[col + 1]++;
			if (problem->edges[pos] > problem->max_edge) {
				problem->max_edge = problem->edges[pos];
			}
		}
	}

	for (row = 1; row <= problem->states; row++) {
		problem->transpose_row_ptr[row + 1] += problem->transpose_row_ptr[row];
	}

	write_cursor = spectral_xmalloc((problem->states + 2) * sizeof(*write_cursor), "transpose CSR write cursors");
	memcpy(write_cursor, problem->transpose_row_ptr, (problem->states + 2) * sizeof(*write_cursor));

	for (row = 1; row <= problem->states; row++) {
		unsigned long int start = problem->row_ptr[row] + 1;
		unsigned long int end = problem->row_ptr[row + 1];
		for (pos = start; pos <= end; pos++) {
			unsigned long int col = problem->out_states[pos];
			unsigned long int out_pos = ++write_cursor[col];
			problem->transpose_rows[out_pos] = row;
			problem->transpose_edges[out_pos] = problem->edges[pos];
		}
	}

	free(write_cursor);
}

void tm_spectral_free(TmSpectralProblem *problem)
{
	free(problem->transpose_row_ptr);
	free(problem->transpose_rows);
	free(problem->transpose_edges);
	problem->transpose_row_ptr = NULL;
	problem->transpose_rows = NULL;
	problem->transpose_edges = NULL;
}

static void multiply_right(
	const TmSpectralProblem *problem,
	const double *powers,
	const double *right_in,
	double *right_out)
{
	unsigned long int row;

	#pragma omp parallel for schedule(static)
	for (row = 1; row <= problem->states; row++) {
		unsigned long int pos;
		unsigned long int start = problem->row_ptr[row] + 1;
		unsigned long int end = problem->row_ptr[row + 1];
		double sum = 0.0;
		for (pos = start; pos <= end; pos++) {
			sum += right_in[problem->out_states[pos]] * powers[problem->edges[pos]];
		}
		right_out[row] = sum;
	}
}

static void multiply_left(
	const TmSpectralProblem *problem,
	const double *powers,
	const double *left_in,
	double *left_out)
{
	unsigned long int col;

	#pragma omp parallel for schedule(static)
	for (col = 1; col <= problem->states; col++) {
		unsigned long int pos;
		unsigned long int start = problem->transpose_row_ptr[col] + 1;
		unsigned long int end = problem->transpose_row_ptr[col + 1];
		double sum = 0.0;
		for (pos = start; pos <= end; pos++) {
			sum += left_in[problem->transpose_rows[pos]] * powers[problem->transpose_edges[pos]];
		}
		left_out[col] = sum;
	}
}

double tm_spectral_max_eval(TmSpectralProblem *problem, double fugacity)
{
	int k;
	unsigned long int i;
	double L0, L1, R0, R1;
	double *fug_pow = build_fugacity_powers(fugacity, problem->max_edge);

	#pragma omp parallel for schedule(static)
	for (i = 1; i <= problem->states; i++) {
		problem->left[0][i] = 1.0;
		problem->right[0][i] = 1.0;
	}

	for (k = 1; k <= 500; k++) {
		L0 = 0.1;
		R0 = 0.1;
		multiply_left(problem, fug_pow, problem->left[0], problem->left[1]);
		multiply_right(problem, fug_pow, problem->right[0], problem->right[1]);

		#pragma omp parallel for reduction(max:L0, R0) schedule(static)
		for (i = 1; i <= problem->states; i++) {
			if (problem->left[1][i] > L0) L0 = problem->left[1][i];
			if (problem->right[1][i] > R0) R0 = problem->right[1][i];
		}
		#pragma omp parallel for schedule(static)
		for (i = 1; i <= problem->states; i++) {
			double nL = problem->left[1][i] / L0;
			double nR = problem->right[1][i] / R0;
			if (problem->damping_enabled) {
				problem->left[1][i] = 0.5 * problem->left[0][i] + 0.5 * nL;
				problem->right[1][i] = 0.5 * problem->right[0][i] + 0.5 * nR;
			} else {
				problem->left[1][i] = nL;
				problem->right[1][i] = nR;
			}
		}

		L1 = 0.1;
		R1 = 0.1;
		multiply_left(problem, fug_pow, problem->left[1], problem->left[0]);
		multiply_right(problem, fug_pow, problem->right[1], problem->right[0]);

		#pragma omp parallel for reduction(max:L1, R1) schedule(static)
		for (i = 1; i <= problem->states; i++) {
			if (problem->left[0][i] > L1) L1 = problem->left[0][i];
			if (problem->right[0][i] > R1) R1 = problem->right[0][i];
		}
		#pragma omp parallel for schedule(static)
		for (i = 1; i <= problem->states; i++) {
			double nL = problem->left[0][i] / L1;
			double nR = problem->right[0][i] / R1;
			if (problem->damping_enabled) {
				problem->left[0][i] = 0.5 * problem->left[1][i] + 0.5 * nL;
				problem->right[0][i] = 0.5 * problem->right[1][i] + 0.5 * nR;
			} else {
				problem->left[0][i] = nL;
				problem->right[0][i] = nR;
			}
		}

		if (fabs(L1 - L0) < problem->convergence_threshold && fabs(R1 - R0) < problem->convergence_threshold) {
			free(fug_pow);
			return R1 - 1.0;
		}
	}

	free(fug_pow);
	spectral_fatal("power method did not converge");
	return 0.0;
}

double tm_spectral_beta(const TmSpectralProblem *problem, double kappa)
{
	double beta = 0.0;
	unsigned long int row;
	double *fug_pow = build_fugacity_powers(kappa, problem->max_edge);

	#pragma omp parallel for reduction(+:beta) schedule(static)
	for (row = 1; row <= problem->states; row++) {
		unsigned long int pos;
		unsigned long int start = problem->row_ptr[row] + 1;
		unsigned long int end = problem->row_ptr[row + 1];
		for (pos = start; pos <= end; pos++) {
			unsigned long int out_state = problem->out_states[pos];
			beta += problem->left[0][row] * problem->edges[pos] * fug_pow[problem->edges[pos]] * problem->right[0][out_state];
		}
	}

	free(fug_pow);
	return beta;
}
