#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mc_runtime.h"
#include "mc_sampler_weights.h"

void mc_sampler_weights_build(const McSamplerWeightInput *input, McSamplerWeights *weights)
{
	unsigned long int ordinal = 0;
	int section;
	int nth;

	memset(weights, 0, sizeof(*weights));
	weights->left_count = input->total_left_endhinges;
	if (weights->left_count == 0) {
		fprintf(stderr, "Fatal: no left endhinges available for sampler weights\n");
		exit(EXIT_FAILURE);
	}
	weights->left_acceptance = mc_xcalloc(weights->left_count, sizeof(*weights->left_acceptance), "left endhinge acceptance weights");
	weights->left_section = mc_xcalloc(weights->left_count, sizeof(*weights->left_section), "left endhinge section map");
	weights->left_endhinge = mc_xcalloc(weights->left_count, sizeof(*weights->left_endhinge), "left endhinge ordinal map");
	weights->right_acceptance = mc_xcalloc((size_t)input->max_tspans + 1, sizeof(*weights->right_acceptance), "right endhinge acceptance weights");

	/* For a left endhinge in section s, the acceptance envelope is the sum of
	   destination right-eigenvector weights over every outgoing two-span. */
	for (section = 1; section <= input->max_sections; section++) {
		for (nth = 1; nth <= input->num_left_endhinges[section]; nth++) {
			int out_idx;
			double value = 0.0;
			if (ordinal >= weights->left_count) {
				fprintf(stderr, "Fatal: left endhinge ordinal overflow while building sampler weights\n");
				exit(EXIT_FAILURE);
			}
			for (out_idx = 1; out_idx <= input->num_outsections[section]; out_idx++) {
				unsigned long int tspan_num = input->t_nrr[section][out_idx];
				double right_value;
				if (tspan_num > (unsigned long int)input->max_tspans) {
					fprintf(stderr, "Fatal: tspan number %lu exceeds max_tspans %d while building left sampler weights\n",
						tspan_num, input->max_tspans);
					exit(EXIT_FAILURE);
				}
				right_value = input->right_evector[tspan_num];
				if (right_value < 0.0 || !isfinite(right_value)) {
					fprintf(stderr, "Fatal: invalid right eigenvector value for tspan %lu: %.17g\n",
						tspan_num, right_value);
					exit(EXIT_FAILURE);
				}
				value += right_value;
			}
			if (!isfinite(value) || value < 0.0) {
				fprintf(stderr, "Fatal: invalid left acceptance weight for section %d endhinge %d: %.17g\n",
					section, nth, value);
				exit(EXIT_FAILURE);
			}
			weights->left_acceptance[ordinal] = value;
			weights->left_section[ordinal] = (unsigned long int)section;
			weights->left_endhinge[ordinal] = nth;
			if (value > weights->max_left_acceptance) {
				weights->max_left_acceptance = value;
			}
			ordinal++;
		}
	}

	if (ordinal != weights->left_count) {
		fprintf(stderr, "Fatal: left endhinge count mismatch while building sampler weights (saw %lu, expected %lu)\n",
			ordinal, weights->left_count);
		exit(EXIT_FAILURE);
	}
	if (weights->max_left_acceptance <= 0.0 || !isfinite(weights->max_left_acceptance)) {
		fprintf(stderr, "Fatal: invalid maximum left acceptance weight: %.17g\n", weights->max_left_acceptance);
		exit(EXIT_FAILURE);
	}

	/* For the final cap, the envelope depends on the chosen transition: it is
	   the number of right endhinges available at the output section divided by
	   that transition's eigenvector weight. */
	for (section = 1; section <= input->max_sections; section++) {
		for (nth = 1; nth <= input->num_outsections[section]; nth++) {
			unsigned long int tspan_num = input->t_nrr[section][nth];
			double right_value;
			double value;
			if (tspan_num > (unsigned long int)input->max_tspans) {
				fprintf(stderr, "Fatal: tspan number %lu exceeds max_tspans %d while building sampler weights\n",
					tspan_num, input->max_tspans);
				exit(EXIT_FAILURE);
			}
			right_value = input->right_evector[tspan_num];
			if (right_value <= 0.0 || !isfinite(right_value)) {
				fprintf(stderr, "Fatal: invalid right eigenvector value for tspan %lu: %.17g\n",
					tspan_num, right_value);
				exit(EXIT_FAILURE);
			}
			value = input->num_right_endhinges[input->t_outsection[section][nth]] / right_value;
			if (!isfinite(value) || value < 0.0) {
				fprintf(stderr, "Fatal: invalid right acceptance weight for section %d tspan %d: %.17g\n",
					section, nth, value);
				exit(EXIT_FAILURE);
			}
			weights->right_acceptance[tspan_num] = value;
			if (value > weights->max_right_acceptance) {
				weights->max_right_acceptance = value;
			}
		}
	}
	if (weights->max_right_acceptance <= 0.0 || !isfinite(weights->max_right_acceptance)) {
		fprintf(stderr, "Fatal: invalid maximum right acceptance weight: %.17g\n", weights->max_right_acceptance);
		exit(EXIT_FAILURE);
	}
}

void mc_sampler_weights_free(McSamplerWeights *weights)
{
	free(weights->left_acceptance);
	free(weights->left_section);
	free(weights->left_endhinge);
	free(weights->right_acceptance);
	memset(weights, 0, sizeof(*weights));
}
