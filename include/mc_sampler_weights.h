#ifndef MC_SAMPLER_WEIGHTS_H
#define MC_SAMPLER_WEIGHTS_H

typedef struct McSamplerWeightInput {
	/* Transition graph and endhinge counts produced by build_transfer_matrix_graph(). */
	int max_sections;
	int max_tspans;
	unsigned long int total_left_endhinges;
	unsigned long int *num_left_endhinges;
	unsigned long int *num_right_endhinges;
	unsigned long int *num_outsections;
	unsigned long int **t_outsection;
	unsigned long int **t_nrr;
	const double *right_evector;
} McSamplerWeightInput;

typedef struct McSamplerWeights {
	/* Dense ordinal -> (section, nth left endhinge) map. The sampler picks a
	   uniform ordinal, then accepts it using left_acceptance / max_left_acceptance. */
	unsigned long int left_count;
	double *left_acceptance;
	unsigned long int *left_section;
	int *left_endhinge;
	double max_left_acceptance;
	/* right_acceptance is indexed by transition number because the final cap is
	   chosen after the last two-span transition is known. */
	double *right_acceptance;
	double max_right_acceptance;
} McSamplerWeights;

void mc_sampler_weights_build(const McSamplerWeightInput *input, McSamplerWeights *weights);
void mc_sampler_weights_free(McSamplerWeights *weights);

#endif
