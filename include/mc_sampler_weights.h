#ifndef MC_SAMPLER_WEIGHTS_H
#define MC_SAMPLER_WEIGHTS_H

typedef struct McSamplerWeightInput {
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
	unsigned long int left_count;
	double *left_acceptance;
	unsigned long int *left_section;
	int *left_endhinge;
	double max_left_acceptance;
	double *right_acceptance;
	double max_right_acceptance;
} McSamplerWeights;

void mc_sampler_weights_build(const McSamplerWeightInput *input, McSamplerWeights *weights);
void mc_sampler_weights_free(McSamplerWeights *weights);

#endif
