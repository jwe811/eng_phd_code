#include "mc_globals.h"
#include "mc_paths.h"
#include "mc_sampler_weights.h"
#include "run_metadata.h"
#include "../include/marsaglia.h"
#include <errno.h>
#include <limits.h>
#include <stdarg.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Monte Carlo sampler for modes 0/1. Modes 2/3 dispatch to the integrated
   2SAP entry points from mc_2sap_*_integrated.c, but they share the same CLI
   validation, RNG seed handling, and output conventions. */

typedef struct McThreadRng {
	unsigned int iran[43];
	unsigned int iran1;
	unsigned int ic;
	unsigned int kounter;
} McThreadRng;

typedef struct McSampleRecord {
	char *text;
	unsigned long int reject_one;
	unsigned long int reject_two;
} McSampleRecord;

static void checked_snprintf(char *buffer, size_t size, const char *fmt, ...) {
	va_list args;
	int written;

	va_start(args, fmt);
	written = vsnprintf(buffer, size, fmt, args);
	va_end(args);

	if (written < 0 || (size_t)written >= size) {
		fprintf(stderr, "Formatted path exceeded its destination buffer.\n");
		exit(EXIT_FAILURE);
	}
}

static void ensure_directory(const char *path) {
	if (mkdir(path, 0775) != 0 && errno != EEXIST) {
		fprintf(stderr, "Could not create directory '%s': %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}
}

static void *sampler_xcalloc(size_t count, size_t size, const char *label)
{
	void *ptr = calloc(count, size);
	if (!ptr) {
		fprintf(stderr, "Fatal: unable to allocate %s (%zu x %zu bytes)\n", label, count, size);
		exit(EXIT_FAILURE);
	}
	return ptr;
}

static unsigned int mix_sample_seed(unsigned int seed, unsigned long int sample_index)
{
	unsigned int x = seed ^ (unsigned int)(sample_index + 0x9e3779b9U);
	x ^= x >> 16;
	x *= 0x7feb352dU;
	x ^= x >> 15;
	x *= 0x846ca68bU;
	x ^= x >> 16;
	return x ? x : 1U;
}

static void mc_thread_rng_init(McThreadRng *rng, unsigned int seed)
{
	int i, i1, j1;
	unsigned int m, mm[3];
	unsigned int bit;

	rng->kounter = 0;
	rng->ic = 0;
	m = seed;
	rng->iran1 = m;
	for (i = 0; i < 3; i++) {
		mm[i] = rng->iran1 - KWEYL;
		rng->iran1 = mm[i];
	}

	for (i1 = 0; i1 < 43; i1++) {
		rng->iran[i1] = 0;
		for (j1 = 1; j1 <= 32; j1++) {
			i = j1 + 32 * i1;
			bit = m * mm[i % 3];
			bit = ((bit % 64) < 32) ? 0U : 1U;
			rng->iran[i1] += bit << (j1 - 1);
			mm[i % 3] = (mm[0] * mm[1] * mm[2]) % 179;
			m = (53 * m + 1) % 169;
		}
	}

	rng->iran1 = 0;
	for (i = 1; i <= 32; i++) {
		bit = m * mm[i % 3];
		bit = ((bit % 64) < 32) ? 0U : 1U;
		rng->iran1 += bit << (i - 1);
		mm[i % 3] = (mm[0] * mm[1] * mm[2]) % 179;
		m = (53 * m + 1) % 169;
	}
}

static double mc_thread_rng_real(McThreadRng *rng)
{
	int ir = rng->iran[(rng->kounter + 43 - 22) % 43] - rng->iran[rng->kounter] - rng->ic;
	double r;

	if (ir < 0) {
		rng->iran[rng->kounter] = (ir + MaxInt5) % MaxInt5;
		rng->ic = 1;
	} else {
		rng->iran[rng->kounter] = ir % MaxInt5;
		rng->ic = 0;
	}
	rng->iran1 -= KWEYL;
	r = (rng->iran[rng->kounter] - rng->iran1) * INVERSE_MAXINT;
	rng->kounter = (rng->kounter + 1) % 43;
	return r;
}

static int **sampler_alloc_int2(int rows, int cols, const char *label)
{
	int **arr = sampler_xcalloc((size_t)rows, sizeof(*arr), label);
	int *data = sampler_xcalloc((size_t)rows * (size_t)cols, sizeof(*data), label);
	int i;

	for (i = 0; i < rows; i++) {
		arr[i] = data + (size_t)i * (size_t)cols;
	}
	return arr;
}

static void sampler_free_int2(int **arr)
{
	if (arr) {
		free(arr[0]);
		free(arr);
	}
}

static void reset_thread_built_walks(void)
{
	int i, j;
	for (i = 0; i <= vM * vL / 2 - 1; i++) {
		for (j = 0; j <= 2; j++) {
			built_walks_start[i][j] = -1;
			built_walks_end[i][j] = -1;
		}
		for (j = 0; j <= vM * vL * (totalspan + 1) - 1; j++) {
			built_walks_direcs[i][j] = 0;
		}
	}
	num_built_walks = 0;
}

static void copy_left_endhinge_to_thread_walks(unsigned long int secnum, int nth_endhinge)
{
	int nth_walk, i;

	num_built_walks = Lend_num_walks[secnum][nth_endhinge];
	for (nth_walk = 0; nth_walk <= num_built_walks - 1; nth_walk++) {
		for (i = 0; i <= 2; i++) {
			built_walks_start[nth_walk][i] = Lend_start[secnum][nth_endhinge][i][nth_walk];
			built_walks_end[nth_walk][i] = Lend_end[secnum][nth_endhinge][i][nth_walk];
		}
		for (i = 0; i <= vM * vL; i++) {
			built_walks_direcs[nth_walk][i] = Lend_walks[secnum][nth_endhinge][i][nth_walk];
		}
		for (i = vM * vL + 1; i <= vM * vL * (totalspan + 1) - 1; i++) {
			built_walks_direcs[nth_walk][i] = 0;
		}
	}
}

static char *render_thread_sample_record(void)
{
	int max_directions = vM * vL * (totalspan + 1);
	size_t cap = 64 + (size_t)(max_directions + 1) * 4U;
	char *record = sampler_xcalloc(cap, sizeof(*record), "sample record");
	size_t used = 0;
	int written;
	int i = 0;

	written = snprintf(record + used, cap - used, "%d %d %d\n",
		built_walks_start[0][0], built_walks_start[0][1], built_walks_start[0][2]);
	if (written < 0 || (size_t)written >= cap - used) {
		fprintf(stderr, "Fatal: sample record buffer too small\n");
		exit(EXIT_FAILURE);
	}
	used += (size_t)written;

	while (built_walks_direcs[0][i] != 0) {
		written = snprintf(record + used, cap - used, "%d\n", built_walks_direcs[0][i]);
		if (written < 0 || (size_t)written >= cap - used) {
			fprintf(stderr, "Fatal: sample record buffer too small\n");
			exit(EXIT_FAILURE);
		}
		used += (size_t)written;
		i++;
	}

	written = snprintf(record + used, cap - used, "-111\n");
	if (written < 0 || (size_t)written >= cap - used) {
		fprintf(stderr, "Fatal: sample record buffer too small\n");
		exit(EXIT_FAILURE);
	}
	return record;
}

static int requested_sample_threads(void)
{
	const char *value = getenv("MC_SAMPLE_THREADS");
	int threads = 1;

	if (value && *value) {
		char *endptr;
		long parsed;
		errno = 0;
		parsed = strtol(value, &endptr, 10);
		if (errno != 0 || *endptr != '\0' || parsed < 1 || parsed > INT_MAX) {
			fprintf(stderr, "Error: MC_SAMPLE_THREADS must be a positive integer (received '%s')\n", value);
			exit(EXIT_FAILURE);
		}
		threads = (int)parsed;
	}
#ifdef _OPENMP
	if (threads > omp_get_max_threads()) {
		threads = omp_get_max_threads();
	}
#else
	threads = 1;
#endif
	if (threads > samplesize) {
		threads = samplesize;
	}
	return threads;
}

double generate_evectors() {
	printf("Starting integrated spectral solver...\n");
	double fugacity = 1.0;
	
	/* Uniform sampling uses the transition matrix at fugacity x = 1. The
	   solver writes into the legacy global eigenvector buffers used by the
	   sampling tables below. */
	extern double max_eval_LRvec(double fugacity);
	double calculated_evalue = max_eval_LRvec(fugacity) + 1.0;
	printf("Calculated dominant eigenvalue: %.15f (Expected: %.15f)\n", calculated_evalue, dom_evalue);
	dom_evalue = calculated_evalue;
	
	/* If an archival eigenvector file already exists, report numerical drift.
	   The sampler still uses the freshly computed vector so CLI runs do not
	   depend on stale files. */
	char r_filename[100];
	if (ham_check) {
		checked_snprintf(r_filename, sizeof(r_filename), "data/MonteCarlo/MC_Evectors/R_EvectorHam_TS_L%dM%d.txt", L, M);
	} else {
		checked_snprintf(r_filename, sizeof(r_filename), "data/MonteCarlo/MC_Evectors/R_Evector_TS_L%dM%d.txt", L, M);
	}

	FILE *r_fp = fopen(r_filename, "r");
	if (r_fp != NULL) {
		printf("Verifying calculated eigenvectors against archival data in %s...\n", r_filename);
		double max_diff = 0.0;
		int diff_count = 0;
		for (int i = 1; i <= (int)num_tspans; i++) {
			double file_val;
			if (fscanf(r_fp, "%lf", &file_val) == 1) {
				double diff = fabs(R_Evector_solve[0][i] - file_val);
				if (diff > max_diff) max_diff = diff;
				if (diff > 1e-9) diff_count++;
			}
		}
		printf("Verification complete. Max difference: %e, Significant differences (>1e-9): %d\n", max_diff, diff_count);
		fclose(r_fp);
	} else {
		printf("No archival file found for verification (%s). Proceeding with calculated values.\n", r_filename);
	}

	/* Keep the traditional file output for parity audits and external scripts. */
	ensure_directory("data");
	ensure_directory("data/MonteCarlo");
	ensure_directory("data/MonteCarlo/MC_Evectors");
	char export_fn[100];
	if (ham_check) {
		checked_snprintf(export_fn, sizeof(export_fn), "data/MonteCarlo/MC_Evectors/R_EvectorHam_TS_L%dM%d.txt", L, M);
	} else {
		checked_snprintf(export_fn, sizeof(export_fn), "data/MonteCarlo/MC_Evectors/R_Evector_TS_L%dM%d.txt", L, M);
	}
	
	FILE *export_fp = fopen(export_fn, "w");
	if (export_fp != NULL) {
		for (int i = 1; i <= (int)num_tspans; i++) {
			fprintf(export_fp, "%.15f\n", R_Evector_solve[0][i]);
			}
			fclose(export_fp);
			run_metadata_write(export_fn, "mc_master", "right_eigenvector", mode, L, M, totalspan, seednum, runnum, calculated_evalue);
			printf("Calculated eigenvectors exported to %s\n", export_fn);
	}
	return calculated_evalue;
}

static void sampler_fatal_selection(const char *stage, unsigned long int section, unsigned long int limit, double sum, double prob)
{
	fprintf(stderr, "Fatal: sampler could not select %s transition from section %lu (outsections=%lu, cumulative=%.17g, target=%.17g)\n",
		stage, section, limit, sum, prob);
	exit(EXIT_FAILURE);
}

static int choose_left_conditioned_tspan(unsigned long int section, double acceptance, double prob)
{
	unsigned long int limit = num_outsections[section];
	unsigned long int nth_tspan = 1;
	double sumofprobs;

	if (limit == 0 || acceptance <= 0.0 || !isfinite(acceptance)) {
		sampler_fatal_selection("left-conditioned", section, limit, acceptance, prob);
	}

	/* After accepting a left endhinge, choose the first two-span from the same
	   section with probability proportional to the destination eigenvector. */
	sumofprobs = R_Evector_solve[0][t_nrr[section][nth_tspan]] / acceptance;
	while (sumofprobs < prob && nth_tspan < limit) {
		nth_tspan++;
		sumofprobs = sumofprobs + R_Evector_solve[0][t_nrr[section][nth_tspan]] / acceptance;
	}
	if (sumofprobs < prob && fabs(sumofprobs - prob) > 1e-12) {
		sampler_fatal_selection("left-conditioned", section, limit, sumofprobs, prob);
	}
	return (int)nth_tspan;
}

static int choose_middle_tspan(unsigned long int section, unsigned long int current_tspan, double total_prob, double prob)
{
	unsigned long int limit = num_outsections[section];
	unsigned long int nth_tspan = 1;
	double denom;
	double sumofprobs;

	if (limit == 0) {
		sampler_fatal_selection("middle", section, limit, 0.0, prob);
	}
	denom = R_Evector_solve[0][current_tspan] * dom_evalue;
	if (denom <= 0.0 || !isfinite(denom) || !isfinite(total_prob) || total_prob <= 0.0) {
		sampler_fatal_selection("middle", section, limit, total_prob, prob);
	}

	/* Interior spans are conditioned on the current transition and normalized
	   by lambda. total_prob is kept for diagnostics because it should be close
	   to one in a correctly normalized row. */
	sumofprobs = R_Evector_solve[0][t_nrr[section][nth_tspan]] / denom;
	while (sumofprobs < prob && nth_tspan < limit) {
		nth_tspan++;
		sumofprobs = sumofprobs + R_Evector_solve[0][t_nrr[section][nth_tspan]] / denom;
	}
	if (sumofprobs < prob && fabs(sumofprobs - prob) > 1e-12) {
		sampler_fatal_selection("middle", section, limit, sumofprobs, prob);
	}
	return (int)nth_tspan;
}

static McSampleRecord build_parallel_sample_record(const McSamplerWeights *weights, unsigned long int sample_index)
{
	McThreadRng rng;
	McSampleRecord record;
	unsigned long int secnum, sec1, sec2;
	int nth_endhinge;
	int curspan;
	int nth_tspan;
	unsigned long int chosenLEH;
	unsigned long int curLEH;
	int chosenREH;
	double prob;
	double sumofprobs;
	unsigned long int cur_tspan_num;

	memset(&record, 0, sizeof(record));
	mc_thread_rng_init(&rng, mix_sample_seed(seednum, sample_index));

	for (;;) {
		reset_thread_built_walks();

		chosenLEH = floor(mc_thread_rng_real(&rng) * tot_left_endhinges);
		if (chosenLEH >= weights->left_count) chosenLEH = weights->left_count - 1;

		secnum = weights->left_section[chosenLEH];
		nth_endhinge = weights->left_endhinge[chosenLEH];
		curLEH = chosenLEH;
		if (mc_thread_rng_real(&rng) < weights->left_acceptance[curLEH] / weights->max_left_acceptance) {
			copy_left_endhinge_to_thread_walks(secnum, nth_endhinge);
			curspan = 1;
			sec1 = secnum;

			prob = mc_thread_rng_real(&rng);
			nth_tspan = choose_left_conditioned_tspan(sec1, weights->left_acceptance[curLEH], prob);
			add_to_built_walks(sec1, nth_tspan);
			sec2 = t_outsection[sec1][nth_tspan];
			cur_tspan_num = t_nrr[sec1][nth_tspan];
			curspan++;

			while (curspan < totalspan) {
				int i;
				sumofprobs = 0.0;
				for (i = 1; i <= (int)num_outsections[sec2]; i++) {
					sumofprobs = sumofprobs + R_Evector_solve[0][t_nrr[sec2][i]] / R_Evector_solve[0][cur_tspan_num] / dom_evalue;
				}

				if (sumofprobs < 0.9 || sumofprobs > 1.1) {
					fprintf(stderr, "Fatal: middle transition probabilities sum to %f from tspan %lu on sample %lu\n",
						sumofprobs, cur_tspan_num, sample_index);
					exit(EXIT_FAILURE);
				}

				prob = mc_thread_rng_real(&rng) * sumofprobs;
				nth_tspan = choose_middle_tspan(sec2, cur_tspan_num, sumofprobs, prob);
				sec1 = sec2;
				add_to_built_walks(sec1, nth_tspan);
				sec2 = t_outsection[sec1][nth_tspan];
				cur_tspan_num = t_nrr[sec1][nth_tspan];
				curspan++;
			}

			if (mc_thread_rng_real(&rng) < weights->right_acceptance[cur_tspan_num] / weights->max_right_acceptance) {
				if (num_right_endhinges[sec2] == 0) {
					fprintf(stderr, "Fatal: no right endhinges available for sampled terminal section %lu\n", sec2);
					exit(EXIT_FAILURE);
				}
				chosenREH = 1 + floor(mc_thread_rng_real(&rng) * num_right_endhinges[sec2]);
				if ((unsigned long int)chosenREH > num_right_endhinges[sec2]) chosenREH = (int)num_right_endhinges[sec2];
				add_right_endhinge(sec2, chosenREH);
				record.text = render_thread_sample_record();
				return record;
			}
			record.reject_two++;
		} else {
			record.reject_one++;
		}
	}
}

static void rotate_sample_file_if_needed(void)
{
	if (filetotal < (unsigned long int)maxpolys) {
		return;
	}

	printf("Finished filling file %lu. It contains %lu polys. Creating a new one.\n", filenum, filetotal);
	fprintf(fp, "-999\n");
	fclose(fp);
	filenum++;
	{
		const char *file_prefix = ham_check ? "MCpolysHam" : "MCpolys";
		checked_snprintf(filename, sizeof(filename), "%s/%sL%dM%dspan%drun%dnum%lu.txt", output_dir, file_prefix, L, M, totalspan, runnum, filenum);
	}
	fp = fopen(filename, "w");
	if (fp != NULL) {
		fprintf(fp, "UofS\n");
		run_metadata_write(filename, "mc_master", "samples_uofs", mode, L, M, totalspan, seednum, runnum, dom_evalue);
		printf("printed UofS in file '%s'\n", filename);
	} else {
		fprintf(stderr, "Could not open '%s' for writing: %s\n", filename, strerror(errno));
		exit(EXIT_FAILURE);
	}
	filetotal = 0;
}

static void write_sample_record_text(const char *text)
{
	rotate_sample_file_if_needed();
	fputs(text, fp);
	filetotal++;
}

static int parse_int_arg(const char *value, const char *label)
{
	char *endptr;
	long parsed;

	errno = 0;
	parsed = strtol(value, &endptr, 10);
	if (errno != 0 || *endptr != '\0' || parsed < INT_MIN || parsed > INT_MAX) {
		fprintf(stderr, "Error: %s must be an integer (received '%s')\n", label, value);
		exit(EXIT_FAILURE);
	}
	return (int)parsed;
}

static unsigned int parse_uint_arg(const char *value, const char *label)
{
	char *endptr;
	unsigned long parsed;

	errno = 0;
	parsed = strtoul(value, &endptr, 10);
	if (errno != 0 || *endptr != '\0' || parsed > UINT_MAX) {
		fprintf(stderr, "Error: %s must be an unsigned integer (received '%s')\n", label, value);
		exit(EXIT_FAILURE);
	}
	return (unsigned int)parsed;
}

void parse_args(int argc, char *argv[]) {
	int opt;
	while ((opt = getopt(argc, argv, "L:M:s:n:r:S:m:")) != -1) {
		switch (opt) {
			case 'L': L = parse_int_arg(optarg, "Lattice Width (-L)"); break;
			case 'M': M = parse_int_arg(optarg, "Lattice Height (-M)"); break;
			case 's': totalspan = parse_int_arg(optarg, "Span (-s)"); break;
			case 'n': samplesize = parse_int_arg(optarg, "Sample Size (-n)"); break;
			case 'r': runnum = parse_int_arg(optarg, "Run Number (-r)"); break;
			case 'S': seednum = parse_uint_arg(optarg, "Seed (-S)"); break;
			case 'm': 
				mode = parse_int_arg(optarg, "Mode (-m)");
				if (mode == 1 || mode == 3) ham_check = 1;
				break;
			default:
				fprintf(stderr, "Usage: %s [-L L] [-M M] [-s totalspan] [-n samplesize] [-r runnum] [-S seednum] [-m mode]\n", argv[0]);
				fprintf(stderr, "  -m mode      0=SAP, 1=HamSAP, 2=2SAP, 3=Ham2SAP\n");
				fprintf(stderr, "  -L L         Lattice Width\n");
				fprintf(stderr, "  -M M         Lattice Height\n");
				fprintf(stderr, "  -s span      Total Span (number of sections)\n");
				fprintf(stderr, "  -n size      Sample size (number of objects)\n");
				exit(EXIT_FAILURE);
		}
	}
}

void build_transfer_matrix_graph() {
	int ordNum[3], side = 0;
	int i, j;

	clock();

	/* ordNum tracks the current recursive numbering:
	   [left boundary order, right boundary order, hinge-edge count]. */
	ordNum[0] = 1;
	ordNum[1] = 1;
	ordNum[2] = 0;

	for (i=1; i<=max_sections; i++){
		current_hinge_span[i] = newhinge(1);		/* free up space of hinge structure */
		first_hinge_span[i] = current_hinge_span[i];	/* make first_hinge_span[i] point to this first one */
		currentendhinge[i] = newendhinge(1);		/* free up space of endhinge structure */
		firstendhinge[i] = currentendhinge[i];		/* make firstendhinge[i] point to this first one */
		sectionkey[i] = 0;
	}
	printf("space created for hinge structures, endhinge structures, and sectionkey initialized\n");

	//initializations
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			hingestatus[i][j]=0;
			alreadyentered[i][j]=0;
		}
	}

	for(i=0; i<=M-1; i++){
		for(j=0; j<=L; j++){
			colhingeedges[i][j]=0;
		}
	}
	for(i=0; i<=M; i++){
		for(j=0; j<=L-1; j++){
			rowhingeedges[i][j]=0;
		}
	}

	num_walks=0;
	for(i=0; i<=vM*vL-1; i++){
		for(j=0; j<=2; j++){
			curstart[j][i]=-1;
			curend[j][i]=-1;
		}
		for(j=0; j<=vM*vL; j++){
			curwalks[j][i]=0;
		}
	}

	//create 2-spans
	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if( !(i==M && j==L) ){
				printf("2-span: entering at i=%d,j=%d\n", i, j);
				enterhinge(&primary_polygon_state, i, j, side, &ordNum, 0);
			}
			/* Once a start vertex has been tried, later starts skip it so the
			   same two-span is not rediscovered with a different root. */
			alreadyentered[i][j]=1;
		}
	}

	printf("%lu 2-spans recorded\n", num_tspans);
	printf("%lu duplicate 2spans\n", num_duplicate_tspans);

	/* Endhinges are the left and right caps used to close a sampled sequence of
	   two-spans into a full polygon. They have their own discovery pass because
	   they are sampled at the beginning and end of the chain. */
	num_walks=0;
	for(i=0; i<=vM*vL-1; i++){
		for(j=0; j<=2; j++){
			curstart[j][i]=-1;
			curend[j][i]=-1;
		}
		for(j=0; j<=vM*vL; j++){
			curwalks[j][i]=0;
		}
	}

	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			hingestatus[i][j]=0;
			alreadyentered[i][j]=0;
		}
	}

	for(i=0; i<=M-1; i++){
		for(j=0; j<=L; j++){
			colhingeedges[i][j]=0;
		}
	}
	for(i=0; i<=M; i++){
		for(j=0; j<=L-1; j++){
			rowhingeedges[i][j]=0;
		}
	}

	EndOrdNum[0]=1;
	EndOrdNum[1]=0;

	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if( !(i==M && j==L) ){
				printf("Endhinge: entering at i=%d,j=%d\n", i, j);
				enterendhinge(i, j, 0);
			}
			alreadyentered[i][j]=1;
		}
	}

	printf("%lu right endhinges recorded\n", tot_right_endhinges);
	printf("%lu left endhinges recorded\n", tot_left_endhinges);
	printf("%lu duplicate right endhinges\n", num_duplicate_right_endhinges);
	printf("%lu duplicate left endhinges\n", num_duplicate_left_endhinges);

	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);

	printf("converting hinge_span linked lists to arrays...\n");
	conv_to_array();
	printf("done converting hinge_spans to arrays.\n\n");

	printf("converting endhinge linked lists to arrays...\n");
	conv_endhinges_to_array();
	printf("done converting endhinges to arrays.\n\n");

	unsigned long int totalsecs = 0;
	for (i = 1; i <= (unsigned long int)max_sections; i++) {
		if (num_outsections[i] > 0) {
			totalsecs++;
		}
	}
	printf("Final: Total sections=%lu, total 2spans=%lu\n", totalsecs, num_tspans);
	printf("max_sections was set = %d\n", max_sections);
	printf("max_tspans was set = %d\n", max_tspans);
	printf("There are a total of %lu potential two-span patterns.\n", num_2_spans - 1);
	printf("There are a total of %lu valid two-span patterns.\n", num_tspans);
	printf("Number of duplicate 2-spans generated that were not recorded=%lu\n", dupcounter);
	printf("Dominant Eigenvalue: %.15f\n", dom_evalue);

	if (ham_check) {
		check_reachability();
	}

	printf("Finished generating all 2-spans and endhinges.\n");
}

void run_rejection_sampler() {
	if(totalspan < 2){
		printf("totalspan=%d. This program only works for span>=2. Exitting\n", totalspan);
		exit(1);
	}

	printf("\nNOW SAMPLING: %d samples from L=%d, M=%d, span=%d%s\n", samplesize, L, M, totalspan, (ham_check ? " (Hamiltonian)" : ""));

	const McModePaths *paths = mc_mode_paths(mode);
	const char *file_prefix = paths ? paths->sample_prefix : "MCpolys";
	ensure_directory("data");
	ensure_directory("data/MonteCarlo");
	ensure_directory(output_dir);
	checked_snprintf(filename, sizeof(filename), "%s/%sL%dM%dspan%drun%dnum%lu.txt", output_dir, file_prefix, L, M, totalspan, runnum, filenum);

	unsigned long int secnum, sec1, sec2;
	int nth_endhinge, nth_walk;
	int curspan;
	int nth_tspan;

	unsigned long int chosenLEH;
	unsigned long int curLEH;
	int chosenREH;

	unsigned long int curSample=0;

	double prob;
	double sumofprobs;
	unsigned long int cur_tspan_num;
	McSamplerWeightInput weight_input;
	McSamplerWeights weights;

	/* Build the spectral data after the transition graph exists, then derive
	   the rejection envelopes used throughout this sampling run. */
	generate_evectors();

	double* R_Evector;
	R_Evector = R_Evector_solve[0];

	memset(&weight_input, 0, sizeof(weight_input));
	weight_input.max_sections = max_sections;
	weight_input.max_tspans = max_tspans;
	weight_input.total_left_endhinges = tot_left_endhinges;
	weight_input.num_left_endhinges = num_left_endhinges;
	weight_input.num_right_endhinges = num_right_endhinges;
	weight_input.num_outsections = num_outsections;
	weight_input.t_outsection = t_outsection;
	weight_input.t_nrr = t_nrr;
	weight_input.right_evector = R_Evector;
	mc_sampler_weights_build(&weight_input, &weights);

	fp = fopen(filename, "w");
	if(fp != NULL){
		fprintf(fp, "UofS\n");
		run_metadata_write(filename, "mc_master", "samples_uofs", mode, L, M, totalspan, seednum, runnum, dom_evalue);
		printf("printed UofS in file '%s'\n", filename);
	} else {
		fprintf(stderr, "Could not open '%s' for writing: %s\n", filename, strerror(errno));
		exit(1);
	}

	int i, j;
	unsigned long int reject_one=0;
	unsigned long int reject_two=0;
	int sample_threads = requested_sample_threads();

	if (sample_threads > 1) {
#ifdef _OPENMP
		McSampleRecord *records = sampler_xcalloc((size_t)samplesize, sizeof(*records), "parallel sample records");
		int **global_built_walks_start = built_walks_start;
		int **global_built_walks_end = built_walks_end;
		int **global_built_walks_direcs = built_walks_direcs;
		int global_num_built_walks = num_built_walks;
		unsigned long int sample_idx;

		printf("Using %d OpenMP sample worker threads. Set MC_SAMPLE_THREADS=1 for legacy sequential RNG output.\n", sample_threads);
		#pragma omp parallel num_threads(sample_threads) reduction(+:reject_one, reject_two)
		{
			int **thread_built_walks_start = sampler_alloc_int2(vM * vL / 2 + 1, 3, "thread built walk starts");
			int **thread_built_walks_end = sampler_alloc_int2(vM * vL / 2 + 1, 3, "thread built walk ends");
			int **thread_built_walks_direcs = sampler_alloc_int2(vM * vL / 2 + 1, vM * vL * (totalspan + 1) + 1, "thread built walk directions");

			built_walks_start = thread_built_walks_start;
			built_walks_end = thread_built_walks_end;
			built_walks_direcs = thread_built_walks_direcs;
			num_built_walks = 0;

			#pragma omp for schedule(dynamic)
			for (sample_idx = 0; sample_idx < (unsigned long int)samplesize; sample_idx++) {
				McSampleRecord record = build_parallel_sample_record(&weights, sample_idx);
				records[sample_idx] = record;
				reject_one += record.reject_one;
				reject_two += record.reject_two;
			}

			sampler_free_int2(thread_built_walks_start);
			sampler_free_int2(thread_built_walks_end);
			sampler_free_int2(thread_built_walks_direcs);
		}

		built_walks_start = global_built_walks_start;
		built_walks_end = global_built_walks_end;
		built_walks_direcs = global_built_walks_direcs;
		num_built_walks = global_num_built_walks;

		for (sample_idx = 0; sample_idx < (unsigned long int)samplesize; sample_idx++) {
			write_sample_record_text(records[sample_idx].text);
			free(records[sample_idx].text);
		}
		curSample = (unsigned long int)samplesize;
		free(records);

		fprintf(fp, "-999\n");
		fclose(fp);
		printf("File %lu contains %lu polys.\n", filenum, filetotal);
		printf("\nSampling Complete.\n");
		printf("\nRESULTS: %lu samples created from L=%d, M=%d, span=%d\n", curSample, L, M, totalspan);
		printf("Rejected at first step %lu times\n", reject_one);
		printf("Rejected at last step %lu times\n", reject_two);
		mc_sampler_weights_free(&weights);
		return;
#else
		fprintf(stderr, "Warning: MC_SAMPLE_THREADS ignored because this binary was built without OpenMP.\n");
#endif
	}

	/* Each accepted polygon is assembled as:
	   left endhinge -> first two-span -> optional interior two-spans ->
	   right endhinge. Rejections at the first/last step preserve the target
	   measure while avoiding per-sample normalization over all endhinges. */
	while(curSample <= samplesize-1){
		for(i=0; i<=vM*vL/2-1; i++){
			for(j=0; j<=2; j++){
				built_walks_start[i][j]=-1;
				built_walks_end[i][j]=-1;
			}
			for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
				built_walks_direcs[i][j]=0;
			}
		}

		chosenLEH = floor(ran1real_()*tot_left_endhinges);
		if (chosenLEH >= weights.left_count) chosenLEH = weights.left_count - 1;

		secnum = weights.left_section[chosenLEH];
		nth_endhinge = weights.left_endhinge[chosenLEH];
		curLEH = chosenLEH;
		/* Left rejection uses a precomputed upper bound over all left endhinges. */
		if( ran1real_() < weights.left_acceptance[curLEH] / weights.max_left_acceptance ){
			num_built_walks = Lend_num_walks[secnum][nth_endhinge];
			for(nth_walk=0; nth_walk<=num_built_walks-1; nth_walk++){
				for(i=0; i<=2; i++){
					built_walks_start[nth_walk][i] = Lend_start[secnum][nth_endhinge][i][nth_walk];
					built_walks_end[nth_walk][i] = Lend_end[secnum][nth_endhinge][i][nth_walk];
				}
				for(i=0; i<=vM*vL; i++){
					built_walks_direcs[nth_walk][i] = Lend_walks[secnum][nth_endhinge][i][nth_walk];
				}
				for(i=vM*vL+1; i<=vM*vL*(totalspan+1)-1; i++){
					built_walks_direcs[nth_walk][i]=0;
				}
			}
			curspan=1;
			sec1=secnum;

			prob = ran1real_();
			nth_tspan = choose_left_conditioned_tspan(sec1, weights.left_acceptance[curLEH], prob);
			add_to_built_walks(sec1, nth_tspan);
			sec2 = t_outsection[sec1][nth_tspan];
			cur_tspan_num = t_nrr[sec1][nth_tspan];
			curspan++;

			while(curspan < totalspan){
				sumofprobs=0.0;
				for(i=1; i<=num_outsections[sec2]; i++){
					sumofprobs = sumofprobs + R_Evector[t_nrr[sec2][i]] / R_Evector[cur_tspan_num] / dom_evalue;
				}

				/* The row probability should be one; keep the legacy wide
				   tolerance because the archival code used the same guard. */
				if(sumofprobs<0.9 || sumofprobs>1.1){
					printf("second sumofprobs=%f\n", sumofprobs);
					printf("problem with second sumofprobs when leaving tspan number %lu, (made up of %lu and %lu)\n", cur_tspan_num, sec1, sec2);
					printf("This was on sample %lu. Exitting\n", curSample);
					exit(1);
				}

				prob = ran1real_();
				prob = prob * sumofprobs;
				nth_tspan = choose_middle_tspan(sec2, cur_tspan_num, sumofprobs, prob);
				sec1=sec2;
				add_to_built_walks(sec1, nth_tspan);
				sec2 = t_outsection[sec1][nth_tspan];
				cur_tspan_num = t_nrr[sec1][nth_tspan];
				curspan++;
			}

			/* Final rejection accounts for the number of right endhinges
			   available at the terminal section. */
			if( ran1real_() < weights.right_acceptance[cur_tspan_num] / weights.max_right_acceptance ){
				if (num_right_endhinges[sec2] == 0) {
					fprintf(stderr, "Fatal: no right endhinges available for sampled terminal section %lu\n", sec2);
					exit(EXIT_FAILURE);
				}
				chosenREH = 1 + floor(ran1real_()*num_right_endhinges[sec2]);
				if ((unsigned long int)chosenREH > num_right_endhinges[sec2]) chosenREH = (int)num_right_endhinges[sec2];
				add_right_endhinge(sec2, chosenREH);
				curSample++;
				printtofile();
			}
			else{
				reject_two++;
			}
		}
		else{
			reject_one++;
		}
	}

	fprintf(fp, "-999\n");
	fclose(fp);
	printf("File %lu contains %lu polys.\n", filenum, filetotal);

	printf("\nSampling Complete.\n");

	printf("\nRESULTS: %lu samples created from L=%d, M=%d, span=%d\n", curSample, L, M, totalspan);
	printf("Rejected at first step %lu times\n", reject_one);
	printf("Rejected at last step %lu times\n", reject_two);
	
	mc_sampler_weights_free(&weights);
}

void free_sampler_memory() {
	/* The process normally exits after one run, but freeing the largest arrays
	   keeps leak checks useful while refactoring. */
	int i;
	if (built_walks_start) {
		free(built_walks_start[0]);
		free(built_walks_start);
	}
	if (built_walks_end) {
		free(built_walks_end[0]);
		free(built_walks_end);
	}
	if (built_walks_direcs) {
		free(built_walks_direcs[0]);
		free(built_walks_direcs);
	}
	
	for(i = 0; i <= 1; i++) {
		if (R_Evector_solve[i]) free(R_Evector_solve[i]);
		if (L_Evector[i]) free(L_Evector[i]);
	}
}

int main(int argc, char *argv[]) {
	parse_args(argc, argv);
	run_metadata_set_command(argc, argv);

	if (L < 0 || M < 0) {
		fprintf(stderr, "Error: L and M must be non-negative (received L=%d, M=%d).\n", L, M);
		exit(EXIT_FAILURE);
	}
	if (L == 0 && M == 0) {
		fprintf(stderr, "Error: At least one lattice dimension must be positive.\n");
		exit(EXIT_FAILURE);
	}
	if (samplesize <= 0) {
		fprintf(stderr, "Error: Sample Size (-n) must be positive (received %d).\n", samplesize);
		exit(EXIT_FAILURE);
	}
	if (totalspan < 2) {
		fprintf(stderr, "Error: Span (-s) must be at least 2 (received %d).\n", totalspan);
		exit(EXIT_FAILURE);
	}
	if (mode < 0 || mode > 3) {
		fprintf(stderr, "Error: Invalid Simulation Mode (-m %d). Valid modes are 0, 1, 2, 3.\n", mode);
		exit(EXIT_FAILURE);
	}

	if (ham_check && ((L + 1) * (M + 1) * (totalspan + 1)) % 2 != 0) {
		fprintf(stderr, "\nError: Impossible to construct a Hamiltonian polygon with an odd number of vertices on a bipartite lattice (L=%d, M=%d, span=%d).\n\n", L, M, totalspan);
		exit(EXIT_FAILURE);
	}

	if (mode == 2 || mode == 3) {
		/* 2SAP modes have their own integrated samplers because they carry two
		   polygon states through the discovery and sampling process. */
		optind = 1;
		return mode == 3
			? run_integrated_2sap_ham_sampler(argc, argv)
			: run_integrated_2sap_sampler(argc, argv);
	}

	set_system_params();
	allocate_globals();

	const McModePaths *paths = mc_mode_paths(mode);
	output_dir = (char *)(paths ? paths->sample_dir : "data/MonteCarlo/SAPs");

	unsigned int seed=seednum;
	initran_(&seed);

	if(system("clear")){}
	if(system("date")){}

	printf("section-section version\n");
	printf("\nThe variables are as follows M: %i L: %i\n", M, L);

	build_transfer_matrix_graph();
	run_rejection_sampler();

	free_sampler_memory();

	if(system("date")){}
	printf("\nProgram Complete.\n\n");

	return 0;
}
