#include "mc_globals.h"
#include "mc_runtime.h"

#include <errno.h>
#include <limits.h>
#include <stdarg.h>

#define CREATOR_DEFAULT_LIMIT 100000UL

typedef struct CreatorOptions {
	unsigned long int limit;
	int force;
} CreatorOptions;

typedef struct CreatorContext {
	FILE *fp;
	unsigned long int total;
	unsigned long int *sequence_sections;
	int *sequence_nths;
} CreatorContext;

static void checked_snprintf(char *buffer, size_t size, const char *fmt, ...)
{
	va_list args;
	int written;

	va_start(args, fmt);
	written = vsnprintf(buffer, size, fmt, args);
	va_end(args);

	if (written < 0 || (size_t)written >= size) {
		fprintf(stderr, "Fatal: formatted path exceeded its destination buffer\n");
		exit(EXIT_FAILURE);
	}
}

static void ensure_directory(const char *path)
{
	if (mkdir(path, 0775) != 0 && errno != EEXIST) {
		fprintf(stderr, "Fatal: could not create directory '%s': %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}
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

static unsigned long int parse_ulong_arg(const char *value, const char *label)
{
	char *endptr;
	unsigned long parsed;

	errno = 0;
	parsed = strtoul(value, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
		fprintf(stderr, "Error: %s must be an unsigned integer (received '%s')\n", label, value);
		exit(EXIT_FAILURE);
	}
	return parsed;
}

static void print_usage(const char *prog)
{
	printf("Usage: %s -L L -M M -m mode -s span [-C max_count] [-f]\n", prog);
	printf("  -m mode      0=SAP, 1=HamSAP, 2=2SAP, 3=Ham2SAP\n");
	printf("  -C max_count fail if enumeration exceeds this count (default %lu)\n", CREATOR_DEFAULT_LIMIT);
	printf("  -f           force output even if count exceeds -C\n");
}

static void parse_args(int argc, char **argv, CreatorOptions *options)
{
	int opt;

	options->limit = CREATOR_DEFAULT_LIMIT;
	options->force = 0;
	opterr = 0;
	while ((opt = getopt(argc, argv, "L:M:m:s:C:fh")) != -1) {
		switch (opt) {
			case 'L': L = parse_int_arg(optarg, "Lattice Width (-L)"); break;
			case 'M': M = parse_int_arg(optarg, "Lattice Height (-M)"); break;
			case 'm':
				mode = parse_int_arg(optarg, "Mode (-m)");
				ham_check = (mode == 1 || mode == 3);
				break;
			case 's': totalspan = parse_int_arg(optarg, "Span (-s)"); break;
			case 'C': options->limit = parse_ulong_arg(optarg, "Max Count (-C)"); break;
			case 'f': options->force = 1; break;
			case 'h': print_usage(argv[0]); exit(EXIT_SUCCESS);
			default:
				print_usage(argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	if (L < 0 || M < 0) {
		fprintf(stderr, "Error: L and M must be non-negative (received L=%d, M=%d).\n", L, M);
		exit(EXIT_FAILURE);
	}
	if (L == 0 && M == 0) {
		fprintf(stderr, "Error: At least one lattice dimension must be positive.\n");
		exit(EXIT_FAILURE);
	}
	if (mode < 0 || mode > 3) {
		fprintf(stderr, "Error: Invalid mode %d. Valid modes are 0, 1, 2, 3.\n", mode);
		exit(EXIT_FAILURE);
	}
	if (totalspan < 2) {
		fprintf(stderr, "Error: Span (-s) must be at least 2 (received %d).\n", totalspan);
		exit(EXIT_FAILURE);
	}
}

static void reset_built_walks(void)
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
}

static void load_left_endhinge(unsigned long int section, int nth_endhinge)
{
	int nth_walk, i;

	num_built_walks = Lend_num_walks[section][nth_endhinge];
	for (nth_walk = 0; nth_walk <= num_built_walks - 1; nth_walk++) {
		for (i = 0; i <= 2; i++) {
			built_walks_start[nth_walk][i] = Lend_start[section][nth_endhinge][i][nth_walk];
			built_walks_end[nth_walk][i] = Lend_end[section][nth_endhinge][i][nth_walk];
		}
		for (i = 0; i <= vM * vL; i++) {
			built_walks_direcs[nth_walk][i] = Lend_walks[section][nth_endhinge][i][nth_walk];
		}
		for (i = vM * vL + 1; i <= vM * vL * (totalspan + 1) - 1; i++) {
			built_walks_direcs[nth_walk][i] = 0;
		}
	}
}

static void print_current_sap(FILE *out)
{
	int i = 0;

	fprintf(out, "%d %d %d\n", built_walks_start[0][0], built_walks_start[0][1], built_walks_start[0][2]);
	while (built_walks_direcs[0][i] != 0) {
		fprintf(out, "%d\n", built_walks_direcs[0][i]);
		i++;
	}
	fprintf(out, "-111\n");
}

static void build_and_maybe_print_sap(
	CreatorContext *ctx,
	unsigned long int left_section,
	int left_endhinge,
	int right_endhinge,
	int sequence_len)
{
	int idx;
	unsigned long int terminal_section;

	reset_built_walks();
	load_left_endhinge(left_section, left_endhinge);
	for (idx = 0; idx < sequence_len; idx++) {
		add_to_built_walks(ctx->sequence_sections[idx], ctx->sequence_nths[idx]);
	}
	terminal_section = t_outsection[ctx->sequence_sections[sequence_len - 1]][ctx->sequence_nths[sequence_len - 1]];
	add_right_endhinge(terminal_section, right_endhinge);

	if (ctx->fp) {
		print_current_sap(ctx->fp);
	}
	ctx->total++;
}

static void enumerate_sap_paths(
	CreatorContext *ctx,
	unsigned long int left_section,
	int left_endhinge,
	unsigned long int current_section,
	int depth,
	int target_depth)
{
	unsigned long int nth;

	if (depth == target_depth) {
		for (nth = 1; nth <= num_right_endhinges[current_section]; nth++) {
			build_and_maybe_print_sap(ctx, left_section, left_endhinge, (int)nth, target_depth);
		}
		return;
	}

	for (nth = 1; nth <= num_outsections[current_section]; nth++) {
		ctx->sequence_sections[depth] = current_section;
		ctx->sequence_nths[depth] = (int)nth;
		enumerate_sap_paths(ctx, left_section, left_endhinge, t_outsection[current_section][nth], depth + 1, target_depth);
	}
}

static unsigned long int enumerate_saps(FILE *out)
{
	CreatorContext ctx;
	unsigned long int section;
	int left;

	memset(&ctx, 0, sizeof(ctx));
	ctx.fp = out;
	ctx.sequence_sections = mc_xcalloc((size_t)totalspan, sizeof(*ctx.sequence_sections), "creator section sequence");
	ctx.sequence_nths = mc_xcalloc((size_t)totalspan, sizeof(*ctx.sequence_nths), "creator transition sequence");

	for (section = 1; section <= (unsigned long int)max_sections; section++) {
		for (left = 1; left <= num_left_endhinges[section]; left++) {
			enumerate_sap_paths(&ctx, section, left, section, 0, totalspan - 1);
		}
	}

	free(ctx.sequence_sections);
	free(ctx.sequence_nths);
	return ctx.total;
}

static void build_transfer_matrix_graph(void)
{
	int ordNum[3], side = 0;
	int i, j;
	unsigned long int section;

	ordNum[0] = 1;
	ordNum[1] = 1;
	ordNum[2] = 0;

	for (i = 1; i <= max_sections; i++) {
		current_hinge_span[i] = newhinge(1);
		first_hinge_span[i] = current_hinge_span[i];
		currentendhinge[i] = newendhinge(1);
		firstendhinge[i] = currentendhinge[i];
		sectionkey[i] = 0;
	}

	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			hingestatus[i][j] = 0;
			alreadyentered[i][j] = 0;
		}
	}
	for (i = 0; i <= M - 1; i++) for (j = 0; j <= L; j++) colhingeedges[i][j] = 0;
	for (i = 0; i <= M; i++) for (j = 0; j <= L - 1; j++) rowhingeedges[i][j] = 0;

	num_walks = 0;
	for (i = 0; i <= vM * vL - 1; i++) {
		for (j = 0; j <= 2; j++) {
			curstart[j][i] = -1;
			curend[j][i] = -1;
		}
		for (j = 0; j <= vM * vL; j++) {
			curwalks[j][i] = 0;
		}
	}

	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if (!(i == M && j == L)) {
				enterhinge(&primary_polygon_state, i, j, side, &ordNum, 0);
			}
			alreadyentered[i][j] = 1;
		}
	}

	num_walks = 0;
	for (i = 0; i <= vM * vL - 1; i++) {
		for (j = 0; j <= 2; j++) {
			curstart[j][i] = -1;
			curend[j][i] = -1;
		}
		for (j = 0; j <= vM * vL; j++) {
			curwalks[j][i] = 0;
		}
	}
	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			hingestatus[i][j] = 0;
			alreadyentered[i][j] = 0;
		}
	}
	for (i = 0; i <= M - 1; i++) for (j = 0; j <= L; j++) colhingeedges[i][j] = 0;
	for (i = 0; i <= M; i++) for (j = 0; j <= L - 1; j++) rowhingeedges[i][j] = 0;

	EndOrdNum[0] = 1;
	EndOrdNum[1] = 0;
	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if (!(i == M && j == L)) {
				enterendhinge(i, j, 0);
			}
			alreadyentered[i][j] = 1;
		}
	}

	conv_to_array();
	conv_endhinges_to_array();

	for (section = 1; section <= (unsigned long int)max_sections; section++) {
		if (num_left_endhinges[section] || num_right_endhinges[section] || num_outsections[section]) {
			continue;
		}
	}
}

static const char *mode_dir_name(void)
{
	switch (mode) {
		case 0: return "SAPs";
		case 1: return "HamSAPs";
		case 2: return "2SAPs";
		case 3: return "Ham2SAPs";
		default: return "Unknown";
	}
}

static int run_sap_creator(const CreatorOptions *options)
{
	char outdir[160];
	char outfile[240];
	FILE *out;
	unsigned long int count;

	ham_check = (mode == 1);
	set_system_params();
	allocate_globals();
	build_transfer_matrix_graph();

	count = enumerate_saps(NULL);
	if (!options->force && count > options->limit) {
		fprintf(stderr,
			"Refusing to write %lu polygons; limit is %lu. Re-run with -C <limit> or -f to force.\n",
			count, options->limit);
		return EXIT_FAILURE;
	}

	ensure_directory("data");
	ensure_directory("data/CreatorAll");
	checked_snprintf(outdir, sizeof(outdir), "data/CreatorAll/%s", mode_dir_name());
	ensure_directory(outdir);
	checked_snprintf(outfile, sizeof(outfile), "%s/All%sL%dM%dspan%d.txt", outdir, mode_dir_name(), L, M, totalspan);

	out = fopen(outfile, "w");
	if (!out) {
		fprintf(stderr, "Fatal: could not open '%s': %s\n", outfile, strerror(errno));
		return EXIT_FAILURE;
	}
	fprintf(out, "UofS\n");
	(void)enumerate_saps(out);
	fprintf(out, "-999\n");
	fclose(out);

	printf("Generated %lu %s of exact span %d in %dx%d tube.\n", count, mode_dir_name(), totalspan, L, M);
	printf("Wrote %s\n", outfile);
	return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
	CreatorOptions options;

	parse_args(argc, argv, &options);
	if (mode == 2 || mode == 3) {
		optind = 1;
		return mode == 3
			? run_integrated_2sap_ham_creator_all(argc, argv, options.limit, options.force)
			: run_integrated_2sap_creator_all(argc, argv, options.limit, options.force);
	}
	return run_sap_creator(&options);
}
