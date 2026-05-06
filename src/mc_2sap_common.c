#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/stat.h>

#include "mc_2sap_common.h"
#include "mc_runtime.h"
#include "mc_spectral.h"
#include "run_metadata.h"

extern int L;
extern int M;
extern int max_sections;
extern int max_tspans;
extern unsigned long int max_keynum;
extern double dom_evalue;
extern int errno;

static void ensure_directory(const char *path)
{
	if (mkdir(path, 0775) != 0 && errno != EEXIST) {
		fprintf(stderr, "Fatal: could not create directory '%s': %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}
}

typedef struct Mc2SapSystemPreset {
	int hamiltonian;
	int M;
	int L;
	int max_sections;
	unsigned long int max_keynum;
	int max_tspans;
	double dom_evalue;
} Mc2SapSystemPreset;

static const Mc2SapModeSpec MC_2SAP_MODE_SPECS[] = {
	{0, "data/MonteCarlo/2SAPs", "MC2SAPs", "data/CreatorAll/All_2SAPs", "All2SAPs", "2SAPs"},
	{1, "data/MonteCarlo/Ham2SAPs", "MC2SAPsHam", "data/CreatorAll/All_Ham2SAPs", "AllHam2SAPs", "Ham2SAPs"},
};

static const Mc2SapSystemPreset MC_2SAP_SYSTEM_PRESETS[] = {
	{0, 1, 2, 73, 152, 1048, 9.455960990693537},
	{0, 1, 3, 742, 4048, 108410, 65.012508345749453},
	{0, 1, 4, 9309, 108386, 9419930, 380.824254675661734},
	{0, 2, 2, 2619, 21546, 1126578, 180.511702878389116},
	{1, 1, 2, 73, 152, 490, 5.534148126030995},
	{1, 1, 3, 742, 4048, 37454, 24.378235811209002},
	{1, 1, 4, 9309, 108386, 2598620, 97.848097677835298},
	{1, 2, 2, 2619, 21540, 495792, 53.686533141903084},
};

void mc_2sap_set_system_params(int is_hamiltonian)
{
	size_t i;

	for (i = 0; i < sizeof(MC_2SAP_SYSTEM_PRESETS) / sizeof(MC_2SAP_SYSTEM_PRESETS[0]); i++) {
		const Mc2SapSystemPreset *preset = &MC_2SAP_SYSTEM_PRESETS[i];

		if (preset->hamiltonian == (is_hamiltonian ? 1 : 0) && preset->M == M && preset->L == L) {
			max_sections = preset->max_sections;
			max_keynum = preset->max_keynum;
			max_tspans = preset->max_tspans;
			dom_evalue = preset->dom_evalue;
			return;
		}
	}
	fprintf(stderr, "Unsupported %s2SAP L and M values (%d, %d).\n",
		is_hamiltonian ? "Hamiltonian " : "", L, M);
	exit(EXIT_FAILURE);
}

const Mc2SapModeSpec *mc_2sap_mode_spec(int is_hamiltonian)
{
	return &MC_2SAP_MODE_SPECS[is_hamiltonian ? 1 : 0];
}

void mc_2sap_ignore_system_result(const char *command)
{
	int status = system(command);
	(void)status;
}

int ***mc_2sap_alloc_int3_table(unsigned long entries, int components, const unsigned int *widths, const char *label)
{
	int ***table;
	int **component_ptrs;
	int *data;
	unsigned long i;
	int j;
	size_t ptr_count = 0;
	size_t data_count = 0;
	size_t ptr_index = 0;
	size_t data_index = 0;

	if (components <= 0) {
		fprintf(stderr, "Fatal: invalid component count for %s\n", label);
		exit(EXIT_FAILURE);
	}

	ptr_count = (size_t)entries * (size_t)components;
	for (i = 1; i <= entries; i++) {
		data_count += (size_t)components * (size_t)widths[i];
	}

	table = (int ***)mc_xcalloc((size_t)entries + 1, sizeof(*table), label);
	component_ptrs = (int **)mc_xcalloc(ptr_count ? ptr_count : 1, sizeof(*component_ptrs), label);
	data = (int *)mc_xcalloc(data_count ? data_count : 1, sizeof(*data), label);

	for (i = 1; i <= entries; i++) {
		table[i] = &component_ptrs[ptr_index];
		ptr_index += (size_t)components;
		for (j = 0; j < components; j++) {
			table[i][j] = &data[data_index];
			data_index += (size_t)widths[i];
		}
	}

	return table;
}

int ***mc_2sap_alloc_int3_fixed(unsigned long entries, int components, unsigned int width, const char *label)
{
	int ***table;
	int **component_ptrs;
	int *data;
	unsigned long i;
	int j;
	size_t ptr_count;
	size_t data_count;
	size_t ptr_index = 0;
	size_t data_index = 0;

	if (components <= 0) {
		fprintf(stderr, "Fatal: invalid component count for %s\n", label);
		exit(EXIT_FAILURE);
	}

	ptr_count = (size_t)entries * (size_t)components;
	data_count = ptr_count * (size_t)width;
	table = (int ***)mc_xcalloc((size_t)entries + 1, sizeof(*table), label);
	component_ptrs = (int **)mc_xcalloc(ptr_count ? ptr_count : 1, sizeof(*component_ptrs), label);
	data = (int *)mc_xcalloc(data_count ? data_count : 1, sizeof(*data), label);

	for (i = 1; i <= entries; i++) {
		table[i] = &component_ptrs[ptr_index];
		ptr_index += (size_t)components;
		for (j = 0; j < components; j++) {
			table[i][j] = &data[data_index];
			data_index += (size_t)width;
		}
	}

	return table;
}

void mc_2sap_init_int_rows(int **rows, size_t row_count, size_t width, const char *label)
{
	int *data;
	size_t row;
	size_t safe_width = width ? width : 1;

	data = (int *)mc_xcalloc(row_count * safe_width, sizeof(*data), label);
	for (row = 0; row < row_count; row++) {
		rows[row] = &data[row * safe_width];
	}
}

int **mc_2sap_alloc_int_rows(size_t row_count, size_t width, const char *label)
{
	int **rows = (int **)mc_xcalloc(row_count ? row_count : 1, sizeof(*rows), label);

	mc_2sap_init_int_rows(rows, row_count, width, label);
	return rows;
}

int mc_2sap_reverse_direction(int direction, const char *caller_name)
{
	switch (direction) {
	case 1:
		return 2;
	case 2:
		return 1;
	case 3:
		return 4;
	case 4:
		return 3;
	case 5:
		return 6;
	case 6:
		return 5;
	default:
		fprintf(stderr, "Fatal: invalid direction %d in %s\n", direction, caller_name);
		exit(EXIT_FAILURE);
	}
}

void mc_2sap_open_sample_file(const Mc2SapSampleWriterConfig *config)
{
	ensure_directory("data");
	ensure_directory("data/MonteCarlo");
	ensure_directory(config->outdir);
	mc_checked_snprintf(config->filename, config->filename_size, "%s/%sL%dM%dspan%drun%dnum%lu.txt",
		config->outdir, config->prefix, config->L, config->M, config->totalspan, config->runnum, *config->filenum);
	*config->fp = fopen(config->filename, "w");
	if (*config->fp) {
		fprintf(*config->fp, "UofS\n");
		run_metadata_write(config->filename, "mc_master", "samples_uofs", config->mode,
			config->L, config->M, config->totalspan, config->seednum, config->runnum, config->dom_evalue);
		printf("printed UofS in file '%s'\n", config->filename);
		return;
	}
	fprintf(stderr, "Fatal: could not open '%s': %s\n", config->filename, strerror(errno));
	exit(EXIT_FAILURE);
}

void mc_2sap_write_sample_pair(const Mc2SapSampleWriterConfig *config)
{
	int i;

	if (*config->filetotal >= (unsigned long int)config->maxpolys) {
		printf("Finished filling file %lu. It contains %lu 2SAPs. Creating a new one.\n",
			*config->filenum, *config->filetotal);
		fprintf(*config->fp, "-999\n");
		fclose(*config->fp);
		(*config->filenum)++;
		mc_2sap_open_sample_file(config);
		*config->filetotal = 0;
	}

	fprintf(*config->fp, "%d %d %d\n", config->poly1_start[0], config->poly1_start[1], config->poly1_start[2]);
	for (i = 0; config->poly1_directions[i] != 0; i++) {
		fprintf(*config->fp, "%d\n", config->poly1_directions[i]);
	}
	fprintf(*config->fp, "-111\n");

	fprintf(*config->fp, "%d %d %d\n", config->poly2_start[0], config->poly2_start[1], config->poly2_start[2]);
	for (i = 0; config->poly2_directions[i] != 0; i++) {
		fprintf(*config->fp, "%d\n", config->poly2_directions[i]);
	}
	fprintf(*config->fp, "-111\n");
	(*config->filetotal)++;
}

static size_t mc_2sap_append_format(char *buffer, size_t capacity, size_t offset, const char *format, ...)
{
	int needed;
	va_list args;

	if (offset >= capacity) {
		fprintf(stderr, "Fatal: 2SAP formatted buffer exceeded capacity\n");
		exit(EXIT_FAILURE);
	}
	va_start(args, format);
	needed = vsnprintf(buffer + offset, capacity - offset, format, args);
	va_end(args);
	if (needed < 0 || (size_t)needed >= capacity - offset) {
		fprintf(stderr, "Fatal: 2SAP formatted buffer exceeded capacity\n");
		exit(EXIT_FAILURE);
	}
	return offset + (size_t)needed;
}

char *mc_2sap_walk_string(const int *start, const int *directions)
{
	char *buffer;
	size_t offset = 0;
	int i = 0;

	buffer = (char *)mc_xcalloc(4096, sizeof(*buffer), "2SAP walk key");
	offset = mc_2sap_append_format(buffer, 4096, offset, "%d,%d,%d;",
		start[0], start[1], start[2]);
	while (directions[i] != 0) {
		offset = mc_2sap_append_format(buffer, 4096, offset, "%d,", directions[i]);
		i++;
	}
	return buffer;
}

char *mc_2sap_unordered_pair_key(
	const int *poly1_start,
	const int *poly1_directions,
	const int *poly2_start,
	const int *poly2_directions)
{
	char *a = mc_2sap_walk_string(poly1_start, poly1_directions);
	char *b = mc_2sap_walk_string(poly2_start, poly2_directions);
	char *key = (char *)mc_xcalloc(8192, sizeof(*key), "2SAP unordered pair key");

	if (strcmp(a, b) <= 0) {
		mc_checked_snprintf(key, 8192, "%s|%s", a, b);
	} else {
		mc_checked_snprintf(key, 8192, "%s|%s", b, a);
	}
	free(a);
	free(b);
	return key;
}

void mc_2sap_print_walk(FILE *out, const int *start, const int *directions)
{
	int i = 0;

	fprintf(out, "%d %d %d\n", start[0], start[1], start[2]);
	while (directions[i] != 0) {
		fprintf(out, "%d\n", directions[i]);
		i++;
	}
	fprintf(out, "-111\n");
}

void mc_2sap_print_unordered_pair(
	FILE *out,
	const int *poly1_start,
	const int *poly1_directions,
	const int *poly2_start,
	const int *poly2_directions)
{
	char *a = mc_2sap_walk_string(poly1_start, poly1_directions);
	char *b = mc_2sap_walk_string(poly2_start, poly2_directions);
	int first = strcmp(a, b) <= 0 ? 0 : 1;

	if (first == 0) {
		mc_2sap_print_walk(out, poly1_start, poly1_directions);
		mc_2sap_print_walk(out, poly2_start, poly2_directions);
	} else {
		mc_2sap_print_walk(out, poly2_start, poly2_directions);
		mc_2sap_print_walk(out, poly1_start, poly1_directions);
	}
	free(a);
	free(b);
}

double mc_2sap_max_eval(
	unsigned long int max_keynum_value,
	unsigned long int max_tspans_value,
	unsigned long int *num_outsections,
	unsigned long int **tspans_outsection,
	unsigned long int **tspans_edges,
	unsigned long int **tspans_nrr,
	double *left[2],
	double *right[2],
	double force,
	int lattice_L,
	int lattice_M,
	int hamiltonian,
	double fugacity)
{
	McTransitionSpectralInput input;

	memset(&input, 0, sizeof(input));
	input.max_keynum = max_keynum_value;
	input.max_tspans = max_tspans_value;
	input.num_outsections = num_outsections;
	input.tspans_outsection = tspans_outsection;
	input.tspans_edges = tspans_edges;
	input.tspans_nrr = tspans_nrr;
	input.left[0] = left[0];
	input.left[1] = left[1];
	input.right[0] = right[0];
	input.right[1] = right[1];
	input.force = force;
	input.L = lattice_L;
	input.M = lattice_M;
	input.hamiltonian = hamiltonian;
	return mc_spectral_max_eval_2sap(&input, fugacity);
}

static size_t mc_2sap_next_power_of_two(size_t value)
{
	size_t capacity = 1;
	while (capacity < value) {
		capacity <<= 1;
	}
	return capacity;
}

static size_t mc_2sap_pair_hash(unsigned long int first, unsigned long int second)
{
	unsigned long int hash = 14695981039346656037ULL;

	hash ^= first;
	hash *= 1099511628211ULL;
	hash ^= second;
	hash *= 1099511628211ULL;
	return (size_t)hash;
}

void mc_2sap_pair_map_init(Mc2SapPairMap *map, unsigned long int max_entries)
{
	size_t capacity = mc_2sap_next_power_of_two(((size_t)max_entries + 1) * 2);

	map->capacity = capacity;
	map->first = (unsigned long int *)mc_xcalloc(capacity, sizeof(*map->first), "2SAP pair map first keys");
	map->second = (unsigned long int *)mc_xcalloc(capacity, sizeof(*map->second), "2SAP pair map second keys");
	map->value = (unsigned long int *)mc_xcalloc(capacity, sizeof(*map->value), "2SAP pair map values");
}

void mc_2sap_pair_map_free(Mc2SapPairMap *map)
{
	free(map->first);
	free(map->second);
	free(map->value);
	map->first = NULL;
	map->second = NULL;
	map->value = NULL;
	map->capacity = 0;
}

void mc_2sap_pair_map_put(Mc2SapPairMap *map, unsigned long int first, unsigned long int second, unsigned long int value)
{
	size_t mask = map->capacity - 1;
	size_t slot = mc_2sap_pair_hash(first, second) & mask;

	if (first == 0 || second == 0 || value == 0) {
		fprintf(stderr, "Fatal: invalid zero key while building 2SAP pair map\n");
		exit(EXIT_FAILURE);
	}

	while (map->value[slot] != 0) {
		if (map->first[slot] == first && map->second[slot] == second) {
			map->value[slot] = value;
			return;
		}
		slot = (slot + 1) & mask;
	}

	map->first[slot] = first;
	map->second[slot] = second;
	map->value[slot] = value;
}

unsigned long int mc_2sap_pair_map_get(const Mc2SapPairMap *map, unsigned long int first, unsigned long int second)
{
	size_t mask = map->capacity - 1;
	size_t slot = mc_2sap_pair_hash(first, second) & mask;

	while (map->value[slot] != 0) {
		if (map->first[slot] == first && map->second[slot] == second) {
			return map->value[slot];
		}
		slot = (slot + 1) & mask;
	}

	return 0;
}

void mc_2sap_pair_map_build(Mc2SapPairMap *map, unsigned long int max_entries, unsigned long int pairs[][2])
{
	unsigned long int key;

	mc_2sap_pair_map_init(map, max_entries);
	for (key = 1; key <= max_entries; key++) {
		if (pairs[key][0] == 0) {
			continue;
		}
		mc_2sap_pair_map_put(map, pairs[key][0], pairs[key][1], key);
	}
}

typedef struct Mc2SapCreatorSeenNode {
	char *key;
	struct Mc2SapCreatorSeenNode *next;
} Mc2SapCreatorSeenNode;

typedef struct Mc2SapCreatorContext {
	const Mc2SapCreatorConfig *config;
	FILE *fp;
	unsigned long int total;
	unsigned long int file_total;
	unsigned long int file_num;
	char last_path[512];
	unsigned long int *sequence_sections;
	int *sequence_nths;
	Mc2SapCreatorSeenNode **seen;
	size_t seen_size;
} Mc2SapCreatorContext;

#define MC_2SAP_CREATOR_HASH_SIZE 262144

static unsigned long int mc_2sap_creator_hash_string(const char *text)
{
	unsigned long int h = 14695981039346656037ULL;
	while (*text) {
		h ^= (unsigned char)*text++;
		h *= 1099511628211ULL;
	}
	return h;
}

static char *mc_2sap_creator_strdup(const char *text)
{
	size_t len = strlen(text) + 1;
	char *copy = (char *)mc_xcalloc(len, sizeof(*copy), "creator key copy");
	memcpy(copy, text, len);
	return copy;
}

static void mc_2sap_creator_seen_clear(Mc2SapCreatorContext *ctx)
{
	size_t i;

	if (!ctx->seen) {
		return;
	}
	for (i = 0; i < ctx->seen_size; i++) {
		Mc2SapCreatorSeenNode *node = ctx->seen[i];
		while (node) {
			Mc2SapCreatorSeenNode *next = node->next;
			free(node->key);
			free(node);
			node = next;
		}
		ctx->seen[i] = NULL;
	}
}

static int mc_2sap_creator_seen_insert(Mc2SapCreatorContext *ctx, const char *key)
{
	unsigned long int bucket = mc_2sap_creator_hash_string(key) % ctx->seen_size;
	Mc2SapCreatorSeenNode *node = ctx->seen[bucket];

	while (node) {
		if (strcmp(node->key, key) == 0) {
			return 0;
		}
		node = node->next;
	}
	node = (Mc2SapCreatorSeenNode *)mc_xcalloc(1, sizeof(*node), "creator seen node");
	node->key = mc_2sap_creator_strdup(key);
	node->next = ctx->seen[bucket];
	ctx->seen[bucket] = node;
	return 1;
}

static void mc_2sap_creator_open_next_file(Mc2SapCreatorContext *ctx)
{
	const Mc2SapCreatorConfig *config = ctx->config;

	if (ctx->fp) {
		fprintf(ctx->fp, "-999\n");
		fclose(ctx->fp);
		ctx->fp = NULL;
	}
	mc_checked_snprintf(ctx->last_path, sizeof(ctx->last_path), "%s/%sL%dM%dspan%dnum%lu.txt",
		config->outdir, config->prefix, config->L, config->M, config->totalspan, ctx->file_num);
	ctx->fp = fopen(ctx->last_path, "w");
	if (!ctx->fp) {
		fprintf(stderr, "Fatal: could not open '%s': %s\n", ctx->last_path, strerror(errno));
		exit(EXIT_FAILURE);
	}
	fprintf(ctx->fp, "UofS\n");
	ctx->file_total = 0;
}

static void mc_2sap_creator_prepare_write(Mc2SapCreatorContext *ctx)
{
	ctx->file_num = 1;
	mc_2sap_creator_open_next_file(ctx);
}

static void mc_2sap_creator_finish_write(Mc2SapCreatorContext *ctx)
{
	if (ctx->fp) {
		fprintf(ctx->fp, "-999\n");
		fclose(ctx->fp);
		ctx->fp = NULL;
	}
}

static void mc_2sap_creator_write_pair(Mc2SapCreatorContext *ctx)
{
	if (!ctx->fp) {
		return;
	}
	if (ctx->file_total >= (unsigned long int)ctx->config->maxpolys) {
		ctx->file_num++;
		mc_2sap_creator_open_next_file(ctx);
	}
	ctx->config->print_pair(ctx->config->user, ctx->fp);
	ctx->file_total++;
}

static void mc_2sap_creator_build_leaf(
	Mc2SapCreatorContext *ctx,
	unsigned long int left_section,
	int left_endhinge,
	int right_endhinge,
	int sequence_len)
{
	const Mc2SapCreatorConfig *config = ctx->config;
	int idx;
	unsigned long int terminal_section;
	char *key;

	config->reset_built_walks(config->user);
	config->load_left_endhinge(config->user, left_section, left_endhinge);
	for (idx = 0; idx < sequence_len; idx++) {
		config->add_transition(config->user, ctx->sequence_sections[idx], ctx->sequence_nths[idx]);
	}
	terminal_section = config->t_outsection[ctx->sequence_sections[sequence_len - 1]][ctx->sequence_nths[sequence_len - 1]];
	config->add_right_endhinge(config->user, terminal_section, right_endhinge);

	key = config->pair_key(config->user);
	if (mc_2sap_creator_seen_insert(ctx, key)) {
		mc_2sap_creator_write_pair(ctx);
		ctx->total++;
	}
	free(key);
}

static void mc_2sap_creator_enumerate_paths(
	Mc2SapCreatorContext *ctx,
	unsigned long int left_section,
	int left_endhinge,
	unsigned long int current_section,
	int depth,
	int target_depth)
{
	const Mc2SapCreatorConfig *config = ctx->config;
	unsigned long int nth;

	if (depth == target_depth) {
		for (nth = 1; nth <= config->num_right_endhinges[current_section]; nth++) {
			mc_2sap_creator_build_leaf(ctx, left_section, left_endhinge, (int)nth, target_depth);
		}
		return;
	}
	for (nth = 1; nth <= config->num_outsections[current_section]; nth++) {
		ctx->sequence_sections[depth] = current_section;
		ctx->sequence_nths[depth] = (int)nth;
		mc_2sap_creator_enumerate_paths(ctx, left_section, left_endhinge, config->t_outsection[current_section][nth], depth + 1, target_depth);
	}
}

static unsigned long int mc_2sap_creator_enumerate(Mc2SapCreatorContext *ctx, int write_output)
{
	const Mc2SapCreatorConfig *config = ctx->config;
	unsigned long int section;
	int left;

	ctx->total = 0;
	ctx->sequence_sections = (unsigned long int *)mc_xcalloc((size_t)config->totalspan, sizeof(*ctx->sequence_sections), "creator section sequence");
	ctx->sequence_nths = (int *)mc_xcalloc((size_t)config->totalspan, sizeof(*ctx->sequence_nths), "creator transition sequence");
	ctx->seen = (Mc2SapCreatorSeenNode **)mc_xcalloc(MC_2SAP_CREATOR_HASH_SIZE, sizeof(*ctx->seen), "creator seen table");
	ctx->seen_size = MC_2SAP_CREATOR_HASH_SIZE;
	if (write_output) {
		mc_2sap_creator_prepare_write(ctx);
	}

	for (section = 1; section <= config->max_keynum; section++) {
		for (left = 1; left <= config->num_left_endhinges[section]; left++) {
			mc_2sap_creator_enumerate_paths(ctx, section, left, section, 0, config->totalspan - 1);
		}
	}

	free(ctx->sequence_sections);
	free(ctx->sequence_nths);
	ctx->sequence_sections = NULL;
	ctx->sequence_nths = NULL;
	mc_2sap_creator_seen_clear(ctx);
	free(ctx->seen);
	ctx->seen = NULL;
	ctx->seen_size = 0;
	return ctx->total;
}

int mc_2sap_run_creator_all(const Mc2SapCreatorConfig *config)
{
	Mc2SapCreatorContext ctx;
	unsigned long int count;

	if (config->totalspan < 2) {
		fprintf(stderr, "Error: Span (-s) must be at least 2 (received %d).\n", config->totalspan);
		return EXIT_FAILURE;
	}

	memset(&ctx, 0, sizeof(ctx));
	ctx.config = config;
	count = mc_2sap_creator_enumerate(&ctx, 0);
	if (!config->force && count > config->limit) {
		fprintf(stderr,
			"Refusing to write %lu unordered %s; limit is %lu. Re-run with -C <limit> or -f to force.\n",
			count, config->object_label, config->limit);
		return EXIT_FAILURE;
	}

	ensure_directory("data");
	ensure_directory("data/CreatorAll");
	ensure_directory(config->outdir);
	memset(&ctx, 0, sizeof(ctx));
	ctx.config = config;
	(void)mc_2sap_creator_enumerate(&ctx, 1);
	mc_2sap_creator_finish_write(&ctx);
	{
		char metadata_target[256];
		mc_checked_snprintf(metadata_target, sizeof(metadata_target), "%s/%sL%dM%dspan%d.summary",
			config->outdir, config->prefix, config->L, config->M, config->totalspan);
		run_metadata_write(metadata_target, "creator_all", "creator_all_summary", config->mode,
			config->L, config->M, config->totalspan, config->seednum, config->runnum, 0.0);
	}

	printf("Generated %lu unordered %s of exact span %d in %dx%d tube.\n",
		count, config->object_label, config->totalspan, config->L, config->M);
	printf("Wrote %lu file(s) under %s, with at most %d entries per file.\n",
		ctx.file_num, config->outdir, config->maxpolys);
	return EXIT_SUCCESS;
}

#undef MC_2SAP_CREATOR_HASH_SIZE
