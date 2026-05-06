#ifndef MC_2SAP_COMMON_H
#define MC_2SAP_COMMON_H

#include <stddef.h>
#include <stdio.h>

typedef struct Mc2SapPairMap {
	unsigned long int *first;
	unsigned long int *second;
	unsigned long int *value;
	size_t capacity;
} Mc2SapPairMap;

typedef struct Mc2SapModeSpec {
	int hamiltonian;
	const char *sample_outdir;
	const char *sample_prefix;
	const char *creator_outdir;
	const char *creator_prefix;
	const char *creator_object_label;
} Mc2SapModeSpec;

typedef struct Mc2SapCreatorConfig {
	unsigned long int max_keynum;
	unsigned long int *num_left_endhinges;
	unsigned long int *num_right_endhinges;
	unsigned long int *num_outsections;
	unsigned long int **t_outsection;
	int totalspan;
	int L;
	int M;
	int mode;
	unsigned int seednum;
	int runnum;
	int maxpolys;
	unsigned long int limit;
	int force;
	const char *outdir;
	const char *prefix;
	const char *object_label;
	void *user;
	void (*reset_built_walks)(void *user);
	void (*load_left_endhinge)(void *user, unsigned long int section, int nth_endhinge);
	void (*add_transition)(void *user, unsigned long int section, int nth_tspan);
	void (*add_right_endhinge)(void *user, unsigned long int section, int nth_endhinge);
	char *(*pair_key)(void *user);
	void (*print_pair)(void *user, FILE *out);
} Mc2SapCreatorConfig;

typedef struct Mc2SapSampleWriterConfig {
	FILE **fp;
	char *filename;
	size_t filename_size;
	unsigned long int *filetotal;
	unsigned long int *filenum;
	const char *outdir;
	const char *prefix;
	const int *poly1_start;
	const int *poly1_directions;
	const int *poly2_start;
	const int *poly2_directions;
	int maxpolys;
	int mode;
	int L;
	int M;
	int totalspan;
	unsigned int seednum;
	int runnum;
	double dom_evalue;
} Mc2SapSampleWriterConfig;

void mc_2sap_set_system_params(int is_hamiltonian);
const Mc2SapModeSpec *mc_2sap_mode_spec(int is_hamiltonian);
void mc_2sap_ignore_system_result(const char *command);
int ***mc_2sap_alloc_int3_table(unsigned long entries, int components, const unsigned int *widths, const char *label);
int ***mc_2sap_alloc_int3_fixed(unsigned long entries, int components, unsigned int width, const char *label);
void mc_2sap_init_int_rows(int **rows, size_t row_count, size_t width, const char *label);
int **mc_2sap_alloc_int_rows(size_t row_count, size_t width, const char *label);
int mc_2sap_reverse_direction(int direction, const char *caller_name);
void mc_2sap_open_sample_file(const Mc2SapSampleWriterConfig *config);
void mc_2sap_write_sample_pair(const Mc2SapSampleWriterConfig *config);
char *mc_2sap_walk_string(const int *start, const int *directions);
char *mc_2sap_unordered_pair_key(
	const int *poly1_start,
	const int *poly1_directions,
	const int *poly2_start,
	const int *poly2_directions);
void mc_2sap_print_walk(FILE *out, const int *start, const int *directions);
void mc_2sap_print_unordered_pair(
	FILE *out,
	const int *poly1_start,
	const int *poly1_directions,
	const int *poly2_start,
	const int *poly2_directions);
double mc_2sap_max_eval(
	unsigned long int max_keynum,
	unsigned long int max_tspans,
	unsigned long int *num_outsections,
	unsigned long int **tspans_outsection,
	unsigned long int **tspans_edges,
	unsigned long int **tspans_nrr,
	double *left[2],
	double *right[2],
	double force,
	int L,
	int M,
	int hamiltonian,
	double fugacity);
void mc_2sap_pair_map_init(Mc2SapPairMap *map, unsigned long int max_entries);
void mc_2sap_pair_map_free(Mc2SapPairMap *map);
void mc_2sap_pair_map_put(Mc2SapPairMap *map, unsigned long int first, unsigned long int second, unsigned long int value);
unsigned long int mc_2sap_pair_map_get(const Mc2SapPairMap *map, unsigned long int first, unsigned long int second);
void mc_2sap_pair_map_build(Mc2SapPairMap *map, unsigned long int max_entries, unsigned long int pairs[][2]);
int mc_2sap_run_creator_all(const Mc2SapCreatorConfig *config);

#endif
