#include <stdio.h>
#include <stdlib.h>

#include "mc_2sap_common.h"
#include "mc_runtime.h"

extern int L;
extern int M;
extern int max_sections;
extern int max_tspans;
extern unsigned long int max_keynum;
extern double dom_evalue;

void mc_2sap_set_system_params(int is_hamiltonian)
{
	if (is_hamiltonian) {
		if (M == 1 && L == 2) {
			max_sections = 73;
			max_keynum = 152;
			max_tspans = 490;
			dom_evalue = 5.534148126030995;
		} else if (M == 1 && L == 3) {
			max_sections = 742;
			max_keynum = 4048;
			max_tspans = 37454;
			dom_evalue = 24.378235811209002;
		} else if (M == 1 && L == 4) {
			max_sections = 9309;
			max_keynum = 108386;
			max_tspans = 2598620;
			dom_evalue = 97.848097677835298;
		} else if (M == 2 && L == 2) {
			max_sections = 2619;
			max_keynum = 21540;
			max_tspans = 495792;
			dom_evalue = 53.686533141903084;
		} else {
			fprintf(stderr, "Unsupported Hamiltonian 2SAP L and M values (%d, %d).\n", L, M);
			exit(EXIT_FAILURE);
		}
	} else {
		if (M == 1 && L == 2) {
			max_sections = 73;
			max_keynum = 152;
			max_tspans = 1048;
			dom_evalue = 9.455960990693537;
		} else if (M == 1 && L == 3) {
			max_sections = 742;
			max_keynum = 4048;
			max_tspans = 108410;
			dom_evalue = 65.012508345749453;
		} else if (M == 1 && L == 4) {
			max_sections = 9309;
			max_keynum = 108386;
			max_tspans = 9419930;
			dom_evalue = 380.824254675661734;
		} else if (M == 2 && L == 2) {
			max_sections = 2619;
			max_keynum = 21546;
			max_tspans = 1126578;
			dom_evalue = 180.511702878389116;
		} else {
			fprintf(stderr, "Unsupported 2SAP L and M values (%d, %d).\n", L, M);
			exit(EXIT_FAILURE);
		}
	}
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
