#ifndef MC_2SAP_COMMON_H
#define MC_2SAP_COMMON_H

#include <stddef.h>

typedef struct Mc2SapPairMap {
	unsigned long int *first;
	unsigned long int *second;
	unsigned long int *value;
	size_t capacity;
} Mc2SapPairMap;

void mc_2sap_set_system_params(int is_hamiltonian);
void mc_2sap_ignore_system_result(const char *command);
int ***mc_2sap_alloc_int3_table(unsigned long entries, int components, const unsigned int *widths, const char *label);
int ***mc_2sap_alloc_int3_fixed(unsigned long entries, int components, unsigned int width, const char *label);
void mc_2sap_pair_map_init(Mc2SapPairMap *map, unsigned long int max_entries);
void mc_2sap_pair_map_free(Mc2SapPairMap *map);
void mc_2sap_pair_map_put(Mc2SapPairMap *map, unsigned long int first, unsigned long int second, unsigned long int value);
unsigned long int mc_2sap_pair_map_get(const Mc2SapPairMap *map, unsigned long int first, unsigned long int second);
void mc_2sap_pair_map_build(Mc2SapPairMap *map, unsigned long int max_entries, unsigned long int pairs[][2]);

#endif
