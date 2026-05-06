#ifndef MC_2SAP_COMMON_H
#define MC_2SAP_COMMON_H

void mc_2sap_set_system_params(int is_hamiltonian);
void mc_2sap_ignore_system_result(const char *command);
int ***mc_2sap_alloc_int3_table(unsigned long entries, int components, const unsigned int *widths, const char *label);
int ***mc_2sap_alloc_int3_fixed(unsigned long entries, int components, unsigned int width, const char *label);

#endif
