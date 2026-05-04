#include "mc_globals.h"

void check_reachability() {
    if (!ham_check) return;
    
    printf("\n--- HAMILTONIAN REACHABILITY CHECK ---\n");
    unsigned long int testsection = 1; 
    unsigned char *testvector = (unsigned char*)calloc(max_sections + 1, sizeof(unsigned char));
    unsigned char *temptestvector = (unsigned char*)calloc(max_sections + 1, sizeof(unsigned char));
    unsigned char *testvector2 = NULL;
    unsigned char *temptestvector2 = NULL;
    
    if (M == 2 && L == 2) {
        testvector2 = (unsigned char*)calloc(max_sections + 1, sizeof(unsigned char));
        temptestvector2 = (unsigned char*)calloc(max_sections + 1, sizeof(unsigned char));
    }

    testvector[testsection] = 1;
    unsigned long int testpower = 200; // Sufficient for most reachable spaces
    
    printf("Performing reachability analysis (power=%lu)...\n", testpower);
    for (unsigned long int k = 1; k <= testpower; k++) {
        for (int i = 1; i <= max_sections; i++) temptestvector[i] = 0;
        for (int i = 1; i <= max_sections; i++) {
            if (testvector[i]) {
                for (int j = 1; j <= num_outsections[i]; j++) {
                    unsigned long int out = t_outsection[i][j];
                    if (out > 0) temptestvector[out] = 1;
                }
            }
        }
        int changed = 0;
        for (int i = 1; i <= max_sections; i++) {
            if (testvector[i] != temptestvector[i]) changed = 1;
            testvector[i] = temptestvector[i];
        }
        if (!changed && k > 10) break; // Converged
    }
    
    if (M == 2 && L == 2) {
        testvector2[testsection] = 1;
        for (unsigned long int k = 1; k <= testpower + 1; k++) {
            for (int i = 1; i <= max_sections; i++) temptestvector2[i] = 0;
            for (int i = 1; i <= max_sections; i++) {
                if (testvector2[i]) {
                    for (int j = 1; j <= num_outsections[i]; j++) {
                        unsigned long int out = t_outsection[i][j];
                        if (out > 0) temptestvector2[out] = 1;
                    }
                }
            }
            int changed = 0;
            for (int i = 1; i <= max_sections; i++) {
                if (testvector2[i] != temptestvector2[i]) changed = 1;
                testvector2[i] = temptestvector2[i];
            }
            if (!changed && k > 10) break;
        }
    }

    unsigned long int numconnect = 0;
    unsigned long int numnoconnect = 0;
    for (int i = 1; i <= max_sections; i++) {
        if (num_outsections[i] > 0) {
            int reachable = 0;
            if (M == 2 && L == 2) reachable = (testvector[i] || testvector2[i]);
            else reachable = testvector[i];
            
            if (!reachable) {
                numnoconnect++;
                tot_left_endhinges -= num_left_endhinges[i];
                tot_right_endhinges -= num_right_endhinges[i];
                // Invalidate transitions
                for (int j = 1; j <= num_outsections[i]; j++) {
                    t_outsection[i][j] = 0;
                    t_num_walks[i][j] = 0;
                }
                num_outsections[i] = 0;
                num_left_endhinges[i] = 0;
                num_right_endhinges[i] = 0;
            } else {
                numconnect++;
            }
        }
    }
    
    // Re-index transitions (t_nrr)
    unsigned long int new_tspan_num = 0;
    for (int i = 1; i <= max_sections; i++) {
        for (int j = 1; j <= num_outsections[i]; j++) {
            t_nrr[i][j] = ++new_tspan_num;
        }
    }
    
    printf("Reachable Hamiltonian sections: %lu (Removed %lu isolated states)\n", numconnect, numnoconnect);
    printf("Total valid Hamiltonian 2-spans: %lu\n", new_tspan_num);
    num_tspans = new_tspan_num;

    free(testvector); free(temptestvector);
    if (testvector2) { free(testvector2); free(temptestvector2); }
}
