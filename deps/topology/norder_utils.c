/* norder_utils.c - Robust State Mapping */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int M, L;
extern int vM, vL;
extern unsigned int ordertemplate[2][17][17];
extern unsigned int ordertemplate2[2][17][17];
unsigned int reordertemplate[2][17][17];
unsigned int reordertemplate2[2][17][17];

typedef struct StateNode {
    unsigned int signature[256];
    unsigned long int id;
    struct StateNode *next;
} StateNode;

static StateNode *state_hash[65536];
static unsigned long int state_counter = 0;

unsigned long int get_state_id(unsigned int *sig, int len) {
    unsigned int hash = 0;
    for (int i = 0; i < len; i++) hash = (hash * 31) + sig[i];
    hash %= 65536;

    StateNode *node = state_hash[hash];
    while (node) {
        if (memcmp(node->signature, sig, len * sizeof(unsigned int)) == 0) return node->id;
        node = node->next;
    }

    StateNode *new_node = (StateNode *)malloc(sizeof(StateNode));
    memcpy(new_node->signature, sig, len * sizeof(unsigned int));
    new_node->id = ++state_counter;
    new_node->next = state_hash[hash];
    state_hash[hash] = new_node;
    return new_node->id;
}

void fillreordertemplate(int ledges, int redges) {
    int i, j, ii, jj, s;
    int firstentry, secondentry;
    int connectingedge;

    for (i = 0; i <= M; i++) {
        for (j = 0; j <= L; j++) {
            reordertemplate[0][i][j] = 0;
            reordertemplate[1][i][j] = 0;
        }
    }

    for (s = 0; s <= 1; s++) {
        firstentry = 1;
        int max_e = (s == 0) ? ledges : redges;
        for (i = 0; i <= M; i++) {
            for (j = 0; j <= L; j++) {
                if (ordertemplate[s][i][j] % 2 != 0 && reordertemplate[s][i][j] == 0) {
                    reordertemplate[s][i][j] = firstentry;
                    if (ordertemplate[s][i][j] == 1) connectingedge = max_e;
                    else connectingedge = ordertemplate[s][i][j] - 1;

                    if (firstentry == 1) secondentry = max_e;
                    else secondentry = firstentry + 1;

                    for (ii = 0; ii <= M; ii++) {
                        for (jj = 0; jj <= L; jj++) {
                            if (ordertemplate[s][ii][jj] == connectingedge) {
                                reordertemplate[s][ii][jj] = secondentry;
                                if (firstentry == 1) firstentry++;
                                else firstentry += 2;
                                goto next_edge;
                            }
                        }
                    }
                } else if (ordertemplate[s][i][j] > 0 && reordertemplate[s][i][j] == 0) {
                    reordertemplate[s][i][j] = firstentry;
                    if ((s == 0 && ordertemplate[s][i][j] == ledges) || (s == 1 && ordertemplate[s][i][j] == redges))
                        connectingedge = 1;
                    else connectingedge = ordertemplate[s][i][j] + 1;

                    if (firstentry == 1) secondentry = max_e;
                    else secondentry = firstentry + 1;

                    for (ii = 0; ii <= M; ii++) {
                        for (jj = 0; jj <= L; jj++) {
                            if (ordertemplate[s][ii][jj] == connectingedge) {
                                reordertemplate[s][ii][jj] = secondentry;
                                if (firstentry == 1) firstentry++;
                                else firstentry += 2;
                                goto next_edge;
                            }
                        }
                    }
                }
                next_edge:;
            }
        }
    }
}

void fillreordertemplate2(int ledges, int redges) {
    int i, j, ii, jj, s;
    int firstentry, secondentry;
    int connectingedge;

    for (i = 0; i <= M; i++) {
        for (j = 0; j <= L; j++) {
            reordertemplate2[0][i][j] = 0;
            reordertemplate2[1][i][j] = 0;
        }
    }

    for (s = 0; s <= 1; s++) {
        firstentry = 1;
        int max_e = (s == 0) ? ledges : redges;
        for (i = 0; i <= M; i++) {
            for (j = 0; j <= L; j++) {
                if (ordertemplate2[s][i][j] % 2 != 0 && reordertemplate2[s][i][j] == 0) {
                    reordertemplate2[s][i][j] = firstentry;
                    if (ordertemplate2[s][i][j] == 1) connectingedge = max_e;
                    else connectingedge = ordertemplate2[s][i][j] - 1;

                    if (firstentry == 1) secondentry = max_e;
                    else secondentry = firstentry + 1;

                    for (ii = 0; ii <= M; ii++) {
                        for (jj = 0; jj <= L; jj++) {
                            if (ordertemplate2[s][ii][jj] == connectingedge) {
                                reordertemplate2[s][ii][jj] = secondentry;
                                if (firstentry == 1) firstentry++;
                                else firstentry += 2;
                                goto next_edge2;
                            }
                        }
                    }
                } else if (ordertemplate2[s][i][j] > 0 && reordertemplate2[s][i][j] == 0) {
                    reordertemplate2[s][i][j] = firstentry;
                    if ((s == 0 && ordertemplate2[s][i][j] == ledges) || (s == 1 && ordertemplate2[s][i][j] == redges))
                        connectingedge = 1;
                    else connectingedge = ordertemplate2[s][i][j] + 1;

                    if (firstentry == 1) secondentry = max_e;
                    else secondentry = firstentry + 1;

                    for (ii = 0; ii <= M; ii++) {
                        for (jj = 0; jj <= L; jj++) {
                            if (ordertemplate2[s][ii][jj] == connectingedge) {
                                reordertemplate2[s][ii][jj] = secondentry;
                                if (firstentry == 1) firstentry++;
                                else firstentry += 2;
                                goto next_edge2;
                            }
                        }
                    }
                }
                next_edge2:;
            }
        }
    }
}

unsigned long int num_section_norder(int side) {
    unsigned int sig[256] = {0};
    int idx = 0;
    for (int i = 0; i <= M; i++) for (int j = 0; j <= L; j++) sig[idx++] = reordertemplate[side][i][j];
    return get_state_id(sig, idx);
}

unsigned long int num_section_norder2(int side) {
    unsigned int sig[256] = {0};
    int idx = 0;
    for (int i = 0; i <= M; i++) for (int j = 0; j <= L; j++) sig[idx++] = reordertemplate2[side][i][j];
    return get_state_id(sig, idx);
}
