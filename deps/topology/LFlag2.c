/* LFlag2.c */
#include <stdio.h>

extern int M, L;
extern int vM, vL;
extern unsigned int ordertemplate[2][17][17];
extern unsigned int ordertemplate2[2][17][17];

unsigned short int LFlag2(int (*pointordNum)[6]) {
    unsigned short int i, j, b;
    if ((*pointordNum)[0] + (*pointordNum)[3] - 2 < vM * vL) {
        return 1;
    } else {
        if (L == 1 && M == 1) {
            if (ordertemplate[0][0][0] != 0 && ordertemplate[0][1][1] == 0) return 1;
            else if (ordertemplate2[0][0][0] != 0 && ordertemplate2[0][1][1] == 0) return 1;
        } else {
            // Check 1st polygon
            if ((*pointordNum)[0] - 1 > 2) {
                for (i = 0; i <= M; i++) {
                    for (j = 0; j < L; j++) {
                        int diff = ordertemplate[0][i][j] - (b = ordertemplate[0][i][j + 1]);
                        if (diff == 1) { if (b != 0 && (b % 2 == 0)) return 1; }
                        else if (diff == -1) { if (b != 1 && (b % 2 == 1)) return 1; }
                        
                        if (ordertemplate[0][i][j] != 0 && ordertemplate[0][i][j + 1] != 0) {
                            if (ordertemplate[0][i][j] - ordertemplate[0][i][j + 1] == (*pointordNum)[0] - 2) return 1;
                            if (ordertemplate[0][i][j] - ordertemplate[0][i][j + 1] == 2 - (*pointordNum)[0]) return 1;
                        }
                    }
                }
                for (i = 0; i < M; i++) {
                    for (j = 0; j <= L; j++) {
                        int diff = ordertemplate[0][i][j] - (b = ordertemplate[0][i + 1][j]);
                        if (diff == 1) { if (b != 0 && (b % 2 == 0)) return 1; }
                        else if (diff == -1) { if (b != 1 && (b % 2 == 1)) return 1; }

                        if (ordertemplate[0][i][j] != 0 && ordertemplate[0][i + 1][j] != 0) {
                            if (ordertemplate[0][i][j] - ordertemplate[0][i + 1][j] == (*pointordNum)[0] - 2) return 1;
                            if (ordertemplate[0][i][j] - ordertemplate[0][i + 1][j] == 2 - (*pointordNum)[0]) return 1;
                        }
                    }
                }
            }
            // Check 2nd polygon
            if ((*pointordNum)[3] - 1 > 2) {
                for (i = 0; i <= M; i++) {
                    for (j = 0; j < L; j++) {
                        int diff = ordertemplate2[0][i][j] - (b = ordertemplate2[0][i][j + 1]);
                        if (diff == 1) { if (b != 0 && (b % 2 == 0)) return 1; }
                        else if (diff == -1) { if (b != 1 && (b % 2 == 1)) return 1; }

                        if (ordertemplate2[0][i][j] != 0 && ordertemplate2[0][i][j + 1] != 0) {
                            if (ordertemplate2[0][i][j] - ordertemplate2[0][i][j + 1] == (*pointordNum)[3] - 2) return 1;
                            if (ordertemplate2[0][i][j] - ordertemplate2[0][i][j + 1] == 2 - (*pointordNum)[3]) return 1;
                        }
                    }
                }
                for (i = 0; i < M; i++) {
                    for (j = 0; j <= L; j++) {
                        int diff = ordertemplate2[0][i][j] - (b = ordertemplate2[0][i + 1][j]);
                        if (diff == 1) { if (b != 0 && (b % 2 == 0)) return 1; }
                        else if (diff == -1) { if (b != 1 && (b % 2 == 1)) return 1; }

                        if (ordertemplate2[0][i][j] != 0 && ordertemplate2[0][i + 1][j] != 0) {
                            if (ordertemplate2[0][i][j] - ordertemplate2[0][i + 1][j] == (*pointordNum)[3] - 2) return 1;
                            if (ordertemplate2[0][i][j] - ordertemplate2[0][i + 1][j] == 2 - (*pointordNum)[3]) return 1;
                        }
                    }
                }
            }
        }
    }
    return 0;
}
