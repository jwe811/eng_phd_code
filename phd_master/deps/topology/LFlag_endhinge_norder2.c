/* LFlag_endhinge_norder2.c */
/* Ports legacy boundary connectivity logic for 2-SAP sampling */

#include <stdio.h>

extern int M, L;
extern int vM, vL;
extern int EndOrdNum[2], EndOrdNum2[2];
extern int endtemplate[16][16];
extern int endtemplate3[16][16];

unsigned short int LFlag_endhinge2() {
    unsigned short int i, j, b;
    if (EndOrdNum[0] + EndOrdNum2[0] - 2 < vM * vL) {
        return 1;
    } else {
        if (L == 1 && M == 1) {
            if (endtemplate[0][0] != 0 && endtemplate[1][1] == 0) return 1;
            else if (endtemplate3[0][0] != 0 && endtemplate3[1][1] == 0) return 1;
        } else {
            // Check 1st polygon
            if (EndOrdNum[0] - 1 > 2) {
                for (i = 0; i <= M; i++) {
                    for (j = 0; j < L; j++) {
                        int diff = endtemplate[i][j] - (b = endtemplate[i][j + 1]);
                        if (diff == 1) { if (b != 0 && (b % 2 == 0)) return 1; }
                        else if (diff == -1) { if (b != 1 && (b % 2 == 1)) return 1; }
                        
                        if (endtemplate[i][j] != 0 && endtemplate[i][j + 1] != 0) {
                            if (endtemplate[i][j] - endtemplate[i][j + 1] == EndOrdNum[0] - 2) return 1;
                            if (endtemplate[i][j] - endtemplate[i][j + 1] == 2 - EndOrdNum[0]) return 1;
                        }
                    }
                }
                for (i = 0; i < M; i++) {
                    for (j = 0; j <= L; j++) {
                        int diff = endtemplate[i][j] - (b = endtemplate[i + 1][j]);
                        if (diff == 1) { if (b != 0 && (b % 2 == 0)) return 1; }
                        else if (diff == -1) { if (b != 1 && (b % 2 == 1)) return 1; }

                        if (endtemplate[i][j] != 0 && endtemplate[i + 1][j] != 0) {
                            if (endtemplate[i][j] - endtemplate[i + 1][j] == EndOrdNum[0] - 2) return 1;
                            if (endtemplate[i][j] - endtemplate[i + 1][j] == 2 - EndOrdNum[0]) return 1;
                        }
                    }
                }
            }
            // Check 2nd polygon
            if (EndOrdNum2[0] - 1 > 2) {
                for (i = 0; i <= M; i++) {
                    for (j = 0; j < L; j++) {
                        int diff = endtemplate3[i][j] - (b = endtemplate3[i][j + 1]);
                        if (diff == 1) { if (b != 0 && (b % 2 == 0)) return 1; }
                        else if (diff == -1) { if (b != 1 && (b % 2 == 1)) return 1; }

                        if (endtemplate3[i][j] != 0 && endtemplate3[i][j + 1] != 0) {
                            if (endtemplate3[i][j] - endtemplate3[i][j + 1] == EndOrdNum2[0] - 2) return 1;
                            if (endtemplate3[i][j] - endtemplate3[i][j + 1] == 2 - EndOrdNum2[0]) return 1;
                        }
                    }
                }
                for (i = 0; i < M; i++) {
                    for (j = 0; j <= L; j++) {
                        int diff = endtemplate3[i][j] - (b = endtemplate3[i + 1][j]);
                        if (diff == 1) { if (b != 0 && (b % 2 == 0)) return 1; }
                        else if (diff == -1) { if (b != 1 && (b % 2 == 1)) return 1; }

                        if (endtemplate3[i][j] != 0 && endtemplate3[i + 1][j] != 0) {
                            if (endtemplate3[i][j] - endtemplate3[i + 1][j] == EndOrdNum2[0] - 2) return 1;
                            if (endtemplate3[i][j] - endtemplate3[i + 1][j] == 2 - EndOrdNum2[0]) return 1;
                        }
                    }
                }
            }
        }
    }
    return 0;
}
