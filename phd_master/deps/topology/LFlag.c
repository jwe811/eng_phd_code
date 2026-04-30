/* LFlag.c */
#include <stdio.h>

extern int M, L;
extern int vM, vL;
extern unsigned int ordertemplate[2][17][17];

unsigned short int LFlag(int (*pointordNum)[6]) {
    if ((*pointordNum)[0] - 1 < vM * vL) {
        return 1;
    } else {
        unsigned short int i, j, b;
        for (i = 0; i <= M; i++) {
            for (j = 0; j < L; j++) {
                int diff = ordertemplate[0][i][j] - (b = ordertemplate[0][i][j + 1]);
                if (diff == 1) { if (b != 0 && (b % 2 == 0)) return 1; }
                else if (diff == -1) { if (b != 1 && (b % 2 == 1)) return 1; }
            }
        }
        for (i = 0; i < M; i++) {
            for (j = 0; j <= L; j++) {
                int diff = ordertemplate[0][i][j] - (b = ordertemplate[0][i + 1][j]);
                if (diff == 1) { if (b != 0 && (b % 2 == 0)) return 1; }
                else if (diff == -1) { if (b != 1 && (b % 2 == 1)) return 1; }
            }
        }
    }
    return 0;
}
