extern unsigned int reordertemplate[2][17][17];

unsigned short int LFlag_endhinge(int (*EndOrdNum)[3]) {
    if ((*EndOrdNum)[0] - 1 < vM * vL) {
        return 1;
    } else {
        unsigned short int i, j, b;

        for (i = 0; i <= M; i++) {
            for (j = 0; j < L; j++) {
                int diff = reordertemplate[0][i][j] - (b = reordertemplate[0][i][j + 1]);
                if (diff == 1) {
                    if (b != 0 && (b % 2 == 0)) return 1;
                } else if (diff == -1) {
                    if (b != 1 && (b % 2 == 1)) return 1;
                } else if (diff == (vM * vL - 1) || diff == (1 - vM * vL)) {
                    return 1;
                }
            }
        }

        for (i = 0; i < M; i++) {
            for (j = 0; j <= L; j++) {
                int diff = reordertemplate[0][i][j] - (b = reordertemplate[0][i + 1][j]);
                if (diff == 1) {
                    if (b != 0 && (b % 2 == 0)) return 1;
                } else if (diff == -1) {
                    if (b != 1 && (b % 2 == 1)) return 1;
                } else if (diff == (vM * vL - 1) || diff == (1 - vM * vL)) {
                    return 1;
                }
            }
        }
    }
    return 0;
}
