unsigned short int
LFlag(int (*pointordNum)[])
{
    int *pON = (int*)pointordNum;
	if (pON[0] - 1 < vM*vL) {
		return (1);
	} else {
		unsigned short int i, j, b;
		for (i = 0; i <= lat_M; i++) {
			for (j = 0; j < lat_L; j++) {
				int diff = (ordertemplate[0][i][j]) - (b=ordertemplate[0][i][j + 1]);
				if (diff == 1) { if (b!= 0 && (b & 1) == 0) return 1; }
				else if (diff == -1) { if (b!=1 && (b & 1) == 1) return 1; }
				else if (diff == vM*vL-1 || diff == 1-vM*vL) return 1;
			}
		}
		for (i = 0; i < lat_M; i++) {
			for (j = 0; j <= lat_L; j++) {
				int diff = (ordertemplate[0][i][j]) - (b=ordertemplate[0][i+1][j]);
				if (diff == 1) { if (b!= 0 && (b & 1) == 0) return 1; }
				else if (diff == -1) { if (b!=1 && (b & 1) == 1) return 1; }
				else if (diff == vM*vL-1 || diff == 1-vM*vL) return 1;
			}
		}
	}
	return (0);
}
