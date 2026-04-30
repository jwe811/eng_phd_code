unsigned short int
LFlag2(int (*pointordNum)[])
{
	unsigned short int i, j, b;
    int *pON = (int*)pointordNum;
	if (pON[0]+pON[3] - 2 < vM*vL) {
		return 1;
	} else {
		if(lat_L==1 && lat_M==1){
			if(ordertemplate[0][0][0] != 0 && ordertemplate[0][1][1] == 0) return 1;
			else if(ordertemplate2[0][0][0] != 0 && ordertemplate2[0][1][1] == 0) return 1;
		} else {
			if(pON[0]-1 > 2){
				for (i = 0; i <= lat_M; i++) {
					for (j = 0; j < lat_L; j++) {
						int diff = (ordertemplate[0][i][j]) - (b=ordertemplate[0][i][j + 1]);
						if (diff == 1) { if (b!= 0 && (b & 1) == 0) return 1; }
						else if (diff == -1) { if (b!=1 && (b & 1) == 1) return 1; }
						if(ordertemplate[0][i][j]!=0 && ordertemplate[0][i][j + 1]!=0){
							if(ordertemplate[0][i][j] - ordertemplate[0][i][j + 1]==pON[0]-1-1) return 1;
							else if(ordertemplate[0][i][j] - ordertemplate[0][i][j + 1]==1-(pON[0]-1)) return 1;
						}
					}
				}
				for (i = 0; i < lat_M; i++) {
					for (j = 0; j <= lat_L; j++) {
						int diff = (ordertemplate[0][i][j]) - (b=ordertemplate[0][i+1][j]);
						if (diff == 1) { if (b!= 0 && (b & 1) == 0) return 1; }
						else if (diff == -1) { if (b!=1 && (b & 1) == 1) return 1; }
						if(ordertemplate[0][i][j]!=0 && ordertemplate[0][i+1][j]!=0){
							if(ordertemplate[0][i][j] - ordertemplate[0][i+1][j]==pON[0]-1-1) return 1;
							else if(ordertemplate[0][i][j] - ordertemplate[0][i+1][j]==1-(pON[0]-1)) return 1;
						}
					}
				}
			}
			if(pON[3]-1 > 2){
				for (i = 0; i <= lat_M; i++) {
					for (j = 0; j < lat_L; j++) {
						int diff = (ordertemplate2[0][i][j]) - (b=ordertemplate2[0][i][j + 1]);
						if (diff == 1) { if (b!= 0 && (b & 1) == 0) return 1; }
						else if (diff == -1) { if (b!=1 && (b & 1) == 1) return 1; }
						if(ordertemplate2[0][i][j]!=0 && ordertemplate2[0][i][j + 1]!=0){
							if(ordertemplate2[0][i][j] - ordertemplate2[0][i][j + 1]==pON[3]-1-1) return 1;
							else if(ordertemplate2[0][i][j] - ordertemplate2[0][i][j + 1]==1-(pON[3]-1)) return 1;
						}
					}
				}
				for (i = 0; i < lat_M; i++) {
					for (j = 0; j <= lat_L; j++) {
						int diff = (ordertemplate2[0][i][j]) - (b=ordertemplate2[0][i+1][j]);
						if (diff == 1) { if (b!= 0 && (b & 1) == 0) return 1; }
						else if (diff == -1) { if (b!=1 && (b & 1) == 1) return 1; }
						if(ordertemplate2[0][i][j]!=0 && ordertemplate2[0][i+1][j]!=0){
							if(ordertemplate2[0][i][j] - ordertemplate2[0][i+1][j]==pON[3]-1-1) return 1;
							else if(ordertemplate2[0][i][j] - ordertemplate2[0][i+1][j]==1-(pON[3]-1)) return 1;
						}
					}
				}
			}
		}
	}
	return (0);
}
