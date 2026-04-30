void fillreordertemplate(int s, int ledges, int redges) {
    int i, j, ii, jj, firstentry, secondentry, connectingedge;
    int cur_j = (s == 0) ? 0 : lat_L;
    for(i=0; i<=lat_M; i++) for(j=0; j<=lat_L; j++) reordertemplate[s][i][j]=0;
    int m = (s == 0) ? ledges : redges; if (m <= 0) return;
    firstentry=1;
    for(i=0; i<=lat_M; i++) {
        if(ordertemplate[s][i][cur_j] > 0 && reordertemplate[s][i][cur_j]==0) {
            reordertemplate[s][i][cur_j] = firstentry;
            connectingedge = (ordertemplate[s][i][cur_j] % 2 != 0) ? 
                             ((ordertemplate[s][i][cur_j] == 1) ? m : ordertemplate[s][i][cur_j] - 1) :
                             ((ordertemplate[s][i][cur_j] == m) ? 1 : ordertemplate[s][i][cur_j] + 1);
            secondentry = firstentry + 1;
            for(ii=0; ii<=lat_M; ii++) {
                if(ordertemplate[s][ii][cur_j] == connectingedge) {
                    reordertemplate[s][ii][cur_j] = secondentry;
                    firstentry += 2;
                    break;
                }
            }
        }
    }
}
