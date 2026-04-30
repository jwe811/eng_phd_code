/*
group 22: 25, 35, 36, 39, 45, 1705, 1716, 1726, 1730, 1742, 2049, 2055, 2070, 2084, 2094, 3476, 3479, 3482, 3488, 3511, 5887, 5910, 5916, 5919, 5922, 7322, 7333, 7347, 7363, 7368, 7674, 7687, 7692, 7701, 7712, 9353, 9359, 9362, 9363, 9373, 

group 145: 206, 219, 222, 589, 596, 605, 827, 833, 1266, 1283, 3366, 3373, 3376, 3789, 3798, 3801, 4215, 4225, 4600, 4617, 5161, 5178, 5530, 5541, 5962, 5965, 5969, 6401, 6402, 6409, 8472, 8490, 8945, 8951, 9172, 9179, 9186, 9541, 9544, 9552, 

group 330: 514, 524, 529, 887, 891, 2002, 2016, 2018, 2238, 2244, 2248, 2523, 2536, 2993, 3008, 3441, 3450, 3457, 4540, 4559, 5219, 5238, 6316, 6321, 6330, 6746, 6760, 7218, 7230, 7517, 7521, 7522, 7747, 7749, 7758, 8887, 8891, 9244, 9247, 9257,
*/

void group330(double Beta,double LmultR,double kappa) {
	//outputs the number of edges in each 2-span of group 330

	int i, j, k;

	int tspans[40];
	tspans[0]=514;
	tspans[1]=524;
	tspans[3]=529;
	tspans[4]=887;
	tspans[5]=891;
	tspans[6]=2002;
	tspans[7]=2016;
	tspans[8]=2018;
	tspans[9]=2238;
	tspans[10]=2244;
	tspans[11]=2248;
	tspans[12]=2523;
	tspans[13]=2536;
	tspans[14]=2993;
	tspans[15]=3008;
	tspans[16]=3441;
	tspans[17]=3450;
	tspans[18]=3457;
	tspans[19]=4540;
	tspans[20]=4559;
	tspans[21]=5219;
	tspans[22]=5238;
	tspans[23]=6316;
	tspans[24]=6321;
	tspans[25]=6330;
	tspans[26]=6746;
	tspans[27]=6760;
	tspans[28]=7218;
	tspans[29]=7230;
	tspans[30]=7517;
	tspans[31]=7521;
	tspans[32]=7522;
	tspans[33]=7747;
	tspans[34]=7749;
	tspans[35]=7758;
	tspans[36]=8887;
	tspans[37]=8891;
	tspans[38]=9244;
	tspans[39]=9257;

	printf("group 330:\n");
	for(i=0; i<=39; i++){
		//find number of edges in tspans[i]
		//find the appropriate 2-span section and outsection number
		for(j=1; j<=max_sections; j++){
			if(num_outsections[j]>0){
				for(k=1; k<=num_outsections[j]; k++){
					if(tspans_nrr[j][k]==tspans[i]){	//found the 2-span
						printf("Number of edges in tspan %d=%d is %d\n", tspans_nrr[j][k], tspans[i], tspans_edges[j][k]);
					}
				}
			}
		}
	}
}










