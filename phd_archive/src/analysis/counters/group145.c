/*
group 22: 25, 35, 36, 39, 45, 1705, 1716, 1726, 1730, 1742, 2049, 2055, 2070, 2084, 2094, 3476, 3479, 3482, 3488, 3511, 5887, 5910, 5916, 5919, 5922, 7322, 7333, 7347, 7363, 7368, 7674, 7687, 7692, 7701, 7712, 9353, 9359, 9362, 9363, 9373, 

group 145: 206, 219, 222, 589, 596, 605, 827, 833, 1266, 1283, 3366, 3373, 3376, 3789, 3798, 3801, 4215, 4225, 4600, 4617, 5161, 5178, 5530, 5541, 5962, 5965, 5969, 6401, 6402, 6409, 8472, 8490, 8945, 8951, 9172, 9179, 9186, 9541, 9544, 9552, 

group 330: 514, 524, 529, 887, 891, 2002, 2016, 2018, 2238, 2244, 2248, 2523, 2536, 2993, 3008, 3441, 3450, 3457, 4540, 4559, 5219, 5238, 6316, 6321, 6330, 6746, 6760, 7218, 7230, 7517, 7521, 7522, 7747, 7749, 7758, 8887, 8891, 9244, 9247, 9257,
*/

void group145(double Beta,double LmultR,double kappa) {
	//outputs the number of edges in each 2-span of group 145

	int i, j, k;

	int tspans[40];
	tspans[0]=206;
	tspans[1]=219;
	tspans[3]=222;
	tspans[4]=589;
	tspans[5]=596;
	tspans[6]=605;
	tspans[7]=827;
	tspans[8]=833;
	tspans[9]=1266;
	tspans[10]=1283;
	tspans[11]=3366;
	tspans[12]=3373;
	tspans[13]=3789;
	tspans[14]=3798;
	tspans[15]=3801;
	tspans[16]=4215;
	tspans[17]=4225;
	tspans[18]=4600;
	tspans[19]=4617;
	tspans[20]=5161;
	tspans[21]=5178;
	tspans[22]=5530;
	tspans[23]=5541;
	tspans[24]=5962;
	tspans[25]=5965;
	tspans[26]=5969;
	tspans[27]=6401;
	tspans[28]=6402;
	tspans[29]=6409;
	tspans[30]=8472;
	tspans[31]=8490;
	tspans[32]=8945;
	tspans[33]=8951;
	tspans[34]=9172;
	tspans[35]=9179;
	tspans[36]=9186;
	tspans[37]=9541;
	tspans[38]=9544;
	tspans[39]=9552;

	printf("group 145:\n");
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










