/*
group 22: 25, 35, 36, 39, 45, 1705, 1716, 1726, 1730, 1742, 2049, 2055, 2070, 2084, 2094, 3476, 3479, 3482, 3488, 3511, 5887, 5910, 5916, 5919, 5922, 7322, 7333, 7347, 7363, 7368, 7674, 7687, 7692, 7701, 7712, 9353, 9359, 9362, 9363, 9373, 

group 145: 206, 219, 222, 589, 596, 605, 827, 833, 1266, 1283, 3366, 3373, 3376, 3789, 3798, 3801, 4215, 4225, 4600, 4617, 5161, 5178, 5530, 5541, 5962, 5965, 5969, 6401, 6402, 6409, 8472, 8490, 8945, 8951, 9172, 9179, 9186, 9541, 9544, 9552, 

group 330: 514, 524, 529, 887, 891, 2002, 2016, 2018, 2238, 2244, 2248, 2523, 2536, 2993, 3008, 3441, 3450, 3457, 4540, 4559, 5219, 5238, 6316, 6321, 6330, 6746, 6760, 7218, 7230, 7517, 7521, 7522, 7747, 7749, 7758, 8887, 8891, 9244, 9247, 9257,
*/

void group22(double Beta,double LmultR,double kappa) {
	//outputs the number of edges in each 2-span of group 22

	int i, j, k;

	int tspans[40];
	tspans[0]=25;
	tspans[1]=35;
	tspans[3]=36;
	tspans[4]=39;
	tspans[5]=45;
	tspans[6]=1705;
	tspans[7]=1716;
	tspans[8]=1726;
	tspans[9]=1730;
	tspans[10]=1742;
	tspans[11]=2049;
	tspans[12]=2055;
	tspans[13]=2070;
	tspans[14]=2084;
	tspans[15]=2094;
	tspans[16]=3476;
	tspans[17]=3479;
	tspans[18]=3482;
	tspans[19]=3488;
	tspans[20]=3511;
	tspans[21]=5887;
	tspans[22]=5916;
	tspans[23]=5919;
	tspans[24]=5922;
	tspans[25]=7322;
	tspans[26]=7333;
	tspans[27]=7347;
	tspans[28]=7363;
	tspans[29]=7368;
	tspans[30]=7674;
	tspans[31]=7687;
	tspans[32]=7692;
	tspans[33]=7701;
	tspans[34]=7712;
	tspans[35]=9353;
	tspans[36]=9359;
	tspans[37]=9362;
	tspans[38]=9363;
	tspans[39]=9373;

	printf("group 22:\n");
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










