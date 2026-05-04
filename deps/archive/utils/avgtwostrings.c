double avgtwostrings(double Beta,double LmultR, double kappa, double fval) {
	//ALGORITHM:
	//loop through all possible sections with two edges
	//find it's section number
	//find that section's exp num occur.

	int i, j;
	int i1, j1, i2, j2, i3, j3, i4, j4, i5, j5, i6, j6;
	int sec_num;
	double total=0.0;
	int num_two_secs=0;
	int num_four_secs=0;
	int num_six_secs=0;

	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			ordertemplate[0][i][j]=0;	//initialize ordertemplate to be empty
		}
	}

	Beta=1;	//=1 if you don't want to divide by beta

	//CASE 1: 2edges in section
	for(i1=0; i1<=M; i1++){
		for(j1=0; j1<=L; j1++){
			ordertemplate[0][i1][j1]=1;	//loop through all possible locations for the first edge
			for(i2=0; i2<=M; i2++){
				for(j2=0; j2<=L; j2++){		//loop through all possible locations for the second edge
					if((i1!=i2) || (j1!=j2)){	//that isn't in the location where the first edge is.
						ordertemplate[0][i2][j2]=2;
						//section generated
						sec_num = num_section(0);
						if(num_outsections[sec_num]>1){
							num_two_secs++;
	//						printf("section %d has 2 edges. Find it's exp. num. occur.\n", sec_num);
							printf("sec_num=%d, total=%f, adding %f\n", sec_num, total, L_Evector[0][sec_num] * R_Evector[0][sec_num] / Beta / LmultR);
							total += L_Evector[0][sec_num] * R_Evector[0][sec_num] / Beta / LmultR;
						}

						ordertemplate[0][i2][j2]=0;	//reset to empty
					}
				}
			}
			ordertemplate[0][i1][j1]=0;
		}
	}









	//CASE 2: 4edges in section
	for(i1=0; i1<=M; i1++){
		for(j1=0; j1<=L; j1++){
			ordertemplate[0][i1][j1]=1;	//loop through all possible locations for the first edge
			for(i2=0; i2<=M; i2++){
				for(j2=0; j2<=L; j2++){		//loop through all possible locations for the second edge
					if((i1!=i2) || (j1!=j2)){	//that isn't in the location where the first edge is.
						ordertemplate[0][i2][j2]=2;
						for(i3=0; i3<=M; i3++){
							for(j3=0; j3<=L; j3++){		//loop through all possible locations for the third edge
								if((i3!=i1 || j3!=j1) && (i3!=i2 || j3!=j2)){	//that isn't the location of the 1st or 2nd edge
									ordertemplate[0][i3][j3]=3;
									for(i4=0; i4<=M; i4++){
										for(j4=0; j4<=L; j4++){
											if((i4!=i1 || j4!=j1) && (i4!=i2 || j4!=j2) && (i4!=i3 || j4!=j3)){	//that isn't the location of the 1st or 2nd edge
												ordertemplate[0][i4][j4]=4;
												//section generated
												sec_num = num_section(0);
												if(num_outsections[sec_num]>1){
													num_four_secs++;
							//						printf("section %d has 4 edges. Find it's exp. num. occur.\n", sec_num);
													printf("sec_num=%d, total=%f, adding %f\n", sec_num, total, L_Evector[0][sec_num] * R_Evector[0][sec_num] / Beta / LmultR);
													total += L_Evector[0][sec_num] * R_Evector[0][sec_num] / Beta / LmultR;
												}

												ordertemplate[0][i4][j4]=0;	//reset to empty
											}
										}
									}
									ordertemplate[0][i3][j3]=0;
								}
							}
						}
						ordertemplate[0][i2][j2]=0;
					}
				}
			}
			ordertemplate[0][i1][j1]=0;
		}
	}

	//CASE 3: 6edges in section
	for(i1=0; i1<=M; i1++){
		for(j1=0; j1<=L; j1++){
			ordertemplate[0][i1][j1]=1;	//loop through all possible locations for the 1st edge
			for(i2=0; i2<=M; i2++){
				for(j2=0; j2<=L; j2++){		//loop through all possible locations for the 2nd edge
					if((i1!=i2) || (j1!=j2)){	//that isn't in the location where the 1st edge is.
						ordertemplate[0][i2][j2]=2;
						for(i3=0; i3<=M; i3++){
							for(j3=0; j3<=L; j3++){		//loop through all possible locations for the 3rd edge
								if((i3!=i1 || j3!=j1) && (i3!=i2 || j3!=j2)){	//that isn't the location of the 1st or 2nd edge
									ordertemplate[0][i3][j3]=3;
									for(i4=0; i4<=M; i4++){
										for(j4=0; j4<=L; j4++){		//loop through all possible locations for the 4th edge
											if((i4!=i1 || j4!=j1) && (i4!=i2 || j4!=j2) && (i4!=i3 || j4!=j3)){	//that isn't the location of the 1st, 2nd, or 3rd edge
												ordertemplate[0][i4][j4]=4;
												for(i5=0; i5<=M; i5++){
													for(j5=0; j5<=L; j5++){	//loop through all possible locations for the 5th edge
														if((i5!=i1 || j5!=j1) && (i5!=i2 || j5!=j2) && (i5!=i3 || j5!=j3) && (i5!=i4 || j5!=j4)){	//that isn't the location of the 1st, 2nd, 3rd, or 4th edge
															ordertemplate[0][i5][j5]=5;
															for(i6=0; i6<=M; i6++){
																for(j6=0; j6<=L; j6++){	//loop through all possible locations for the 6th edge
																	if((i6!=i1 || j6!=j1) && (i6!=i2 || j6!=j2) && (i6!=i3 || j6!=j3) && (i6!=i4 || j6!=j4) && (i6!=i5 || j6!=j5)){	//not in spot of others
												//section generated
												sec_num = num_section(0);
												if(num_outsections[sec_num]>1){
													num_six_secs++;
							//						printf("section %d has 6 edges. Find it's exp. num. occur.\n", sec_num);
													printf("sec_num=%d, total=%f, adding %f\n", sec_num, total, L_Evector[0][sec_num] * R_Evector[0][sec_num] / Beta / LmultR);
													total += L_Evector[0][sec_num] * R_Evector[0][sec_num] / Beta / LmultR;
												}

												ordertemplate[0][i6][j6]=0;	//reset to empty
																	}
																}
															}
															ordertemplate[0][i5][j5]=0;
														}
													}
												}
												ordertemplate[0][i4][j4]=0;
											}
										}
									}
									ordertemplate[0][i3][j3]=0;
								}
							}
						}
						ordertemplate[0][i2][j2]=0;
					}
				}
			}
			ordertemplate[0][i1][j1]=0;
		}
	}


	printf("number of 2-strings: %d\n", num_two_secs);
//	printf("number of 4-strings: %d\n", num_four_secs);
//	printf("number of 6-strings: %d\n", num_six_secs);

	return total;
}
