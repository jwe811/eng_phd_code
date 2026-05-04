double expoccur(double Beta,double LmultR,double kappa, int tspannum) {
	//this function will return the expected number of occurrences of tspannum (which is an input)
	//only works in TS mode

	double Enum_Enum = (L_Evector[0][tspannum]) * (R_Evector[0][tspannum]) / LmultR;

	return(Enum_Enum / Beta);

}
