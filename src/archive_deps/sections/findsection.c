void findsection(double Beta,double LmultR,double kappa)
//this function find the section which has the greatest expected number of occurrences.
{
	printf("Finding section with the highest expected number of occurrences.\n");
	int i;
	double densityi, maxdensity;
	int maxsection;
	double Enum_Enum;
	maxdensity=0;
	maxsection=1;

	for(i=1; i<=max_sections; i++){
		Enum_Enum = (L_Evector[0][i]) * (R_Evector[0][i]) / LmultR;
		densityi = Enum_Enum/ Beta;
		if(densityi>maxdensity){
			maxdensity=densityi;
			maxsection=i;
		}
	}
	printf("Section %d has the highest expected number of occurrences (%f) = (%e).\n", maxsection, maxdensity, maxdensity);
	printf("without dividing by beta, (prob of this section=%f=%e)\n", maxdensity*Beta, maxdensity*Beta);

	printsection(maxsection);


}
