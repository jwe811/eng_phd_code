#include <stdio.h>
#include <stdlib.h>
#include <string.h>

main()
{
        FILE *g;
        char filename[16];
        int  i=1,j=1,k,l,m;

		for(j=1;j<=400;j++)

		{
                strcpy(filename,"");
                sprintf(filename,"script_%d.pbs",j);
                g=fopen(filename,"w");

		fprintf(g, "#PBS -S /bin/bash\n");
		fprintf(g, "#PBS -N file_%d\n",i);
		fprintf(g, "#PBS -l nodes=1\n");
		fprintf(g, "#PBS -l walltime=0:05:00,mem=3000mb\n");
		fprintf(g, "#PBS -m ea\n");
		fprintf(g, "cd $PBS_O_WORKDIR\n");
	
                //fprintf(g, "#!/bin/bash\n");
			//printf("on: %d\n", j);
                        fprintf(g, "./t_run CreatePrimePolys2M1span8num%d.txt . \n", j);

                	fclose(g);
		}


}

