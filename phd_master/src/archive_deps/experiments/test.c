#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>

main(void){
	int maxpolysize=10000;
	int x,y,z, xx, yy, zz;
	x=1;
	y=1;
	z=0;
	xx=0;
	yy=0;
	zz=0;

	int i, minlength;
	int polygon[maxpolysize];
	int polygon2[maxpolysize];
	for(i=0; i<=maxpolysize-1; i++){
		polygon[i]=0;
		polygon2[i]=0;
	}

	//read in poly
	char filename[100];
	FILE* fp;
	sprintf(filename, "PolyToTest_UofS.txt");
	fp = fopen(filename, "r");
	if(fp==NULL){
	  printf("Error opening '%s'. Exitting\n", filename);
	  exit(1);
	}

	char line[20];
	fgets(line, 20, fp);//get UofS out of the way
	fgets(line, 20, fp);//x,y,z starting coordinates
	sscanf(line, "%d %d %d", &x, &y, &z);

	i=0;
	while(fgets(line, 20, fp)!=NULL){
	  sscanf(line, "%d", &polygon[i]);
	  //	  printf("%d\n", polygon[i]);
	  if(polygon[i]==-111){
	    break;
	  }
	  i++;
	  if(i>maxpolysize){
	    printf("Error. Need to increase maxpolysize. Exitting\n");
	    exit(1);
	  }
	}
	printf("read in %d edges\n", i);


	/*
	polygon[0]=1;
	polygon[1]=1;
	polygon[2]=1;
	polygon[3]=3;
	polygon[4]=3;
	polygon[5]=2;
	polygon[6]=2;
	polygon[7]=2;
	polygon[8]=4;
	polygon[9]=4;

	polygon2[0]=5;
	polygon2[0]=1;
	polygon2[0]=1;
	polygon2[0]=1;
	polygon2[0]=6;
	polygon2[0]=6;
	polygon2[0]=2;
	polygon2[0]=2;
	polygon2[0]=2;
	polygon2[0]=5;
	*/


	/*
	polygon[0]=5;
	polygon[1]=5;
	polygon[2]=5;
	polygon[3]=1;
	polygon[4]=1;
	polygon[5]=6;
	polygon[6]=6;
	polygon[7]=3;
	polygon[8]=2;
	polygon[9]=2;
	polygon[10]=5;
	polygon[11]=5;
	polygon[12]=1;
	polygon[13]=1;
	polygon[14]=1;
	polygon[15]=4;
	polygon[16]=6;
	polygon[17]=6;
	polygon[18]=6;
	polygon[19]=2;
	polygon[20]=2;
	polygon[21]=5;
	polygon[22]=5;
	polygon[23]=3;
	polygon[24]=1;
	polygon[25]=1;
	polygon[26]=6;
	polygon[27]=6;
	polygon[28]=2;
	polygon[29]=2;
	polygon[30]=2;
	polygon[31]=4;
	*/
	i=0;
	while(polygon[i]>0){
		printf("polygon[%d]=%d\n", i, polygon[i]);
		i++;
	}
	printf("polygon length is %d\n", i);
	/*	j=0;
	while(polygon2[j]>0){
		printf("polygon2[%d]=%d\n", j, polygon[j]);
		j++;
	}
	printf("polygon2 length is %d\n", j);

	printf("Feeding polygons into bfacf_shrink_links()\n");
	*/

/////////////////////////////////////
//	minlength=bfacf_shrink_links(x,y,z,polygon, i);
/////////////////////////////////////

	printf("Feeding polygon into bfacf_shrink()\n");

	minlength = bfacf_shrink(x,y,z,polygon, i);

	printf("bfacf_shrink() is finished\n");

	printf("New poly:\n");
	i=0;
	while(polygon[i]>0){
		printf("polygon[%d]=%d\n", i, polygon[i]);
		i++;
	}
	printf("New polygon length is %d\n", i);

	/*	printf("New poly2:\n");
	j=0;
	while(polygon2[j]>0){
		printf("polygon2[%d]=%d\n", j, polygon[j]);
		j++;
	}
	printf("New polygon2 length is %d\n", j);
	*/
	printf("minlength is %d\n", minlength);

	printf("\nTest.c Complete.\n\n");
}

//#include "bfacf_shrink_links.c"
#include "../generators/bfacf_shrink_v2.c"
