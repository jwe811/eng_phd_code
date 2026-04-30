#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>

int main(void){
	int i,j,k;
	int section_num, nth_outsection;

	int ****start;
	
	int max_sections=10;
	int vec_length = max_sections+1;

	int num_outsections[vec_length];
	for(i=1; i<=max_sections; i++){
		num_outsections[i]=i;
	}


	start = (int ****) malloc(sizeof(int ****)*(max_sections+1));
	for (i = 0; i <=max_sections; i++) {
		start[i] = (int ***) malloc(sizeof(int***)*(num_outsections+1));
		for (j = 0; j <= num_outsections[i]; j++) {
			array[i][j] = (int **) malloc(sizeof(int*)*cols);
		}
	}


	
	for(section_num=1; section_num<=max_sections; section_num++){
		for(nth_outsection=1; nth_outsection<=num_outsections[section_num]; nth_outsection++){
			start[section_num][
		}
	}







	printf("\ntest2.c Complete.\n\n");
	return 0;
}

#include "../utils/int_vecalloc.c"
