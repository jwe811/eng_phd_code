#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>

void convertUofStoPrintable(char*);

main(void)
{
	printf("\nmain starting\n\n");
	convertUofStoPrintable("Interest_Locally_KnottedL3M1span5.txt");
	printf("\nmain finished\n\n");
}

#include "convertPrintable.c"
