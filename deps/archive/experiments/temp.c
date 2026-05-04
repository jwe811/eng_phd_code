#include <stdio.h>
unsigned long int	getval_mnmx(unsigned long int min, unsigned long int max);
void submain(void);

main(void)
{
submain();
}

void
submain(void)
{
unsigned long int n;
printf("submain\n");
printf("type a number:\n");
n=getval_mnmx(0,10);
printf("you typed %i\n",n);
}


unsigned long int
getval_mnmx(unsigned long int min, unsigned long int max)
/* get unsigned long int from user between the values min and max */
{
	char            input[20];
	unsigned long int x;
	gets(input);
	x = atoi(input);
	if ((min <= x) && (x <= max)) {
		return (x);
	} else {
		return (min);
	}
}
