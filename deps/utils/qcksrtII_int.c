/***************************************************************************/
/*From Numerical Recipes in C section 8.2*/


#define M 7
#define NSTACK 50
#define FM 7875
#define FA 211
#define FC 1663

/*Here M is the size of subarrays sorted by straight insertion */
/*and NSTACK is the required auxiliary storage.*/

void qcksrtII(unsigned long int n,unsigned long int arr[],unsigned long int brr[])
/*Sorts an array arr[1..n] into ascending numerical order using */
/*the Quicksort algorithm. n is input; arr is replaced on output */
/*by its sorted rearrangement.*/
{
	int l=1,jstack=0,j,ir,iq,i;
	int istack[NSTACK+1];
	long int fx=0L;
	unsigned long int a,b;

	ir=n;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				b=brr[j];/*user mod*/
				for (i=j-1;arr[i]>a && i>0;i--) {
					arr[i+1]=arr[i];
					brr[i+1]=brr[i]; /*user mod*/
				}
				arr[i+1]=a;
				brr[i+1]=b;/*user mod*/
			}
			if (jstack == 0) return;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			i=l;
			j=ir;
			fx=(fx*FA+FC) % FM;
			iq=l+((ir-l+1)*fx)/FM;
			a=arr[iq];
			b=brr[iq];/*user mod*/
			arr[iq]=arr[l];
			brr[iq]=brr[l];/*user mod*/
			for (;;) {
				while (j > 0 && a < arr[j]) j--;
				if (j <= i) {
					arr[i]=a;
					brr[i]=b;/*user mod*/
					break;
				}
				arr[i++]=arr[j];
				brr[i-1]=brr[j];/*user mod*/
				while (a > arr[i] && i <= n) i++;
				if (j <= i) {
					arr[(i=j)]=a;
					brr[(i=j)]=b;/*user mod*/
					break;
				}
				arr[j--]=arr[i];
				brr[j+1]=brr[i];/*user mod*/
			}
			if (ir-i >= i-l) {
				istack[++jstack]=i+1;
				istack[++jstack]=ir;
				ir=i-1;
			} else {
				istack[++jstack]=l;
				istack[++jstack]=i-1;
				l=i+1;
			}
			if (jstack > NSTACK) {
				fprintf(stderr,"NSTACK too small in QCKSRT\n");
			}
		}
	}
}

#undef M
#undef NSTACK
#undef FM
#undef FA
#undef FC

/***************************************************************************/
