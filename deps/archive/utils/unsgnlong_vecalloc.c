/***************************************************************************/
unsigned long int *
unsgnlong_vecalloc(int low, int high)
{
unsigned long int *x;

x=(unsigned long int *)calloc((unsigned)(high-low+1),sizeof(unsigned long int));
if (x==NULL){
fprintf(stderr,"unable to allocate memory in function unsgnlong_vecalloc\n");
exit(1);
}
return (x-low);
}
/***************************************************************************/
