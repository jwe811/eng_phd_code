/***************************************************************************/
int *
int_vecalloc(int low, int high)
{
int *x;

x=(int *)calloc((unsigned)(high-low+1),sizeof(int));
if (x==NULL){
fprintf(stderr,"unable to allocate memory in function int_vecalloc\n");
exit(1);
}
return (x-low);
}
/***************************************************************************/
