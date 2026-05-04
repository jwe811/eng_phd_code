/***************************************************************************/
double *
vecalloc(int low, int high)
{
double *x;

x=(double *)calloc((unsigned)(high-low+1),sizeof(double));
if (x==NULL){
fprintf(stderr,"unable to allocate memory in function vecalloc\n");
exit(1);
}
return (x-low);
}
/***************************************************************************/
