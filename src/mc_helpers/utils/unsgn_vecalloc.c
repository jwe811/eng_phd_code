/***************************************************************************/
unsigned int *
unsgn_vecalloc(int low, int high)
{
unsigned int *x;

x=(unsigned int *)calloc((unsigned)(high-low+1),sizeof(unsigned int));
if (x==NULL){
fprintf(stderr,"unable to allocate memory in function unsgn_vecalloc\n");
exit(1);
}
return (x-low);
}
/***************************************************************************/
