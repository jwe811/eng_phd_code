/*	This program will create ALL polygons in a LxM tube, with chosen span, that are "patterns of interest". This means their
start and end sections only have two edges. Polygon must go the full span. ALSO-polys must have a min number of edges. 	*/
//32 for 2x1, 26 for 3x1?

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>

#define		bits_16	unsigned short int	/* for occstatebin,hnodebin need one bit for ever vertex   */
						/* on hammer unsigned short int uses 16 bits		   */
						/* occstatebin denotes that the binary representation      */
						/* of the integer "occstatebin" indicates which vertices   */
						/* are occupied (i.e. occstatebin = occupancy state binary */
						/* hnodebin = header node binary			   */
						/* The "header node" is the location of the vertex on      */
						/* the leading edge, as it where, which "walks" through    */
						/* the M by L by 2 lattice to create a potential 2-span	   */
						/* M and L are defined below				   */

/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED	**/
#define	L 2				/* number of horizontal edges							*/
#define	M 1				/* number of vertical edges							*/
#define	totalspan 4		/* span of the polygon that will be generated (must be >=3)	*/

#define 	maxpolys 900000	/* maximum number of polygons in each .txt file	*/





#define	vM (M+1)				/* number of vertices in vertical direction   = M + 1      */
#define	vL (L+1)				/* number of vertices in horizontal direction = L + 1      */

//automatically set maxOrdNum, max_sections, vec_length, depending on the prism size
#if vM*vL<7
#define	maxOrdNum  139			/* Maximum number assigned in num_section.c divided by eight*/
					/* and rounded up plus one to have array[0..number/8] */
					/* V6 => 139  ; V8 => 7778  ; V10 => 699942  */
#define	max_sections 1110		/* Maximum number assigned in num_section.c  */
					/* V6 => 1110 ; V8 => 62216 ; V10 => 5599530 */
					/* 11 bits 16 bits 23 bits */
                                        /* (0,2)->3,8,20,50,126,322,834,2187,5797 in slit cases */
                                        /* 274,725 */
					/* This is a function of the number of vertices !!!insert function here!!!*/

 
#define	vec_length   1111		/* vector dimension for section matrix = max_sections + 1		*/
					/* This is because C assigns arrays of dimension n as [0..(n-1)]	*/
#elif vM*vL<9
#define maxOrdNum	7778
#define max_sections	62216
#define vec_length	62217

#elif vM*vL<12
#define maxOrdNum	699942
#define max_sections	5599530
#define vec_length	5599531
#endif

#define	makeheaderij (1 << (j + i * vL))	/* given i,j yields appropriate hnodebin */

#define		vec_ent double		/*vector entries will be of type defined here*/
#define		mat_ent double		/*matrix entries will be of type defined here*/

#define	newline		printf("\n")	/* Allows the use of newline; in place of printf("\n"); 		  */
#define VALID		1		/* This simply allows statements such as "if(something == VALID) {...}    */
#define NOTVALID 	!VALID		/* This simply allows statements such as "if(something == NOTVALID) {...} */

#define hingeedges (L + 2*M*L + M)


//Calculated these maxendhinge numbers just by running generate_endhinges and counting how many endhinges there were.

#if (M==1 && L==1)
#define maxendhinges 15
#elif (M==1 && L==2)
#define maxendhinges 127
#elif (M==1 && L==3)
#define maxendhinges 1023
#elif (M==2 && L==2)
#define maxendhinges 4095
#endif



/***************************************************************************/
/****************** Structures used in this program ************************/
/***************************************************************************/


struct hinge_span {	/* This data structure is used to store information about two-spans */

	/*unsigned long int 	inorder */		/* no longer needed as it is implicitely defined */
	unsigned long int 	outorder	:26;	/* the section number of the righthand section   */	
	unsigned long int 	edgecount	:6;	/* the number of edges in the hinge		 */
	struct hinge_span 	*nexthinge;		/* pointer to another hinge_span data structure  */
	int				walk[vM*vL*8];	/* this will hold the random walks that the 2span creates (walks are seperated by -222 or -333)	*/
	int				startpoint[2];	//contains the (x,y) coordinates of where the walk is starting (we know z=0)

	unsigned long int	num_con_endhinges_left;	// the number of connecting endhinges that connect to the left of section "inorder"
	unsigned long int*	con_endhinges_left;		// array that holds the connecting endhinge numbers of the endhinges that connect
										// to the left of section "inorder"
	unsigned long int	num_con_endhinges_right;	// the number of connecting endhinges that connect to the right of section "outorder"
	unsigned long int*	con_endhinges_right;	// array that holds the connecting endhinge numbers of the endhinges that connect
										// to the right of section "outorder"

}	*first_hinge_span[vec_length], *current_hinge_span[vec_length]; /* these are global */

	/* first_hinge_span is an array of pointers to hinge_span data structures where each pointer in */
	/* in the array points to the first element of a linked list of hinge_span data structures	 */
	/* The array index implicely defines the "inorder" (or the lefthand section of the two_span)     */

        /* current_hinge_span is a dummy pointer used for traversing the linked lists mentioned above    */
	/* an array of dummy pointers is needed because there is more than one linked list to consider   */


struct bytestruct {

	unsigned char	bit_1:1;
	unsigned char	bit_2:1;
	unsigned char	bit_3:1;
	unsigned char	bit_4:1;
	unsigned char	bit_5:1;
	unsigned char	bit_6:1;
	unsigned char	bit_7:1;
	unsigned char	bit_8:1;

}	onespans[maxOrdNum];  	/* this is global */
				/* if section_num = 8*byteloc + bit_b  where 1 < bit_b < 8	*/
				/* then (onespans[byteloc]).bit_b = 1 if section_num is valid   */
				/* and  (onespans[byteloc]).bit_b = 0 if it is not valid        */
				/* hence there is one bit for each possible section_num         */

struct point {			/* point structure has 2 coordinates, y and x	*/
	unsigned short int x;	/* x ranges from 0-L (horizontal coord)		*/
	unsigned long int y;	/* y ranges from 0-M (vertical coord)		*/
};

struct edge {
	struct point f;		/* first point in the edge	*/
	struct point l;		/* last point in the edge	*/
};

struct endhinge {			/* an endhinge is a linked list of edges			*/
	struct edge edge;		/* the edge that is in this part of the linked list		*/
	struct endhinge* nextendhinge;	/* pointer to the next endhinge structure in the linked list	*/
};

struct endhinge* firstendhinge[maxendhinges];	/* Array of pointers to endhinge structures, where each	*/
						/* pointer in the array points to the first element in	*/
						/* the linked list of endhinge structures		*/
struct endhinge* currentendhinge[maxendhinges];	/* Array of dummy pointers that are used to traverse	*/
						/* the endhinge linked lists.				*/

struct endhinge* firsttemp[maxendhinges];	/* Array of pointers to endhinge structures, where each	*/
						/* pointer in the array points to the first element in	*/
						/* the linked list of endhinge structures		*/
struct endhinge* currenttemp[maxendhinges];	/* Array of dummy pointers that are used to traverse	*/
						/* the endhinge linked lists.				*/

			    


/***************************************************************************/
/****************** end of structures **************************************/
/***************************************************************************/


/***************************************************************************/
/*************** Global variables used in this program *********************/
/***************************************************************************/

int             	ordertemplate[2][vM][vL];	/* The first index is for side (left = 0 ; right = 1) */
							/* the second index is for the location of the vertex from the top */
							/* the third index is for the location of the vertex from the left */
							/* This global array is used as a "template" to store information  */
							/* on the edges in the two sections that make up a two span        */
							/* The array is initiallized to zero and numbers are stored at     */
							/* the appropriate locations (corresponding to given vertices) as  */
							/* the potential two spans are built up algorithmically by the     */
							/* program. The stored numbers represent the order in which the    */
							/* edge is traversed when the two span is part of a polygonal walk */

unsigned long int 	num_2_spans = 1;	/* This variable keeps track of the total number of potential two spans   */

unsigned long int	valid_2_spans = 0;	/* This is used to count the number of two spans which are actually valid */

int				curwalk[vM*vL*8];	/* this will hold the current walks, to be recorded when the 2span is complete. (-222 or -333 separates walks)*/
int				startingx=0, startingy=0;	//starting coordinates for x and y
int				total_end_hinges=0;		/*the total number of possible endhinges	*/
unsigned long int	curindex=0;			//index of which endhinge is being created for ordertemplate (initially zero)
unsigned long int	edgesin_temp = 0;	//will count the number of edges in temp[]
unsigned short int	occupied[vM][vL];	//this will indicate which verticies are occupied (i,j)
int				leftside;				//leftside=1 if testing endhinges on the left; otherwise (testing endhinges on right) leftside=0.

int				pwalksize = vM*vL*(totalspan+1)*4;
int				filetotal=0;
char 			filename[36];			//this will hold the filename of the text file that we want to write to
FILE*			fp;	//file pointer. will point to the file that is to be written to.
int				filenum=1;
int				passedM=M;
int				passedL=L;
int curious=0;
int grand=0;
int grandtotal=0;
int minedges=26;
unsigned long int num_millions=0;
unsigned long int num_uptomill=0;
int smallcounter=0;

/* note: The following are also global: onespans[], first_hinge_span, current_hinge_span */
/* They are discussed in the structures section above

/***************************************************************************/
/******************* end of global variables  ******************************/
/***************************************************************************/

/***************************************************************************/
/******************* functions used in this program  ***********************/
/***************************************************************************/

int			system(const char *string);	/* called by main */
			/* used solely for printing the date and timing the program */

/* conv_to_array	is only called by main and is declared in main */

/****  two_span builder functions  **********************************************************/

void            	enterhinge(bits_16, bits_16, int, int, int, int (*)[3]);
void            	leavehinge(bits_16, bits_16, int, int, int, int (*)[3]);
void            	rowedges(bits_16, bits_16, int, int, int (*)[3]);
void            	coledges(bits_16, bits_16, int, int, int (*)[3]);

/****  two_span validity verifier functions  ******************* *****************************/

unsigned short int	LFlag( int (*)[3]);
unsigned short int	RFlag( int (*)[3]);
unsigned short int	noncrossing(int, int, int, int, int);	/* used in 0_L cases only   */

/****  one_span/section recorder functions  *************************************************/

void				record_1_span(void);
unsigned short int	one_span_byte_bit(unsigned long int, unsigned long int);

void            	recordtemplate( int (*)[3]);
unsigned long int	num_section(int);

/****  Memory allocation functions  ********************************************************/

struct hinge_span	*newhinge(void);
struct endhinge*	newendhinge(void);
unsigned long int	*unsgnlong_vecalloc(int low, int high);
vec_ent			*vecalloc(int low, int high);
mat_ent			**matalloc(int rowlow, int highrow, int lowcol, int highcol);


/**** Added Functions  ********************************************************************/

unsigned long int	choose(unsigned long int, unsigned long int);
unsigned long int 	factorial(unsigned long int n);


void generate_endhinges(void);
void iterate(int k, unsigned short int *edgestoadd, int counter, unsigned short int totedges, unsigned long int* endnum, struct edge* edges, unsigned long int* endhingesadded);
void printendhinge(int k);	//prints endhinge k
void printtempendhinge(int k);	//prints temp[k]

unsigned long int findnumconendhingesleft(void);
unsigned long int findnumconendhingesright(void);
unsigned short int testmatch(int i, int j);

void connect(unsigned long int small, unsigned long int big, unsigned long int cindex);
void goup(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int uindex);
void godown(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int dindex);
void goleft(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int lindex);
void goright(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int rindex);


void printcurwalk(void);
int search(int prevdirec, int curx, int cury, int endhingenumber);
void addtowalk(int printablewalk[], int nextwalk[]);
void recursivegenerate(int curspan, int printablewalk[], struct hinge_span* tspantoadd, int firstsec, int x, int y);
int dupcheck(int x, int y, int dup, int firstsection);
int test2secs(int sectionnum);

/***************************************************************************/
/******************* end of functions used in this program *****************/
/***************************************************************************/


main(void)
{
/**** variable used in main ************************************************/

	/*int	           (*pointordNum)[3]; *//*pointer to an array of int*/
	int				ordNum[3],side = 0;
	int             	i, j, k;	   /*used in for loops */

	unsigned long int	valid_sections=0;
	unsigned long int	byteloc;
	unsigned long int	bitloc;

#if !defined(connectivity_inv)	/* if connectivity_inv is not defined declare it here */

	double			connectivity_inv;

#endif
	double			kappa;
	double			Alpha;		/* Alpha as defined in the notes provided by C. Soteros */
	double			Beta;		/* Beta as defined in the notes provided by C. Soteros */
	double			Ephi_Ephi;
	double			Enum_Enum;	/* Warning: enum is a keyword in C */
	double			LmultR;		/* inner product of left and right eigenvectors */

/**** functions called by main *********************************************/


	clock_t         clock(void);


/**************************************************************************/
/**** start of main program ***********************************************/
/**************************************************************************/

	system("clear"); /* clears the screen */
	system("date");  /* prints the date and time */

	#if defined(TS)
	printf("2 span - 2 span version\n");
	#endif

	#if defined(SS)
	printf("section-section version\n");
	#endif

	printf("\nThe variables are as follows M: %i L: %i\n", M, L);
	#if defined(force)
	printf("\nforce=%f\n", fval);
	#endif


	clock();	/* Start clock to time program */

	ordNum[0] = 1;	/* first edge in section 0 gets numbered 1 (left side of 2-span)	*/
	ordNum[1] = 1;	/* first edge in section 1 gets numbered 1 (right side of 2-span)	*/
	ordNum[2] = 0;	/* to start there are no edges in the 2-span (edges in hinge)		*/

	for (i=1;i<=max_sections;i++){
		current_hinge_span[i] = newhinge();		/* free up space of hinge structure */
		first_hinge_span[i] = current_hinge_span[i];	/* make first_hinge_span[i] point to this first one */
	}
	printf("space created for hinge structures\n");

	//initialize curwalk[] to be all zeros
	for(i=0; i<=vM*vL*(totalspan+2)-1; i++){
		curwalk[i]=0;
	}

	//generate the endhinges
	generate_endhinges();
	//firstendhinge[] (array of endhinges) is now filled with all possible endhinges
	printf("\nALL ENDHINGES GENERATED!!!!! (there were a total of %d endhinges)\n\n", total_end_hinges);

	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			printf("entering hinge at i=%d, j=%d\n", i, j);
			startingx=j;
			startingy=i;
			enterhinge(makeheaderij, makeheaderij, i, j, side, &ordNum);
			printf("finished at i=%d, j=%d\n", i, j);
		}
	}

/*************************************************************************************************/
	printf("\nDone creating all endhinges and 2-spans in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	newline; /* results are printed to the terminal */
	printf("There are a total of %9i potential two-span patterns.\n", --num_2_spans);
	printf("There are a total of %9i   valid   two-span patterns.\n", valid_2_spans);
	for (byteloc = 0; byteloc < maxOrdNum; byteloc++) {
		for (bitloc = 0; bitloc < 8; bitloc++) {
			if (one_span_byte_bit(byteloc,bitloc)==1) {valid_sections++;}
		}
	}
	printf("There are a total of %9i   valid   section orderings.\n\n", valid_sections);
/************************************************************************************************/



	//TIME TO CREATE SOME POLYGONS! YAYYYYY  ^^;;
	printf("TIME TO CREATE SOME POLYGONS!!!!!\n");

	sprintf(filename, "CreatorAll_InterestPKPolysL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);
	fp = fopen(filename, "w");	//create or overwrite the file "filename"

	if(fp != NULL){
		fprintf(fp, "UofS\n");	//first line in file is always "UofS"
		printf("printed UofS in file '%s'\n", filename);
	}
	else{
		printf("file pointer is pointing to NULL\n");
		exit(1);
	}



	int printablewalk[pwalksize];	//this will store the walk that is to be printed
	for(i=0; i<=pwalksize-1; i++){	//initialize printable walk to be all zeros
		printablewalk[i]=0;
	}

	int x, y, z=0;							//this will be the coordinates of the starting point
	int curx, cury, curendhinge;				//going to use these variables to aid in adding the endhinge to printablewalk.
	int newdirec, prevdirec, direcsadded;
	int temp, temp2;
	int counter;
	int dummy;
	int firstsection;
	int nextsection;
	int memorywalk[pwalksize];	//this will remember what printablewalk is before recursive generate is ran
	int dup;



	//loop through all 2-spans(next 2 lines), then check if it is a valid starting 2-span
	for(firstsection=1; firstsection<=max_sections; firstsection++){
		if(test2secs(firstsection)==1){	//ONLY WANT POLYGONS WITH 2 SECTIONS IN THE FIRST SECTION
			current_hinge_span[firstsection] = first_hinge_span[firstsection];	//get current pointing to first
			while((*current_hinge_span[firstsection]).nexthinge != NULL){
				current_hinge_span[firstsection] = (*current_hinge_span[firstsection]).nexthinge;	//current is now pointing at a possible starting 2-span
				//		printf("secondsection=%d\n", current_hinge_span[firstsection]->outorder);
				if((*current_hinge_span[firstsection]).num_con_endhinges_left > 0){	//this tests to make sure we have a VALID starting 2-span (has a connecting endhinge on the left)
					//we have a valid starting 2-span

					//printf("WE HAVE A VALID STARTING 2-SPAN IT STARTS HERE! ");
					//printf("LEFT ENDHINGES = %d.\n", (*current_hinge_span[firstsection]).num_con_endhinges_left);
					x = (*current_hinge_span[firstsection]).startpoint[0];	//record x-coord of startpoint
					y = (*current_hinge_span[firstsection]).startpoint[1];	//record y-coord of startpoint
					//if starting point isn't (0,0), need to check if it's a duplicate.
					//THIS IS THE DUPLICATE TEST
					dup=0;
					if(x>0 || y> 0){
					  //printf("no longer starting at (0,0), bye bye!\n");
						//exit(1);
	//					printf("starting dup test. startpoint=(%d,%d). \n", x, y);
						dup = 0; //assume initially it is NOT a duplicate
						dup = dupcheck(x, y, dup, firstsection);
						if(dup==1){
	//						printf("found a duplicate. not printing it.\n");
						}
					}

					if(dup==0){

						//put walk[] into printablewalk[]
						i=0;
						while((*current_hinge_span[firstsection]).walk[i] != 0){
							printablewalk[i] = (*current_hinge_span[firstsection]).walk[i];
							i++;
						}

						//so we have the startpoint (x and y) and we have the start of printablewalk[].
	/*					printf("STARTPOINT = (%d,%d,%d). printablewalk[] thusfar: ", x, y, z);
						i=0;
						while(printablewalk[i] != 0){
							printf("%d, ", printablewalk[i]);
							i++;
						}

						printf("\n");
	*/



						//next, add appropriate steps to printablewalk[] for all possible endhinges
						for(curendhinge=1; curendhinge <= (*current_hinge_span[firstsection]).num_con_endhinges_left; curendhinge++){	//for each endhinge that connects

							//THIS IS A KEYPOINT. START OF THE HUGE LOOP.

							if(curendhinge>1){	//if we're looking at a 2nd or 3rd or so on endhinge, then we must first reset printablewalk.
								//put walk[] into printablewalk[]
								i=0;
								while((*current_hinge_span[firstsection]).walk[i] != 0){
									printablewalk[i] = (*current_hinge_span[firstsection]).walk[i];
									i++;
								}
							}

		//					printendhinge((*current_hinge_span[firstsection]).con_endhinges_left[curendhinge-1]);
							i=0;
		//					printf("starting from beginning of printable walk.\n");
							while(printablewalk[i]!=0){
								curx=x;
								cury=y;
								while(printablewalk[i] != -333){
									if(printablewalk[i]==1 || printablewalk[i]==-1){
										curx++;
									}
									else if(printablewalk[i]==2 || printablewalk[i]==-2){
										curx--;
									}
									else if(printablewalk[i]==3 || printablewalk[i]==-3){
										cury++;
									}
									else if(printablewalk[i]==4 || printablewalk[i]==-4){
										cury--;
									}
									i++;
								}

								//at this point, printablewalk[i] is -333.
								if(printablewalk[i]==-333){	//if walk has left to the left, need to insert part of the endhinge
		//							printf("need to insert for endhinge\n");
									prevdirec=0;	//will be zero if this is the first insertion
									direcsadded=0;	//initialized to zero appropriately
									newdirec=1;	//initialized to 1 so it will run atleast once
							
									while(newdirec!=0){	//while we previously added a direction,
										//search for edge in endhinge that has curx and cury as a point in it
										//newdirec is the direction that is to be added. newdirec will be zero if there are
										//no more directions to be added.
		//								printf("calling search for the next direction (prevdirec=%d, curx=%d, cury=%d, endhingenumber=%d)\n", prevdirec, curx, cury, (*current_hinge_span[firstsection]).con_endhinges_left[curendhinge-1]);
										newdirec = search(prevdirec, curx, cury, (*current_hinge_span[firstsection]).con_endhinges_left[curendhinge-1]);
		//								printf("the search came up with direction %d\n", newdirec);
										//make changes to curx and cury.
										if(newdirec==1){
											curx++;
										}
										else if(newdirec==2){
											curx--;
										}
										else if(newdirec==3){
											cury++;
										}
										else if(newdirec==4){
											cury--;
										}
										//add newdirec into printablewalk if it isn't 0.
										if(newdirec!=0){
											j=i;
											temp=printablewalk[j];
											printablewalk[j]=newdirec;
											j++;
											while(temp != 0){
												temp2=printablewalk[j];
												printablewalk[j] = temp;
												temp=temp2;
												j++;
											}
											//newdirec is now inserted into printablewalk in spot i and everything else
											//is shifted to the right.
											direcsadded++;
											prevdirec=newdirec;
										}
										i++;
									}
									i--;
									//printablewalk[i] should be -333 at this point
		//							printf("!!!!!!!!!!!printablewalk[%d]=%d,\n", i, printablewalk[i]);
									//if printablewalk[i+1] is zero, then we are at the end. replace the -333 by -111.
									//elseif we are not at the end, remove the -333 and all negative entries after it
									if(printablewalk[i+1]==0){
		//								printf("we are at the end. replace -333 with -111\n");
										printablewalk[i]=-111;
										i++;

										//prints "finished" printable walk
		//								printf("finished adding left endhinge, ");
		//								printf("STARTPOINT = (%d,%d,%d). printablewalk[]: ", x, y, z);

		/*								dummy=0;
										while(printablewalk[dummy] != 0){
											printf("%d, ", printablewalk[dummy]);
											dummy++;
										}
										printf("\n");
		*/

									}
									else{//remove -333 and all negative numbers afterwards (until a positive)
										j=i;
										counter=0;	//counts the number of entries that will have to be removed after the -333
										while(printablewalk[j+1]<0){
											counter++;
											j++;
										}
		//								printf("counter=%d\n", counter);
										//remove the -333 and the next "counter" entries.
										while(printablewalk[i+counter+1]!=0){
											printablewalk[i]=printablewalk[i+counter+1];
											i++;
										}

										//fill everything afterwards with zeros
										j=i+counter+1;
										while(i<=j){
											printablewalk[i]=0;
											i++;
										}
								
										//since this isnt the end of the walk, restart loop with i=0 to look for more needed insertions
										i=0;

									}
								}
							}
							//ALL STEPS FROM THE LEFT ENDHINGE HAVE BEEN ADDED TO THE CURRENT WALK THAT IS BEING LOOKED AT.
							//NOW NEED TO ADD INTERMEDIATE 2-spans (between the first and end 2-spans)





							//DO THE CHECK FOR DUPILICATES HERE, SINCE THE LEFT ENDHINGE AND THE FIRST SECTION ARE NOW CREATED AND STORED IN PRINTABLEWALK[].
					












							//recursively do this? need to keep track of depth (span).
					
							/*add next 2-span (add appropriate steps into printable walk.)
								-what do we need to know to do this? loop through all "outorders", and this way we can look at all
								2-spans (each outorder has a set of 2-spans associated with it)
								-probably make a function that takes in two walks, (one for the left 2-span (printable walk) and one for the right).
								Function will add the appropriate steps from the right walk to the left walk.

							*/

							//copy printablewalk into memorywalk
							for(i=0; i<=pwalksize-1; i++){
								memorywalk[i] = printablewalk[i];
							}

							nextsection=current_hinge_span[firstsection]->outorder;
							struct hinge_span* nexttspan;
							nexttspan = first_hinge_span[nextsection];	//point nextsection's current at first
							while(nexttspan->nexthinge != NULL){		//loop through all of nextsection's outorders.
								nexttspan = nexttspan->nexthinge;
		//						printf("\n\nfirst section=%d, second section=%d, third section =%d\n", firstsection, nextsection, nexttspan->outorder);
		//						printf("thirdsection=%d\n", nexttspan->outorder);
								recursivegenerate(2, printablewalk, nexttspan, nextsection, x, y);	//all polygons now generated with this starting 2-span.
								//reset printablewalk to memorywalk
								for(i=0; i<=pwalksize-1; i++){
									printablewalk[i]=memorywalk[i];
								}
							}



							//reset printable walk to all zeros
							for(i=0; i<=pwalksize-1; i++){	//inialize printable walk to be all zeros
								printablewalk[i]=0;
							}
						}
					}
				}
			}
		}
	}








	fprintf(fp, "-999\n");

	printf("L=%d, M=%d, span=%d, filetotal==%d, filenum=%d\n", L, M, totalspan, filetotal, filenum);
	printf("Only recording polygons with >=%d edges\n", minedges);
	if(filenum==1){
		printf("Total Number of Polys OF INTEREST = %d = %f = %e\n", filetotal + smallcounter, (float)filetotal + smallcounter, (float)filetotal + smallcounter);
		printf("Number of Polys OF INTEREST with >=%d edges = %d = %f = %e!\n", minedges, filetotal, (float)filetotal, (float)filetotal);
	
	}
	else{
		printf("Total Number of Polys OF INTEREST= %d = %f = %e\n", filetotal + ((filenum-1)*(maxpolys+1)) + smallcounter, (float)filetotal + ((filenum-1)*(maxpolys+1)) + smallcounter, (float)filetotal + ((filenum-1)*(maxpolys+1)) + smallcounter);
		printf("Number of Polys OF INTEREST with >=%d edges = %d = %f = %e!\n", minedges, filetotal + ((filenum-1)*(maxpolys+1)), (float)filetotal + ((filenum-1)*(maxpolys+1)), (float)filetotal + ((filenum-1)*(maxpolys+1)));
	}

	printf("grand=%d\n", grand);

	printf("\nDone generating all polygons of interest with >= %d edges in %f seconds.\n\n", minedges, (double)clock()/(double)CLOCKS_PER_SEC);
	fclose(fp);

	printf("curious=%d (number of fake times thought we were at a dup point)\n", curious);

	printf("COUNTING BY MILLIONS: Number of millions=%d, with %d as a remainder.\n", num_millions, num_uptomill);

	system("date");  /* prints the date and time */
	printf("\nProgram Complete.\n\n");


}				/* end of main */


/***************************************************************************/

void
enterhinge(bits_16 occstatebin, bits_16 hnodebin, int i, int j, int side, int (*pointordNum)[3])
{

	//add direction 5 or 6 to curwalk
	int index=0;
	while(curwalk[index]!=0){
		index++;
	}
	//if side==0, then we are entering from the left, so we are in direction posz (5)
	//else, we are on side==1, and entering from the right (negz, 6)
	if(side==0){
		curwalk[index]=5;
	}
	else{
		curwalk[index]=6;
	}

	ordertemplate[side][i][j] = (*pointordNum)[side];
		/* section edge through which the SAW enters the hinge is given a number */
		/* this number represents the order in which it occurs in the section    */

	(*pointordNum)[side]++;
		/* the edge number through which the hinge is entered */
		/* is incremented so that the next edge gets a higher number */

	/* the SAW can now do one of three things */
	/* it may leave the hinge in which case leavehinge(...) is called */
	/* it may explore edges in the horizontal direction => rowedges (...) is called */
	/* it may explore edges in the vertical direction => coleges (...) is called    */
	/* in each case the currrent state of the hinge is passed via the parameters    */


	leavehinge(occstatebin, hnodebin, i, j, side ^ 1, pointordNum);
		/* the SAW leaves the hinge out the opposite side, hence (side ^ 1) is passed */
		/* if side = 1 then side ^ 1 = 0 */
		/* if side = 0 then side ^ 1 = 1 */
	rowedges(occstatebin, hnodebin, i, j, pointordNum);
		/* horizontal edges are explored */
	coledges(occstatebin, hnodebin, i, j, pointordNum);
		/* vertical edges are explored */

	(*pointordNum)[side]--;
		/* the edge number is decremented since, control is being passed back to  */
		/* the function that called enterhinge. That is to say, the program will  */
		/* now explore alternate ways of entering the hinge (i.e. from different  */
		/* locations). In essense the intial incremention above is nullified so	  */
		/* when enterhinge is called again the entering edge is numbered properly */
		
	ordertemplate[side][i][j] = 0;
		/* the edge that was recorded must now be erased as this edge may not exist */
		/* as part of the SAW when alternate possibilities are explored. In the case */
		/* that it does not exist it needs to be numbered as "0" (zero).             */

	//delete last entry in curwalk
	index=0;
	while(curwalk[index+1] != 0){
		index++;
	}
	curwalk[index]=0;

	return;

}				/* end of function enterhinge */

/***************************************************************************/

void
leavehinge(bits_16 occstatebin, bits_16 hnodebin, int i, int j, int side, int (*pointordNum)[3])
/* called by: rowedges,coledges,enterhinge */
{

	//add direction 5 or 6 to curwalk
	int index=0;
	while(curwalk[index]!=0){
		index++;
	}
	//if side==0, then we are exiting to the left, so we are in direction negz (6) (also add -333 for exiting left)
	//else, we are exiting to side==1 (right) (posz, 5) (also add -222 for exiting right)
	if(side==0){
		curwalk[index]=6;
		curwalk[index+1]=-333;
	}
	else{
		curwalk[index]=5;
		curwalk[index+1]=-222;
	}


	int       ii;		/*vertex number in vertical direction*/
	int		jj;		/*vertex number in horizontal direction */
	int		shift = 0;	/*keeps track of how far the headernode is shifted */
					/*from the 0,0 position going left to right first */
					/*then down as each row of vertices is exhausted  */
	int loopi, loopj, looper;
	int fullhinge=1;
	int totadded;



	ordertemplate[side][i][j] = (*pointordNum)[side];	/* record the section
								 * edge (i.e. the
								 * leaving edge) */
	(*pointordNum)[side]++;
		/* the edge number is incremented so that the next edge gets a higher number */

	if (side == 0 && (*pointordNum)[1] > 1) {
 	/* if the leaving edge is on the left and there is at least one edge on the right section */
	/* this may be a valid two span */
		if (LFlag(pointordNum) && RFlag(pointordNum)){
		/* if the two span connects to phi on the left and on the right then do the following*/

			valid_2_spans++; /*This is a valid 2 span so it is counted*/
			record_1_span(); /*The component one spans are therefore valid so need to be recorded*/
			recordtemplate(pointordNum); /*the information contained in ordertemplate is recorded*/

/*				printf("recorded\n");
			for(loopi=0; loopi<=M; loopi++){
				for(loopj=0; loopj<=L; loopj++){
					printf("ordertemplate[0][%d][%d] = %d\n", loopi, loopj, ordertemplate[0][loopi][loopj]);
				}
			}
			for(loopi=0; loopi<=M; loopi++){
				for(loopj=0; loopj<=L; loopj++){
					printf("ordertemplate[1][%d][%d] = %d\n", loopi, loopj, ordertemplate[1][loopi][loopj]);
				}
			}
*/
		}
		num_2_spans++; /*Keep track of all 2 spans, valid and non-valid, just for interest sake*/
	}
	/* since a call to leavehinge means the "headernode" is no longer in the two span     */
	/* we can imagine it doing a walk outside the two span and eventually coming back     */
	/* and re-entering the two span into the hinge. The following amounts to exploring    */
	/* the possible ways the walk can get back into the hinge			      */
	/* Note: there may actually exist a walk via which the header node re-enters at       */
	/* a given point, but all must be considered. Eventually non-valid ones are discarded */
	/* Program looks at each vertex in the hinge through the ii and jj for loops          */
	for (ii = 0; ii <= M; ii++) {
	/* these are the "vertical" locations of the vertices being considered */
		for (jj = 0; jj <= L; jj++) {
		/* these are the horizontal locations of the vertices being considered */
			if ((occstatebin >> (shift) & 1) == 0) {
			/* if the vertex is not occupied then explore entering the hinge at this vertex */
			/* this maintains the self avoiding nature of the walk */

				//need to add steps into curwalk (as negative numbers since the walk is outside the hinge)
				//exited hinge at (j,i), re-entering at (jj,ii).
				index=0;
				while(curwalk[index]!=0){
					index++;
				}

				totadded=0;	//will keep track of how many steps outside of the hinge are added, so we can delete them later

				//add horizontal steps (-1's or -2's)
				if(jj>j){		//if reenterj>originalj, add -1's
					for(loopj=j+1; loopj<=jj; loopj++){
						curwalk[index]=-1;
						index++;
						totadded++;
					}

				}
				else if(jj<j){	//add -2's
					for(loopj=jj+1; loopj<=j; loopj++){
						curwalk[index]=-2;
						index++;
						totadded++;
					}
				}

				//add vertical steps (-3's or -4's)
				if(ii>i){		//if reenteri>originali, add -3's
					for(loopi=i+1; loopi<=ii; loopi++){
						curwalk[index]=-3;
						index++;
						totadded++;
					}
				}
				else if(ii<i){	//add -4's
					for(loopi=ii+1; loopi<=i; loopi++){
						curwalk[index]=-4;
						index++;
						totadded++;
					}
				}



				enterhinge(occstatebin + (1 << (shift)), 1 << (shift), ii, jj, side, pointordNum);
				
				//delete the last "totadded" from curwalk
				index=0;
				while(curwalk[index]!=0){
					index++;
				}
				for(looper=1; looper<=totadded; looper++){
					curwalk[index-looper]=0;
				}

			}
			shift++;
			/* the order in which the ii and jj are incremented means that at each incrementation */
			/* of either ii or jj the binary representation of the header node is shifted by one  */
			/* that is 0001 shifts to 0010 then 0100 then 1000				      */
			/* Warning: when printed for the user 0001 is printed as 1000				*/
			/*                                    0010 is printed as 0100		      		*/
			/*				      0100 is printed as 0010				*/
			/*				      1000 is printed as 0001				*/
			/*											*/
			/* Here is another example: 								*/
			/*				suppose that hnodebin = 16 and M = 1 and L = 3		*/
			/*											*/
			/* then the header node is located at (ii,jj) = (1,0) = 				*/

			/*		0000									*/
			/*		1000									*/

			/* As the ii and jj loops are executed here is what happens				*/
			/*											*/
			/*	1000 ->	0100 ->	0010 ->	0001 ->	0000 ->	0000 ->	0000 ->	0000			*/
			/*	0000	0000	0000	0000	1000	0100	0010	0001			*/

			/* The respective values of shift are 0,1,2,3,4,5,6, and 7				*/

			/* The advantage to doing this is that it helped in remembering which way was "up"	*/
			/* "down" "left" or "right" (not that these concepts are important to a computer)	*/
			/* In short: a positive change in ii means a move "down" whereas negative is "up"	*/
			/*	     a positive dhange in jj means a move "right" whereas negative is "left"	*/

			/* Note: I chose this convention to correspond with the way in which matrix enteries	*/
			/*	 are indexed (aside from the fact that ii, and jj start from zero)		*/			
		}
	}
	(*pointordNum)[side]--;
		/* the edge number is decremented since, control is being passed back to  */
		/* the function that called leavehinge. That is to say, the program will  */
		/* now explore alternate ways of leaving the hinge (i.e. via different    */
		/* paths. this may be as simple as going out the other side). In essense  */
		/* the intial incremention above is nullified so when leavehinge is 	  */
		/* called again the leaving edge is numbered properly  			  */
		
	ordertemplate[side][i][j] = 0;
		/* the edge that was recorded must now be erased as this edge may not exist  */
		/* as part of the SAW when alternate possibilities are explored. In the case */
		/* that it does not exist it needs to be numbered as "0" (zero).             */

	//delete last 2 entries in curwalk
	index=0;
	while(curwalk[index+2] != 0){
		index++;
	}
	curwalk[index]=0;
	curwalk[index+1]=0;


	return;
}
/***************************************************************************/

void
rowedges(bits_16 occstatebin, bits_16 hnodebin, int i, int j, int (*pointordNum)[3])
{
	int             hmove;	/* -1 corresponds to moving "left"  */
				/* +1 corresponds to moving "right" */
	int		jmoved;	/* once the header node moves it receives a new j co-ordinate */
				/* for convienience it is called jmoved */

	/* the following uses a conditional expression " expr1 ? expr2 : expr3 " */
	/* "The expression expr1 is evaluated first. If it is non-zero (true), then the */
	/* expression expr2 is evaluated, and that is the value of the conditional expression. */
	/* Otherwise expr3 is evaluated, and that is the value." from K & R pg.51 */

	/* the expression (occstatebin >> j + i * (vL)) & 1  evaluates to 0 if the i,j vertex is empty */
	/* otherwise it evaluates to 1 */


	for (hmove = ((j == 0) ? 0 : (occstatebin >> (j - 1 + i * (vL)) & 1) - 1 ); \
	     hmove <= ((j == L) ? 0 : 1 - (occstatebin >> (j + 1 + i * (vL)) & 1) ); \
	     hmove++) {

	/* if j = 0 start loop at 0   That is, a move to the left  (-1 move) is disallowed      */
	/* if the vertex to the left is empty then (occstatebin >> (j - 1 + i * (vL)) & 1) - 1 evaluates to -1 */
	/* so that a left move is considered, otherwise the loop starts at 0			*/

	/* if j = L finish loop at 0  That is, a move to the right (+1 move) is disallowed      */
	/* if the vertex to the right is empty then 1 - (occstatebin >> (j + 1 + i * (vL)) & 1) evaluates to +1 */
	/* so that a right move is considered, otherwise the loop ends at 0			*/

		if (hmove != 0) { /* if a move is made then do this stuff */

			//add direction 1 or 2 to curwalk
			int index=0;
			while(curwalk[index]!=0){
				index++;
			}
			//if hmove==-1, add direction negx (2)
			//elseif hmove=1, add direction posx (1)
			if(hmove==-1){
				curwalk[index]=2;
			}
			else if(hmove==1){
				curwalk[index]=1;
			}


			bits_16         occstatebinbuff, hnodebinbuff;	/* allocate space for these buffers */
			occstatebinbuff = occstatebin + (hnodebinbuff = ((hmove == 1) ? hnodebin << 1 : hnodebin >> 1));
				/* this updates the occupancy state */
				/* if hmove = 1 then hnodebin << 1 does a left bit shift */
				/* which corresponds to moving to the right */
				/* if hmove =-1 then hnodebin >> 1 is a right bit shift */
				/* which corresponds to moving to the left */
			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			leavehinge(occstatebinbuff, hnodebinbuff, i, jmoved = (j + hmove), 0, pointordNum);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(occstatebinbuff, hnodebinbuff, i, jmoved, 1, pointordNum);
				/* exit the hinge on side 1 with j moved*/
			rowedges(occstatebinbuff, hnodebinbuff, i, jmoved, pointordNum);
				/* explore horizontal moves within the hinge */
			coledges(occstatebinbuff, hnodebinbuff, i, jmoved, pointordNum);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */

			//delete last entry in curwalk
			index=0;
			while(curwalk[index+1] != 0){
				index++;
			}
			curwalk[index]=0;


		}		/* end of if hmove <> 0 */
	}			/* end of for loop for horizontal moves */
	return;
}				/* end of function rowedges */

/***************************************************************************/

void
coledges(bits_16 occstatebin, bits_16 hnodebin, int i, int j, int (*pointordNum)[3])
/* called by rowedges,coledges */
{
	int             vmove;	/* -1 corresponds to moving "up"   */
				/* +1 corresponds to moving "down" */
	int		imoved; /* once the header node moves it receives a new i co-ordinate */
				/* for convienience it is called imoved */

	/* the following uses a conditional expression " expr1 ? expr2 : expr3 " */
	/* "The expression expr1 is evaluated first. If it is non-zero (true), then the */
	/* expression exprs2 is evaluated, and that is the value of the conditional expression. */
	/* Otherwise expr3 is evaluated, and that is the value." from K & R pg.51 */

	/* the expression (occstatebin >> j + i * (vL)) & 1  evaluates to 0 if the i,j vertex is empty */
	/* otherwise it evaluates to 1 */

	for (vmove = ((i == 0) ? 0 : (occstatebin >> (j + (i - 1) * (vL)) & 1) - 1 ); \
	     vmove <= ((i == M) ? 0 : 1 - (occstatebin >> (j + (i + 1) * (vL)) & 1) ); \
	     vmove++) {

	/* if i = 0 start loop at 0   That is, a move to upwards  (-1 move) is disallowed      */
	/* if the vertex to the top is empty then (occstatebin >> (j + (i - 1) * (vL)) & 1) - 1 evaluates to -1 */
	/* so that an upward move is considered, otherwise the loop starts at 0			*/

	/* if i = M finish loop at 0  That is, a move to the bottom (+1 move) is disallowed      */
	/* if the vertex to the bottom is empty then  1 - (occstatebin >> (j + (i + 1) * (vL)) & 1) evaluates to +1 */
	/* so that a downward move is considered, otherwise the loop ends at 0			*/

		if (vmove != 0) {

			//add direction 3 or 4 to curwalk
			int index=0;
			while(curwalk[index]!=0){
				index++;
			}
			//if vmove==-1, add direction negy (4)
			//elseif vmove=1, add direction posy (3)
			if(vmove==-1){
				curwalk[index]=4;
			}
			else if(vmove==1){
				curwalk[index]=3;
			}


			bits_16         occstatebinbuff, hnodebinbuff; /* allocate space for these buffers */
			occstatebinbuff = occstatebin + (hnodebinbuff = ((vmove == 1) ? hnodebin << (vL) : hnodebin >> (vL)));
				/* this updates the occupancy state */
				/* if vmove = 1 then hnodebin << (vL) does a left bit shift by vL */
				/* which corresponds to moving to the down */
				/* if vmove =-1 then hnodebin >> (vL) is a right bit shift by vL */
				/* which corresponds to moving up */
			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			leavehinge(occstatebinbuff, hnodebinbuff, imoved = (i + vmove), j, 0, pointordNum);
				/* exit the hinge on side 0 with i moved */
			leavehinge(occstatebinbuff, hnodebinbuff, imoved, j, 1, pointordNum);
				/* exit the hinge on side 1 with i moved */
			rowedges(occstatebinbuff, hnodebinbuff, imoved, j, pointordNum);
				/* explore horizontal moves within the hinge */
			coledges(occstatebinbuff, hnodebinbuff, imoved, j, pointordNum);
				/* explore vertical moves withing the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */

			//delete last entry in curwalk
			index=0;
			while(curwalk[index+1] != 0){
				index++;
			}
			curwalk[index]=0;

		}		/* end of if vmove <> 0 */
	}			/* end of for loop for vertical moves*/
	return;
}				/* end of function coledges */

/***************************************************************************/

/***************************************************************************/
/* This function is simply for recording which one spans are valid */
/* where one bit is used for each possible section number */
/* the byte number is determined by the quotient of the section number divided by 8 */
/* the bit number is determined by the remainder of the section number divided by 8 */

void
record_1_span(void)
{
unsigned long int	byteloc,bitloc;
byteloc = (num_section(0)/8);
bitloc  = (num_section(0)%8);
switch(bitloc) {
		case 0:{(onespans[byteloc]).bit_1 = 1;break;}
		case 1:{(onespans[byteloc]).bit_2 = 1;break;}
		case 2:{(onespans[byteloc]).bit_3 = 1;break;}
		case 3:{(onespans[byteloc]).bit_4 = 1;break;}
		case 4:{(onespans[byteloc]).bit_5 = 1;break;}
		case 5:{(onespans[byteloc]).bit_6 = 1;break;}
		case 6:{(onespans[byteloc]).bit_7 = 1;break;}
		case 7:{(onespans[byteloc]).bit_8 = 1;break;}

		default:{	printf("Warning: unexpected switch to default in record_1_span");
				break;
				/* This should never actually execute if the program is error free */
			}
		}
return;
}
/***************************************************************************/
/* This returns 1 if (onespan[byteloc]).bit_bitloc is a valid one span     */
/*  and returns 0 if (onespan[byteloc]).bit_bitloc is not a valid one span */
/* provided this function is called after 	*/
/* the valid one spans have been recorded using */
/* the function record_1_span 			*/ 

unsigned short int
one_span_byte_bit(unsigned long int byteloc, unsigned long int bitloc)
{
switch(bitloc) {
		case 0:{return((onespans[byteloc]).bit_1);break;}
		case 1:{return((onespans[byteloc]).bit_2);break;}
		case 2:{return((onespans[byteloc]).bit_3);break;}
		case 3:{return((onespans[byteloc]).bit_4);break;}
		case 4:{return((onespans[byteloc]).bit_5);break;}
		case 5:{return((onespans[byteloc]).bit_6);break;}
		case 6:{return((onespans[byteloc]).bit_7);break;}
		case 7:{return((onespans[byteloc]).bit_8);break;}
		default:{	printf("Warning: unexpected switch to default in one_span_byte_bit");
				return(0);
				break;
				/* This should never actually execute if the program is error free */
			}
		}
}
/***************************************************************************/

/***************************************************************************/
void recordtemplate( int (*pointordNum)[3])
/* This function records the pertinent information about a two-span */
/* which has been verified as being valid into a linked list */
{
	int i, j;
	unsigned long int	inNum = num_section(0);
	unsigned short int	Num_edges= (*pointordNum)[0] + (*pointordNum)[2] - 1;
				/*number of edges = edges in first section + edges in hinge */
				/*the correction (-1) is necessary since (*pointordNum)[0] */
				/*is the number that the next edge in the first section will */
				/*get */

	(*current_hinge_span[inNum]).nexthinge = newhinge(); /* add a newhinge to the linked list */
	current_hinge_span[inNum] = (*current_hinge_span[inNum]).nexthinge; /* let current_hinge_span[inNum] point to newhinge */
	(*current_hinge_span[inNum]).outorder = num_section(1);	/* record the ordering on the out side */
	(*current_hinge_span[inNum]).edgecount = Num_edges;	/* record the number of edges in first section and in hinge */



//	printf("inNum=%d, outNum=%d. ", inNum, (*current_hinge_span[inNum]).outorder);
//	printf("recording curwalk[]: ");
	//record startpoint (then print starting point)
	(*current_hinge_span[inNum]).startpoint[0]=startingx;
	(*current_hinge_span[inNum]).startpoint[1]=startingy;
//	printf("(%d, %d, 0), ", (*current_hinge_span[inNum]).startpoint[0], (*current_hinge_span[inNum]).startpoint[1]);

	//copy curwalk into the hinge_span structure (and print the walk)
	i=0;
	while(curwalk[i]!=0){
		(*current_hinge_span[inNum]).walk[i] = curwalk[i];
//		printf("%d, ", curwalk[i]);
		i++;
	}


	//calculate num_con_endhinges_left
	//find which endhinges connect to left side of ordertemplate[0] (left side)
	//record their "number" in con_edges_left[num_con_endhinges_left]

	//Then do the same for endhinges that connect to the right side of the section (*current_hinge_span[inNum]).outorder
	//store this num_con_endhinges_right
	//record their "number" in con_edges_right[num_con_endhinges_right]


	unsigned long int num_con_endhinges_left = findnumconendhingesleft();	//finds the number of endhinges that ordertemplate[0] connects to on the left
//	printf("there are %d endhinges that connect to that left side of ordertemplate\n", num_con_endhinges_left);

	if(num_con_endhinges_left==0){
		current_hinge_span[inNum]->num_con_endhinges_left = num_con_endhinges_left;	//number of endhinges that connect to the left of section inNum
	}
	else{
		unsigned long int con_edges_left[num_con_endhinges_left];		//will contain the endhinge "numbers" that connect to the left of section inNum
		unsigned long int index=0;		//index of con_edges_left
		unsigned short int matchresult=0;	//0 if no match, 1 if match
	
		for(i=0; i<=num_con_endhinges_left-1; i++){	//loop through all temp[i=0...num_con_endhinges_left-1]
			//find the corresponding endhinge that contains the same edges as temp[i], then add the endhinge number to con_edges_left[index]
			for(j=0; j<=total_end_hinges-1; j++){	//loop through all endhinges
				matchresult = testmatch(i,j);	//tests if temp[i] matches endhinge[j] (linked lists, doesn't have to be in the same order)
				if(matchresult==1){	//if temp[i] matches endhinge[j]
					con_edges_left[index] = j;
					index++;
					break;
				}
			}
		}
	
		//so con_edges_left contains the indicies of endhinge[] that connect to the left of ordertemplate[0]
	
		current_hinge_span[inNum]->num_con_endhinges_left = num_con_endhinges_left;	//number of endhinges that connect to the left of section inNum
		current_hinge_span[inNum]->con_endhinges_left = unsgnlong_vecalloc(0, num_con_endhinges_left-1);	//list of endhinge "numbers" that connect to the left of section inNum
		for(i=0; i<=num_con_endhinges_left-1; i++){
			current_hinge_span[inNum]->con_endhinges_left[i] = con_edges_left[i];
		}
	}

	curindex=0;	//reset curindex
	edgesin_temp=0;	//reset edgesin_temp


	//free firsttemp[] and all of it's linked components
	for(i=0; i<=maxendhinges-1; i++){
		currenttemp[i] = firsttemp[i];
		while(currenttemp[i]->nextendhinge != NULL){
			currenttemp[i] = currenttemp[i]->nextendhinge;
			free(firsttemp[i]);
			firsttemp[i] = currenttemp[i];
		}
		free(firsttemp[i]);
	}


//	printf("# of endhinges that connect to the left: %d. ", current_hinge_span[inNum]->num_con_endhinges_left);
	if(current_hinge_span[inNum]->num_con_endhinges_left > 0){
//		printf("Endhinge numbers: ");
		for(i=0; i<=num_con_endhinges_left-1; i++){
//			printf("%d, ", current_hinge_span[inNum]->con_endhinges_left[i]);
//			printendhinge(current_hinge_span[inNum]->con_endhinges_left[i]);
		}
	}



//	DO THE SAME THING, BUT NOW FOR THE RIGHTSIDE OF SECTION OUTORDER (right side of ordertemplate)		/////////////////////////////////////////////



	unsigned long int num_con_endhinges_right = findnumconendhingesright();	//finds the number of endhinges that ordertemplate[1] connects to on the right
//	printf("there are %d endhinges that connect to that right side of ordertemplate\n", num_con_endhinges_right);

	if(num_con_endhinges_right==0){
		current_hinge_span[inNum]->num_con_endhinges_right = num_con_endhinges_right;	//number of endhinges that connect to the right of section inNum
	}
	else{
		unsigned long int con_edges_right[num_con_endhinges_right];		//will contain the endhinge "numbers" that connect to the right of section inNum
		unsigned long int index=0;		//index of con_edges_right
		unsigned short int matchresult=0;	//0 if no match, 1 if match
	
		for(i=0; i<=num_con_endhinges_right-1; i++){	//loop through all temp[i=0...num_con_endhinges_right-1]
			//find the corresponding endhinge that contains the same edges as temp[i], then add the endhinge number to con_edges_right[index]
			for(j=0; j<=total_end_hinges-1; j++){	//loop through all endhinges
				matchresult = testmatch(i,j);	//tests if temp[i] matches endhinge[j] (linked lists, doesn't have to be in the same order)
				if(matchresult==1){	//if temp[i] matches endhinge[j]
					con_edges_right[index] = j;
					index++;
					break;
				}
			}
		}
	
		//so con_edges_right contains the indicies of endhinge[] that connect to the right of ordertemplate[1]
//		printf("the endhinges that connect to the right of that ordertemplate are:\n");
		for(i=0; i<=num_con_endhinges_right-1; i++){
//			printendhinge(con_edges_right[i]);
		}
	
		current_hinge_span[inNum]->num_con_endhinges_right = num_con_endhinges_right;	//number of endhinges that connect to the right of section outorder
		current_hinge_span[inNum]->con_endhinges_right = unsgnlong_vecalloc(0, num_con_endhinges_right-1);	//list of endhinge "numbers" that connect to the right of section outorder
		for(i=0; i<=num_con_endhinges_right-1; i++){
			current_hinge_span[inNum]->con_endhinges_right[i] = con_edges_right[i];
		}
	}
	curindex=0;	//reset curindex
	edgesin_temp=0;	//reset edgesin_temp


	//free firsttemp[] and all of it's linked components
	for(i=0; i<=maxendhinges-1; i++){
		currenttemp[i] = firsttemp[i];
		while(currenttemp[i]->nextendhinge != NULL){
			currenttemp[i] = currenttemp[i]->nextendhinge;
			free(firsttemp[i]);
			firsttemp[i] = currenttemp[i];
		}
		free(firsttemp[i]);
	}


//	printf("# of endhinges that connect to the right: %d. ", current_hinge_span[inNum]->num_con_endhinges_right);
	if(current_hinge_span[inNum]->num_con_endhinges_right > 0){
//		printf("Endhinge numbers: ");
		for(i=0; i<=num_con_endhinges_right-1; i++){
//			printf("%d, ", current_hinge_span[inNum]->con_endhinges_right[i]);
//			printendhinge(current_hinge_span[inNum]->con_endhinges_right[i]);

		}
	}


//	printf("\n");
	return;
}


#if defined(CS)
#include "../topology/LFlag_0.c"		/* include if M = 0  */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the left and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#include "../topology/RFlag_0.c"		/* include if M = 0  */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the right and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#include "../utils/noncrossing.c"	/* include if M = 0 */
/* This function takes arguements (int side, int a, int b, int c, int d) */
/* and returns 1 if the walk connecting a and b does not cross */
/* the walk connecting c and d it returns 0 otherwise*/
#include "../sections/cstatenum.c"
#else
#include "../topology/LFlag.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the left and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#include "../topology/RFlag.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the right and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#if vM*vL<7
#include "../sections/Num_section_6V.c"  /*unsigned long int num_section(int side)*/
/* This function takes the arguement (int side) and assigns a unique number */
/* to the section on side = side */
#elif vM*vL<9
#include "../sections/Num_section_8V.c"
#elif vM*vL<11
#include "../sections/Num_section_10V.c"
#endif



#endif


/***************************************************************************/
struct hinge_span *
newhinge(void)		/* returns a pointer to a hinge_span structure */
/* this is a memory allocation function */
{
	struct hinge_span *nextnewhinge;

	nextnewhinge = ((struct hinge_span *) calloc(1, sizeof(struct hinge_span)));
	if (nextnewhinge == NULL) {
		fprintf(stderr, "unable to allocate memory");
		exit(1);
	}
	int i;
	for(i=0; i<=vM*vL*4-1; i++){
		(*nextnewhinge).walk[i]=0;		//initialize to zero
	}
	return nextnewhinge;
}
/***************************************************************************/




void printcurwalk(void){
	printf("curwalk: ");
	int i=0;
	while(curwalk[i]!=0){
		printf("%d, ");
	}
	printf("\n");
}

unsigned long int choose(unsigned long int n, unsigned long int r) {
	if(r>n) {
		return 0;
	}
	return factorial(n)/(factorial(r)*factorial(n-r));
}

unsigned long int factorial(unsigned long int n) {
	if(n==0) {
		return 1;
	}
	if(n==1){
		return 1;
	}
	else {
		return n*factorial(n-1);
	}
}

void generate_endhinges(void){
	unsigned short int i, j, k;		//looping variables
	unsigned short int index=0;		//initially zero

	struct edge edges[hingeedges];		//array of edges.
	//Give each edge in the hinge a number (0 to (hingeedges-1))

	//Label all of the horizontal edges first (left to right, move down, left to right, move down, etc....)
	//By "label", I mean the index number of the array
	for(j=0; j<=M; j++){	//for each row
		for(i=0; i<=L-1; i++){		//for each edge in the row
			edges[index].f.x = i;
			edges[index].f.y = j;
			edges[index].l.x = i+1;
			edges[index].l.y = j;
			index++;
		}
	}

	//Next, label all of the vertical edges (left to right, move down, left to right, move down, etc....)
	//By "label", I mean the index number of the array (continue indexing from after all horizontal edges
	//were generated.
	for(j=0; j<=M-1; j++){		//for each edge in the column
		for(i=0; i<=L; i++){	//for each row
			edges[index].f.x = i;
			edges[index].f.y = j;
			edges[index].l.x = i;
			edges[index].l.y = j+1;
			index++;
		}
	}
//	printf("total number of edges = %d\n", index);

	//print all edges
/*	for(i=0; i<=index-1; i++){
		printf("EDGE %d: (%d,%d) to (%d,%d)\n", i, edges[i].f.x, edges[i].f.y, edges[i].l.x, edges[i].l.y);
	}
*/

	//use this to easily add edges to endhinge....
	//need an algorithm that goes through all combinations of edges in the end hinge.

	//NON-FULL HINGE VERSION:
	
	unsigned short int totedges=1;	//total edges in the endhinge. loop from (1 to hingeedges)
	unsigned long int endnum=0;		//the number/index for the current endhinge
	unsigned long int endhingestoadd;	//given totedges, endhingestoadd is the number of possible endhinges (with totedges edges) that need to be added
	unsigned long int endhingesadded;	//endhinges that have been added (will reset to zero everytime we increase totedges)
	unsigned short int edgestoadd[totedges];	//will hold the edgenumbers that need to be added to a endhinge

	while(totedges <= hingeedges){
//		printf("Will now add all endhinges that contain exactly %d totedge(s)\n", totedges);
		//add all endhinges with "totedges" edges
		//basically there are (hingeedges choose totedges) endhinges that need to be added.
		endhingestoadd = choose(hingeedges, totedges);	//number of possible combinations of endhinges that need to be added
//		printf("there are %d endhinges to be added with %d totedges\n", endhingestoadd, totedges);
		endhingesadded = 0;
		unsigned short int edgestoadd[totedges];	//will contain the edge numbers that will be added

		//get edge numbers that need to be added. store them in edgestoadd[] (find all possible combinations of totedges numbers from (1 to hingeedges)
		iterate(0, edgestoadd, 1, totedges, &endnum, edges, &endhingesadded);

		totedges++;
	}
	total_end_hinges = endnum;
}

void iterate(int k, unsigned short int *edgestoadd, int counter, unsigned short int totedges, unsigned long int* endnum, struct edge* edges, unsigned long int* endhingesadded) {
//	printf("iterate called. k=%d, counter=%d, endnum=%d, endhingesadded=%d\n", k, counter, *endnum, *endhingesadded);
	int i, j, z, a;
	for(i=k; i<=hingeedges-totedges-1+counter; i++){
		edgestoadd[counter-1] = i;	//add edge i
		if(counter==totedges){	//if counter==totedges, then we dont need to recurse any more. (we're at the right depth).
//			printf("counter=%d==totedges=%d\n", counter, totedges);

			//add the totedges to a endhinge structure
			firstendhinge[*endnum] = newendhinge();
			currentendhinge[*endnum] = firstendhinge[*endnum];
			for(j=0; j<=totedges-1; j++){
				//add edgestoadd[j] to endhinge
				currentendhinge[*endnum]->edge.f.x = edges[edgestoadd[j]].f.x;
				currentendhinge[*endnum]->edge.f.y = edges[edgestoadd[j]].f.y;
				currentendhinge[*endnum]->edge.l.x = edges[edgestoadd[j]].l.x;
				currentendhinge[*endnum]->edge.l.y = edges[edgestoadd[j]].l.y;
//				printf("just added edge %d to endhinge #%d\n", edgestoadd[j], *endnum);
				if(j<totedges-1) {	//if there are more edges that need to be added to currentendhinge
//						printf("more edges to add to endhinge #%d\n", *endnum);
					currentendhinge[*endnum]->nextendhinge = newendhinge();
					currentendhinge[*endnum] = currentendhinge[*endnum]->nextendhinge;
				}
				else{
//						printf("NO more edges to add to endhinge #%d\n", *endnum);
				}
			}
//				printf("just FINISHED filling endhinge #%d\n", *endnum);
			(*endnum)++;
			(*endhingesadded)++;
//				printf("endnum now is %d, endhingesadded is %d\n", *endnum, *endhingesadded);

		}
		else{
			iterate(i+1, edgestoadd, counter+1, totedges, endnum, edges, endhingesadded);
		}
	}
}

void printendhinge(int k){
	printf("printing endhinge #%d: \n", k);
	currentendhinge[k] = firstendhinge[k];
	while(currentendhinge[k] != NULL){
		printf("(%d,%d) to (%d,%d),  ", currentendhinge[k]->edge.f.x, currentendhinge[k]->edge.f.y, currentendhinge[k]->edge.l.x, currentendhinge[k]->edge.l.y);
		currentendhinge[k] = currentendhinge[k]->nextendhinge;
	}
	printf("\n");
}
void printtempendhinge(int k){
	printf("printing temp[%d]: \n", k);
	currenttemp[k] = firsttemp[k];
	while(currenttemp[k]->nextendhinge != NULL){
		printf("(%d,%d) to (%d,%d),  ", currenttemp[k]->edge.f.x, currenttemp[k]->edge.f.y, currenttemp[k]->edge.l.x, currenttemp[k]->edge.l.y);
		currenttemp[k] = currenttemp[k]->nextendhinge;
	}
	printf("\n");
}

struct endhinge* newendhinge(void) {	/* this is a memory allocation function */	
					/* returns a pointer to a endhinge structure */
	struct endhinge* nextnewendhinge;

	nextnewendhinge = ((struct endhinge*) calloc(1, sizeof(struct endhinge)));

	if (nextnewendhinge == NULL) {
		fprintf(stderr, "unable to allocate memory for a new endhinge");
		exit(1);
	}
	return nextnewendhinge;
}


unsigned long int findnumconendhingesleft(void){
	// finds and counts how many endhinges connect to left of ordertemplate[0] ("left" side of ordertemplate)
	// returns the count
	unsigned long int i, j;
	unsigned long int connect_count = 0;

	for(i=0; i<=M; i++){		// 0 if empty, 1 if occupied by an ordertemplate edge, 2 if occupied by a hingeedge
		for(j=0; j<=L; j++){
			occupied[i][j] = 0;			//initialize every vertex to be empty
			if(ordertemplate[0][i][j] != 0){	//if there's an edge at (i,j)
				edgesin_temp++;			//count that edge and mark
				occupied[i][j] = 1;			//vertex (i,j) as occupied by an ordertemplate edge
			}
		}
	}

/*	printf("initialized occupied:  ");
	for(i=0; i<=M; i++){		// 0 if empty, 1 if occupied by an ordertemplate edge, 2 if occupied by a hingeedge
		for(j=0; j<=L; j++){
			printf("occupied[%d][%d]=%d, ", i, j, occupied[i][j]);
		}
	}
	printf("\n");
*/
	for (i=0;i<=maxendhinges-1;i++){
		currenttemp[i] = newendhinge();
		firsttemp[i] = currenttemp[i];
	}
	leftside=1;
	connect(1, edgesin_temp, curindex);
	return curindex;
	
}





unsigned long int findnumconendhingesright(void){
	// finds and counts how many endhinges connect to the right of ordertemplate[1] ("right" side of ordertemplate)
	// returns the count
	unsigned long int i, j;
	unsigned long int connect_count = 0;

	for(i=0; i<=M; i++){		// 0 if empty, 1 if occupied by an ordertemplate edge, 2 if occupied by a hingeedge
		for(j=0; j<=L; j++){
			occupied[i][j] = 0;			//initialize every vertex to be empty
			if(ordertemplate[1][i][j] != 0){	//if there's an edge at (i,j)
				edgesin_temp++;			//count that edge and mark
				occupied[i][j] = 1;			//vertex (i,j) as occupied by an ordertemplate edge
			}
		}
	}

/*	printf("initialized occupied:  ");
	for(i=0; i<=M; i++){		// 0 if empty, 1 if occupied by an ordertemplate edge, 2 if occupied by a hingeedge
		for(j=0; j<=L; j++){
			printf("occupied[%d][%d]=%d, ", i, j, occupied[i][j]);
		}
	}
	printf("\n");
*/
	for (i=0;i<=maxendhinges-1;i++){
		currenttemp[i] = newendhinge();
		firsttemp[i] = currenttemp[i];
	}
	leftside=0;

	connect(1, 2, curindex);

	return curindex;
	
}




void connect(unsigned long int small, unsigned long int big, unsigned long int cindex){

	struct endhinge* firstmemory = newendhinge();	//use these to remember what the endhinge was when the function was originally called
	struct endhinge* currentmemory = firstmemory;	//this way, we can go back and use what the endhinge was when the function was called
	currenttemp[cindex] = firsttemp[cindex];

	while(currenttemp[cindex]->nextendhinge != NULL){
		currentmemory->edge.f.x = currenttemp[cindex]->edge.f.x;
		currentmemory->edge.f.y = currenttemp[cindex]->edge.f.y;
		currentmemory->edge.l.x = currenttemp[cindex]->edge.l.x;
		currentmemory->edge.l.y = currenttemp[cindex]->edge.l.y;

		currenttemp[cindex] = currenttemp[cindex]->nextendhinge;
		currentmemory->nextendhinge = newendhinge();
		currentmemory = currentmemory->nextendhinge;
	}

	//adds all endhinges that connect the point where small is to the point where big is.
//	printf("CONNECT called. connecting edge %d to %d. cindex=%d\n", small, big, cindex);

	unsigned long int i, j;
	struct point start;	//starting point (where small is)
	struct point end;	//end point (where big is)

	//assign the proper x and y values for the start and end points
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			if(ordertemplate[leftside ^ 1][i][j] == small){
				start.x = j;
				start.y = i;
			}
			else if(ordertemplate[leftside ^ 1][i][j] == big){
				end.x = j;
				end.y = i;
			}
		}
	}

	// explore first going up
	if(start.y > 0){
//		printf("there is room to go up\n");
		if(occupied[start.y-1][start.x] == 0){
//			printf("can go up. it is empty\n");
			goup(start.x, start.y, end.x, end.y, small, big, curindex);
//			printf("finished going up\n");
			occupied[start.y-1][start.x] = 0;	//change occupied back to empty at that point. will now explore other options
		}
		else if(start.y-1==end.y && start.x==end.x){
//			printf("can go up. it will connect you to the end\n");
			goup(start.x, start.y, end.x, end.y, small, big, curindex);
//			printf("finished going up\n");
		}
		else{
//			printf("but can't go up because it's already occupied\n");
		}
	}
	else{
//		printf("no room to go up\n");
	}

	//copy memory into temp[curindex]
	currenttemp[curindex] = firsttemp[curindex];
	currentmemory = firstmemory;

	while(currentmemory->nextendhinge != NULL){
		currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
		currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
		currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
		currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;

		currentmemory = currentmemory->nextendhinge;
		currenttemp[curindex]->nextendhinge = newendhinge();
		currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
	}

	// explore first going down
	if(start.y < M){
//		printf("there is room to go down\n");
		if(occupied[start.y+1][start.x] == 0){
//			printf("can go down. it is empty\n");
			godown(start.x, start.y, end.x, end.y, small, big, curindex);
//			printf("finished going down\n");
			occupied[start.y+1][start.x] = 0;	//change occupied back to empty at that point. will now explore other options
		}
		else if(start.y+1==end.y && start.x==end.x){
//			printf("can go down. it will connect you to the end\n");
			godown(start.x, start.y, end.x, end.y, small, big, curindex);
//			printf("finished going down\n");
		}
		else{
//			printf("but can't go down because it's occupied by a hingeedge\n");
		}
	}
	else{
//		printf("no room to go down\n");
	}

	//copy memory into temp[curindex]
	currenttemp[curindex] = firsttemp[curindex];
	currentmemory = firstmemory;

	while(currentmemory->nextendhinge != NULL){
		currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
		currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
		currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
		currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;

		currentmemory = currentmemory->nextendhinge;
		currenttemp[curindex]->nextendhinge = newendhinge();
		currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
	}


	// explore first going left
	if(start.x>0){
//		printf("there is room to go left\n");
		if(occupied[start.y][start.x-1] == 0){
//			printf("can go left. it is empty\n");
			goleft(start.x, start.y, end.x, end.y, small, big, curindex);
//			printf("finished going left\n");
			occupied[start.y][start.x-1] = 0;	//change occupied back to empty at that point. will now explore other options
		}
		else if(start.y==end.y && start.x-1==end.x){
//			printf("can go left. it will connect you to the end\n");
			goleft(start.x, start.y, end.x, end.y, small, big, curindex);
//			printf("finished going left\n");
		}
		else{
//			printf("but can't go left because it's occupied by a hingeedge\n");
		}
	}
	else{
//		printf("no room to go left\n");
	}

	//copy memory into temp[curindex]
	currenttemp[curindex] = firsttemp[curindex];
	currentmemory = firstmemory;

	while(currentmemory->nextendhinge != NULL){
		currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
		currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
		currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
		currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;

		currentmemory = currentmemory->nextendhinge;
		currenttemp[curindex]->nextendhinge = newendhinge();
		currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
	}


	// explore first going right
	if(start.x<L){
//		printf("there is room to go right\n");
		if(occupied[start.y][start.x+1] == 0){
//			printf("can go right. it is empty\n");
			goright(start.x, start.y, end.x, end.y, small, big, curindex);
//			printf("finished going right\n");
			occupied[start.y][start.x+1] = 0;	//change occupied back to empty at that point. will now explore other options
		}
		else if(start.y==end.y && start.x+1==end.x){
//			printf("can go right. it will connect you to the end\n");
			goright(start.x, start.y, end.x, end.y, small, big, curindex);
//			printf("finished going right\n");
		}
		else{
//			printf("but can't go right because it's occupied by a hingeedge\n");
		}
	}
	else{
//		printf("no room to go right\n");
	}

	//free firstmemory and all of it's linked components
	currentmemory=firstmemory;
	while(currentmemory->nextendhinge != NULL){
		currentmemory=currentmemory->nextendhinge;
		free(firstmemory);
		firstmemory=currentmemory;
	}
	free(firstmemory);

}



void goup(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int uindex){

//	printf("going up called. prevx=%d, prevy=%d, endx=%d, endy=%d, firstconnum=%d, secondconnum=%d\n", prevx, prevy, endx, endy, firstconnum, secondconnum);
//	printf("GOING UP. ");
	//add edge to temp[curindex]
	currenttemp[curindex]->edge.f.x = prevx;
	currenttemp[curindex]->edge.f.y = prevy-1;
	currenttemp[curindex]->edge.l.x = prevx;
	currenttemp[curindex]->edge.l.y = prevy;
	currenttemp[curindex]->nextendhinge = newendhinge();
	currenttemp[curindex] = currenttemp[curindex]->nextendhinge;

//	printf("temp[%d] is now: ", curindex);
//	printtempendhinge(curindex);

	struct endhinge* firstmemory = newendhinge();	//use these to remember what the endhinge was when the function was originally called
	struct endhinge* currentmemory = firstmemory;	//this way, we can go back and use what the endhinge was when the function was called
	currenttemp[uindex] = firsttemp[uindex];

	while(currenttemp[uindex]->nextendhinge != NULL){
		currentmemory->edge.f.x = currenttemp[uindex]->edge.f.x;
		currentmemory->edge.f.y = currenttemp[uindex]->edge.f.y;
		currentmemory->edge.l.x = currenttemp[uindex]->edge.l.x;
		currentmemory->edge.l.y = currenttemp[uindex]->edge.l.y;

		currenttemp[uindex] = currenttemp[uindex]->nextendhinge;
		currentmemory->nextendhinge = newendhinge();
		currentmemory = currentmemory->nextendhinge;
	}

//	printf("edge added to temp[%d]\n", curindex);
	if(occupied[prevy-1][prevx]!=1){	//if vertex prevy-1, prevx was empty
		occupied[prevy-1][prevx] = 2;	//vertex prevy-1, prevx is now occupied by a hingeedge, unless it was already occupied by an ordertemplate edge
	}
//	printf("occupied[%d][%d]=%d now\n", prevy-1, prevx, occupied[prevy-1][prevx]);
	//check if we're at the vertex where end point is.
	if(prevx==endx && (prevy-1)==endy){	//if firstconnum connects to secondconnum
//		printf("we've connected %d to %d. check if there are more points that we need to connect\n", firstconnum, secondconnum);

		//check if there is another set of 2 points that need to be connected
		if((leftside==1) && ((secondconnum+1==edgesin_temp) || (firstconnum==1&&secondconnum==2))   ||    (leftside==0 && secondconnum==edgesin_temp)){	//no more points to connect

//			printf("no more points need to be connected. increase curindex from %d to %d. temp[%d] is now a valid endhinge\n", curindex, curindex+1, curindex);
//			printtempendhinge(curindex);
			curindex++;			//finished temp[curindex], look for other endhinges that connect ordertemplate.
		}

		else{	//more points to connect
//			printf("there are still more points to connect\n");
			if(leftside==1&&firstconnum==1){				//if just connected the "first and last" points
				connect(firstconnum+1, firstconnum+2, curindex);	//connect sequential points
			}
			else{
				connect(firstconnum+2, secondconnum+2, curindex);
			}
		}
	}
	else{	//we're not at the endpoint	//explore going UP again, LEFT, OR RIGHT (can't go down)
//		printf("we're not connected yet. try going up again, left, and right\n");
		// explore going up
		if(prevy-1 > 0){
//			printf("there is room to go up\n");
			if(occupied[prevy-1-1][prevx] == 0){
//				printf("can go up. it is empty\n");
				goup(prevx, prevy-1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going up\n");
				occupied[prevy-1-1][prevx] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy-1-1, prevx, occupied[prevy-1-1][prevx]);
			}
			else if(prevy-1-1==endy && prevx==endx){
//				printf("can go up. it will connect you to the end\n");
				goup(prevx, prevy-1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going up\n");
			}
			else{
//				printf("but can't go up because it's already occupied\n");
			}
		}
		else{
//			printf("no room to go up\n");
		}
	

		//copy memory into temp[curindex]
		currenttemp[curindex] = firsttemp[curindex];
		currentmemory = firstmemory;
	
		while(currentmemory->nextendhinge != NULL){
			currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
			currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
			currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
			currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;
	
			currentmemory = currentmemory->nextendhinge;
			currenttemp[curindex]->nextendhinge = newendhinge();
			currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
		}


		// explore going left
		if(prevx>0){
//			printf("there is room to go left\n");
			if(occupied[prevy-1][prevx-1] == 0){
//				printf("can go left. it is empty\n");
				goleft(prevx, prevy-1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going left\n");
				occupied[prevy-1][prevx-1] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy-1, prevx-1, occupied[prevy-1][prevx-1]);
			}
			else if(prevy-1==endy && prevx-1==endx){
//				printf("can go left. it will connect you to the end\n");
				goleft(prevx, prevy-1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going left\n");
			}
			else{
//				printf("but can't go left because it's occupied by a hingeedge\n");
			}
		}
		else{
//			printf("no room to go left\n");
		}
	
		//copy memory into temp[curindex]
		currenttemp[curindex] = firsttemp[curindex];
		currentmemory = firstmemory;
	
		while(currentmemory->nextendhinge != NULL){
			currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
			currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
			currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
			currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;
	
			currentmemory = currentmemory->nextendhinge;
			currenttemp[curindex]->nextendhinge = newendhinge();
			currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
		}


		// explore going right
		if(prevx<L){
//			printf("there is room to go right\n");
			if(occupied[prevy-1][prevx+1] == 0){
//				printf("can go right. it is empty\n");
				goright(prevx, prevy-1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going right\n");
				occupied[prevy-1][prevx+1] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy-1, prevx+1, occupied[prevy-1][prevx+1]);
			}
			else if(prevy-1==endy && prevx+1==endx){
//				printf("can go right. it will connect you to the end\n");
				goright(prevx, prevy-1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going right\n");
			}
			else{
//				printf("but can't go right because it's occupied by a hingeedge\n");
			}
		}
		else{
//			printf("no room to go right\n");
		}
	}
}

void godown(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int dindex){

//	printf("going down called. prevx=%d, prevy=%d, endx=%d, endy=%d, firstconnum=%d, secondconnum=%d\n", prevx, prevy, endx, endy, firstconnum, secondconnum);
//	printf("GOING DOWN. ");
	//add edge to temp[curindex]
	currenttemp[curindex]->edge.f.x = prevx;
	currenttemp[curindex]->edge.f.y = prevy;
	currenttemp[curindex]->edge.l.x = prevx;
	currenttemp[curindex]->edge.l.y = prevy+1;
	currenttemp[curindex]->nextendhinge = newendhinge();
	currenttemp[curindex] = currenttemp[curindex]->nextendhinge;

//	printf("temp[%d] is now: ", curindex);
//	printtempendhinge(curindex);

	struct endhinge* firstmemory = newendhinge();	//use these to remember what the endhinge was when the function was originally called
	struct endhinge* currentmemory = firstmemory;	//this way, we can go back and use what the endhinge was when the function was called
	currenttemp[dindex] = firsttemp[dindex];

	while(currenttemp[dindex]->nextendhinge != NULL){
		currentmemory->edge.f.x = currenttemp[dindex]->edge.f.x;
		currentmemory->edge.f.y = currenttemp[dindex]->edge.f.y;
		currentmemory->edge.l.x = currenttemp[dindex]->edge.l.x;
		currentmemory->edge.l.y = currenttemp[dindex]->edge.l.y;

		currenttemp[dindex] = currenttemp[dindex]->nextendhinge;
		currentmemory->nextendhinge = newendhinge();
		currentmemory = currentmemory->nextendhinge;
	}

//	printf("edge added to temp[%d]\n", curindex);
	if(occupied[prevy+1][prevx]!=1){	//if vertex prevy+1, prevx was empty
		occupied[prevy+1][prevx] = 2;	//vertex prevy+1, prevx is now occupied by a hingeedge, unless it was already occupied by an ordertemplate edge
	}
//	printf("occupied[%d][%d]=%d now\n", prevy+1, prevx, occupied[prevy+1][prevx]);
	//check if we're at the vertex where end point is
	if(prevx==endx && (prevy+1)==endy){	//if firstconnum connects to secondconnum
//		printf("we've connected %d to %d. check if there are more points that we need to connect\n", firstconnum, secondconnum);
		//check if there is another set of 2 points that need to be connected
		if((leftside==1) && ((secondconnum+1==edgesin_temp) || (firstconnum==1&&secondconnum==2))   ||    (leftside==0 && secondconnum==edgesin_temp)){	//no more points to connect

//			printf("no more points need to be connected. increase curindex from %d to %d. temp[%d] is now a valid endhinge\n", curindex, curindex+1, curindex);
//			printtempendhinge(curindex);
			curindex++;			//finished temp[curindex], look for other endhinges that connect ordertemplate.

		}
		else{	//more points to connect
//			printf("there are still more points to connect\n");
			if(leftside==1&&firstconnum==1){				//if just connected the "first and last" points
				connect(firstconnum+1, firstconnum+2, curindex);	//connect sequential points
			}
			else{
				connect(firstconnum+2, secondconnum+2, curindex);
			}
		}
	}
	else{	//we're not at the endpoint	//explore going DOWN again, LEFT, OR RIGHT (can't go down)
//		printf("we're not connected yet. try going down again, left, and right\n");

		// explore going down
		if(prevy+1 < M){
//			printf("there is room to go down\n");
			if(occupied[prevy+1+1][prevx] == 0){
//				printf("can go down. it is empty\n");
				godown(prevx, prevy+1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going down\n");
				occupied[prevy+1+1][prevx] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy+1+1, prevx, occupied[prevy+1+1][prevx]);
			}
			else if(prevy+1+1==endy && prevx==endx){
//				printf("can go down. it will connect you to the end\n");
				godown(prevx, prevy+1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going down\n");
			}
			else{
//				printf("but can't go down because it's already occupied\n");
			}
		}
		else{
//			printf("no room to go down\n");
		}
	

		//copy memory into temp[curindex]
		currenttemp[curindex] = firsttemp[curindex];
		currentmemory = firstmemory;
	
		while(currentmemory->nextendhinge != NULL){
			currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
			currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
			currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
			currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;
	
			currentmemory = currentmemory->nextendhinge;
			currenttemp[curindex]->nextendhinge = newendhinge();
			currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
		}


		// explore going left
		if(prevx>0){
//			printf("there is room to go left\n");
			if(occupied[prevy+1][prevx-1] == 0){
//				printf("can go left. it is empty\n");
				goleft(prevx, prevy+1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going left\n");
				occupied[prevy+1][prevx-1] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy+1, prevx-1, occupied[prevy+1][prevx-1]);
			}
			else if(prevy+1==endy && prevx-1==endx){
//				printf("can go left. it will connect you to the end\n");
				goleft(prevx, prevy+1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going left\n");
			}
			else{
//				printf("but can't go left because it's occupied by a hingeedge\n");
			}
		}
		else{
//			printf("no room to go left\n");
		}
	
		//copy memory into temp[curindex]
		currenttemp[curindex] = firsttemp[curindex];
		currentmemory = firstmemory;
	
		while(currentmemory->nextendhinge != NULL){
			currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
			currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
			currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
			currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;
	
			currentmemory = currentmemory->nextendhinge;
			currenttemp[curindex]->nextendhinge = newendhinge();
			currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
		}


		// explore going right
		if(prevx<L){
//			printf("there is room to go right\n");
			if(occupied[prevy+1][prevx+1] == 0){
//				printf("can go right. it is empty\n");
				goright(prevx, prevy+1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going right\n");
				occupied[prevy+1][prevx+1] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy+1, prevx+1, occupied[prevy+1][prevx+1]);
			}
			else if(prevy+1==endy && prevx+1==endx){
//				printf("can go right. it will connect you to the end\n");
				goright(prevx, prevy+1, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going right\n");
			}
			else{
//				printf("but can't go right because it's occupied by a hingeedge\n");
			}
		}
		else{
//			printf("no room to go right\n");
		}
	}
}

void goleft(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int lindex){


//	printf("going left called. prevx=%d, prevy=%d, endx=%d, endy=%d, firstconnum=%d, secondconnum=%d\n", prevx, prevy, endx, endy, firstconnum, secondconnum);
//	printf("GOING LEFT. ");
	//add edge to temp[curindex]
	currenttemp[curindex]->edge.f.x = prevx-1;
	currenttemp[curindex]->edge.f.y = prevy;
	currenttemp[curindex]->edge.l.x = prevx;
	currenttemp[curindex]->edge.l.y = prevy;
	currenttemp[curindex]->nextendhinge = newendhinge();
	currenttemp[curindex] = currenttemp[curindex]->nextendhinge;

//	printf("temp[%d] is now: ", curindex);
//	printtempendhinge(curindex);

	struct endhinge* firstmemory = newendhinge();	//use these to remember what the endhinge was when the function was originally called
	struct endhinge* currentmemory = firstmemory;	//this way, we can go back and use what the endhinge was when the function was called
	currenttemp[lindex] = firsttemp[lindex];

	while(currenttemp[lindex]->nextendhinge != NULL){
		currentmemory->edge.f.x = currenttemp[lindex]->edge.f.x;
		currentmemory->edge.f.y = currenttemp[lindex]->edge.f.y;
		currentmemory->edge.l.x = currenttemp[lindex]->edge.l.x;
		currentmemory->edge.l.y = currenttemp[lindex]->edge.l.y;

		currenttemp[lindex] = currenttemp[lindex]->nextendhinge;
		currentmemory->nextendhinge = newendhinge();
		currentmemory = currentmemory->nextendhinge;
	}

//	printf("edge added to temp[%d]\n", curindex);
	if(occupied[prevy][prevx-1] != 1){	//if vertex prevy, prevx-1 was empty
		occupied[prevy][prevx-1] = 2;	//vertex prevy, prevx-1 is now occupied by a hingeedge, unless it was already occupied by an ordertemplate edge
	}
//	printf("occupied[%d][%d]=%d now\n", prevy, prevx-1, occupied[prevy][prevx-1]);
	//check if we're at the vertex where end point is
	if(prevx-1==endx && prevy==endy){	//if firstconnum connects to secondconnum
//		printf("we've connected %d to %d. check if there are more points that we need to connect\n", firstconnum, secondconnum);
		//check if there is another set of 2 points that need to be connected
		if((leftside==1) && ((secondconnum+1==edgesin_temp) || (firstconnum==1&&secondconnum==2))   ||    (leftside==0 && secondconnum==edgesin_temp)){ //no more points to connect

//			printf("no more points need to be connected. increase curindex from %d to %d. temp[%d] is now a valid endhinge\n", curindex, curindex+1, curindex);
//			printtempendhinge(curindex);
			curindex++;			//finished temp[curindex], look for other endhinges that connect ordertemplate.

		}
		else{	//more points to connect
//			printf("there are still more points to connect\n");
			if(leftside==1&&firstconnum==1){				//if just connected the "first and last" points
				connect(firstconnum+1, firstconnum+2, curindex);	//connect sequential points
			}
			else{
				connect(firstconnum+2, secondconnum+2, curindex);
			}
		}
	}
	else{	//we're not at the endpoint	//explore going LEFT again, UP, OR DOWN (can't go RIGHT)
//		printf("we're not connected yet. try going left again, up, or down.\n");

		// explore going left
		if(prevx-1>0){
//			printf("there is room to go left\n");
			if(occupied[prevy][prevx-1-1] == 0){
//				printf("can go left. it is empty\n");
				goleft(prevx-1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going left\n");
				occupied[prevy][prevx-1-1] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy, prevx-1-1, occupied[prevy][prevx-1-1]);
			}
			else if(prevy==endy && prevx-1-1==endx){
//				printf("can go left. it will connect you to the end\n");
				goleft(prevx-1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going left\n");
			}
			else{
//				printf("but can't go left because it's occupied by a hingeedge\n");
			}
		}
		else{
//			printf("no room to go left\n");
		}

	

		//copy memory into temp[curindex]
		currenttemp[curindex] = firsttemp[curindex];
		currentmemory = firstmemory;
	
		while(currentmemory->nextendhinge != NULL){
			currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
			currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
			currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
			currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;
	
			currentmemory = currentmemory->nextendhinge;
			currenttemp[curindex]->nextendhinge = newendhinge();
			currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
		}


		// explore going up
		if(prevy > 0){
//			printf("there is room to go up\n");
			if(occupied[prevy-1][prevx-1] == 0){
//				printf("can go up. it is empty\n");
				goup(prevx-1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going up\n");
				occupied[prevy-1][prevx-1] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy-1, prevx-1, occupied[prevy-1][prevx-1]);
			}
			else if(prevy-1==endy && prevx-1==endx){
//				printf("can go up. it will connect you to the end\n");
				goup(prevx-1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going up\n");
			}
			else{
//				printf("but can't go up because it's already occupied\n");
			}
		}
		else{
//			printf("no room to go up\n");
		}

		//copy memory into temp[curindex]
		currenttemp[curindex] = firsttemp[curindex];
		currentmemory = firstmemory;
	
		while(currentmemory->nextendhinge != NULL){
			currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
			currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
			currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
			currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;
	
			currentmemory = currentmemory->nextendhinge;
			currenttemp[curindex]->nextendhinge = newendhinge();
			currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
		}


		// explore going down
		if(prevy < M){
//			printf("there is room to go down\n");
			if(occupied[prevy+1][prevx-1] == 0){
//				printf("can go down. it is empty\n");
				godown(prevx-1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going down\n");
				occupied[prevy+1][prevx-1] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy+1, prevx-1, occupied[prevy+1][prevx-1]);
			}
			else if(prevy+1==endy && prevx-1==endx){
//				printf("can go down. it will connect you to the end\n");
				godown(prevx-1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going down\n");
			}
			else{
//				printf("but can't go down because it's already occupied\n");
			}
		}
		else{
//			printf("no room to go down\n");
		}
	}
}

void goright(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int rindex){


//	printf("going right called. prevx=%d, prevy=%d, endx=%d, endy=%d, firstconnum=%d, secondconnum=%d\n", prevx, prevy, endx, endy, firstconnum, secondconnum);
//	printf("GOING RIGHT. ");
	//add edge to temp[curindex]
	currenttemp[curindex]->edge.f.x = prevx;
	currenttemp[curindex]->edge.f.y = prevy;
	currenttemp[curindex]->edge.l.x = prevx+1;
	currenttemp[curindex]->edge.l.y = prevy;
	currenttemp[curindex]->nextendhinge = newendhinge();
	currenttemp[curindex] = currenttemp[curindex]->nextendhinge;

//	printf("temp[%d] is now: ", curindex);
//	printtempendhinge(curindex);

	struct endhinge* firstmemory = newendhinge();	//use these to remember what the endhinge was when the function was originally called
	struct endhinge* currentmemory = firstmemory;	//this way, we can go back and use what the endhinge was when the function was called
	currenttemp[rindex] = firsttemp[rindex];

	while(currenttemp[rindex]->nextendhinge != NULL){
		currentmemory->edge.f.x = currenttemp[rindex]->edge.f.x;
		currentmemory->edge.f.y = currenttemp[rindex]->edge.f.y;
		currentmemory->edge.l.x = currenttemp[rindex]->edge.l.x;
		currentmemory->edge.l.y = currenttemp[rindex]->edge.l.y;

		currenttemp[rindex] = currenttemp[rindex]->nextendhinge;
		currentmemory->nextendhinge = newendhinge();
		currentmemory = currentmemory->nextendhinge;
	}

//	printf("edge added to temp[%d]\n", curindex);
	if(occupied[prevy][prevx+1] != 1){	//if vertex prevy, prevx+1 was empty
		occupied[prevy][prevx+1] = 2;	//vertex prevy, prevx+1 is now occupied by a hingeedge, unless it was already occupied by an ordertemplate edge
	}
//	printf("occupied[%d][%d]=%d now\n", prevy, prevx+1, occupied[prevy][prevx+1]);
	//check if we're at the vertex where end point is
	if(prevx+1==endx && prevy==endy){	//if firstconnum connects to secondconnum
//		printf("we've connected %d to %d. check if there are more points that we need to connect\n", firstconnum, secondconnum);
		//check if there is another set of 2 points that need to be connected
		if((leftside==1) && ((secondconnum+1==edgesin_temp) || (firstconnum==1&&secondconnum==2))   ||    (leftside==0 && secondconnum==edgesin_temp)){	//no more points to connect


//			printf("no more points need to be connected. increase curindex from %d to %d. temp[%d] is now a valid endhinge\n", curindex, curindex+1, curindex);
//			printtempendhinge(curindex);
			curindex++;			//finished temp[curindex], look for other endhinges that connect ordertemplate.

		}
		else{	//more points to connect
//			printf("there are still more points to connect\n");
			if(leftside==1&&firstconnum==1){				//if just connected the "first and last" points
				connect(firstconnum+1, firstconnum+2, curindex);	//connect sequential points
			}
			else{
				connect(firstconnum+2, secondconnum+2, curindex);
			}
		}
	}
	else{	//we're not at the endpoint	//explore going RIGHT again, UP, OR DOWN (can't go LEFT)
//		printf("we're not connected yet. try going left again, up, or down.\n");

		// explore going right
		if(prevx+1<L){
//			printf("there is room to go right\n");
			if(occupied[prevy][prevx+1+1] == 0){
//				printf("can go right. it is empty\n");
				goright(prevx+1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going right\n");
				occupied[prevy][prevx+1+1] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy, prevx+1+1, occupied[prevy][prevx+1+1]);
			}
			else if(prevy==endy && prevx+1+1==endx){
//				printf("can go right. it will connect you to the end\n");
				goright(prevx+1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going right\n");
			}
			else{
//				printf("but can't go right because it's occupied by a hingeedge\n");
			}
		}
		else{
//			printf("no room to go right\n");
		}
	

		//copy memory into temp[curindex]
		currenttemp[curindex] = firsttemp[curindex];
		currentmemory = firstmemory;
	
		while(currentmemory->nextendhinge != NULL){
			currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
			currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
			currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
			currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;
	
			currentmemory = currentmemory->nextendhinge;
			currenttemp[curindex]->nextendhinge = newendhinge();
			currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
		}


		// explore going up
		if(prevy > 0){
//			printf("there is room to go up\n");
			if(occupied[prevy-1][prevx+1] == 0){
//				printf("can go up. it is empty\n");
				goup(prevx+1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going up\n");
				occupied[prevy-1][prevx+1] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy-1, prevx+1, occupied[prevy-1][prevx+1]);
			}
			else if(prevy-1==endy && prevx+1==endx){
//				printf("can go up. it will connect you to the end\n");
				goup(prevx+1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going up\n");
			}
			else{
//				printf("but can't go up because it's already occupied\n");
			}
		}
		else{
//			printf("no room to go up\n");
		}

		//copy memory into temp[curindex]
		currenttemp[curindex] = firsttemp[curindex];
		currentmemory = firstmemory;
	
		while(currentmemory->nextendhinge != NULL){
			currenttemp[curindex]->edge.f.x = currentmemory->edge.f.x;
			currenttemp[curindex]->edge.f.y = currentmemory->edge.f.y;
			currenttemp[curindex]->edge.l.x = currentmemory->edge.l.x;
			currenttemp[curindex]->edge.l.y = currentmemory->edge.l.y;
	
			currentmemory = currentmemory->nextendhinge;
			currenttemp[curindex]->nextendhinge = newendhinge();
			currenttemp[curindex] = currenttemp[curindex]->nextendhinge;
		}


		// explore going down
		if(prevy < M){
//			printf("there is room to go down\n");
			if(occupied[prevy+1][prevx+1] == 0){
//				printf("can go down. it is empty\n");
				godown(prevx+1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going down\n");
				occupied[prevy+1][prevx+1] = 0;	//change occupied back to empty at that point. will now explore other options
//				printf("occupied[%d][%d]=%d now\n", prevy+1, prevx+1, occupied[prevy+1][prevx+1]);
			}
			else if(prevy+1==endy && prevx+1==endx){
//				printf("can go down. it will connect you to the end\n");
				godown(prevx+1, prevy, endx, endy, firstconnum, secondconnum, curindex);
//				printf("finished going down\n");
			}
			else{
//				printf("but can't go down because it's already occupied\n");
			}
		}
		else{
//			printf("no room to go down\n");
		}
	}
}
unsigned short int testmatch(int i, int j){
	//returns 1 if temp[i] has the same edges as endhinge[j] (don't have to be in the same order)
	//returns 0 if they don't have the same edges

//	printf("testing if temp[%d] matches endhinge[%d]:\n", i, j);
//	printtempendhinge(i);
//	printendhinge(j);
	currenttemp[i] = firsttemp[i];
	currentendhinge[j] = firstendhinge[j];

	//first, make sure they have the same number of edges
	unsigned long int tempedges=0;
	unsigned long int endhingeedges=0;
	while(currenttemp[i]->nextendhinge != NULL){	//loop through all edges
		tempedges++;
		currenttemp[i] = currenttemp[i]->nextendhinge;
	}
	while(currentendhinge[j] != NULL){	//loop through all edges
		endhingeedges++;
		currentendhinge[j] = currentendhinge[j]->nextendhinge;
	}
//	printf("temp has %d edges, endhinge has %d edges\n", tempedges, endhingeedges);
	if(tempedges==endhingeedges){	//possible match. now test if they have the same edges
		unsigned short int foundedge=0;

		currenttemp[i] = firsttemp[i];
		currentendhinge[j] = firstendhinge[j];
		while(currenttemp[i]->nextendhinge != NULL){	//go through all edges in temp[i]
			foundedge=0;//haven't found an edge that matches currenttemp[i]'s edge
			currentendhinge[j]=firstendhinge[j];
			while(currentendhinge[j] != NULL){	//go through all edges in endhinge[j]
				if((currenttemp[i]->edge.f.x==currentendhinge[j]->edge.f.x) && (currenttemp[i]->edge.f.y==currentendhinge[j]->edge.f.y) && (currenttemp[i]->edge.l.x==currentendhinge[j]->edge.l.x) && \
				(currenttemp[i]->edge.l.y==currentendhinge[j]->edge.l.y)){	//if currenttemp[i]'s edge matches currentendhinge's edge
					foundedge=1;
					break;
				}
				currentendhinge[j] = currentendhinge[j]->nextendhinge;
			}
			if(foundedge==0){	//if no edge in currentendhinge[j] matches currenttemp[i]
//				printf("temp[%d] DOESN'T MATCH endhinge[%d]!!\n", i, j);
				return 0;	//temp[i] doesn't match endhinge[j]
			}
			currenttemp[i] = currenttemp[i]->nextendhinge;
		}
		//if you made it out of the loop, then all edges in temp[i] matched some edge in endhinge[j]
//		printf("temp[%d] MATCHES endhinge[%d]!!\n", i, j);
		return 1;
		
	}
	else{
//		printf("temp[%d] doesn't have the same number of edges as endhinge[%d]. they don't match\n", i, j);
		return 0;
	}
	
}




int dupcheck(int x, int y, int dup, int firstsection){

	//will recursively check if (*current_hinge_span[firstsection]).walk[i] will be a duplicate
	//check previous "starting" point prior to (x,y). (move there, then check if first step is a 5 or last step is a 6).
	if(dup==1){
		return dup;
	}

	//recursively check through all previous starting points before (x,y)   ordering of starting points is (0,0) (1,0) .... (L,0) (0,1) (1,1) ... (L,1) ... (0,2)... (L,2) ... (L,M)
	//move to new point, then check if first step is 5, or last step is 6. if it is, then we have a duplicate. (it must've already been generated)
	if(x>0){
		x--;
	}
	else if(y>0){
		x=L;
		y--;
	}
	else{//x==0 and y==0;
		return dup;
	}
	//now found proper x,y we need to test.


	int curx=(*current_hinge_span[firstsection]).startpoint[0];
	int cury=(*current_hinge_span[firstsection]).startpoint[1];
	int curz=0;

//	printf("checking for duplicate at point (%d, %d), startpoint is (%d,%d)\n", x, y, curx, cury);

	//put walk[] into dupcheckwalk[]
	
	int i=0;
	int dupcheckwalk[pwalksize];	//this will store the walk that is being checked if it's a duplicate

	for(i=0; i<=pwalksize-1; i++){	//initialize dupcheckwalk to be all zeros
		dupcheckwalk[i]=0;
	}

	//copy walk[] into dupcheck[] walk
	i=0;
//	printf("dupcheck walk is: ");
	while((*current_hinge_span[firstsection]).walk[i] != 0){
		dupcheckwalk[i] = (*current_hinge_span[firstsection]).walk[i];
//		printf("%d, ", dupcheckwalk[i]);
		i++;
	}
//	printf("\n");

	//move along dupcheckwalk[] until we get to (x,y,0)
//	printf("going to start moving through dupcheckwalk[]\n");
	i=0;
	while(((curx!=x || cury!=y) || curz!=0) || (dupcheckwalk[i-1]<0 && dupcheckwalk[i]<0)){	//while not at (x,y,0),
//		printf("curx=%d, cury=%d, curz=%d. Not at (%d, %d, 0) yet!!!\n", curx, cury, curz, x, y);
		if(!((curx!=x || cury!=y) || curz!=0) && (dupcheckwalk[i-1]<0 && dupcheckwalk[i]<0)){
			curious++;
//			printf("At check point (%d,%d,0), but pretty sure it's a fake. curious=%d.\n", x, y, curious);
//			if(curious==4){
//				exit(1);
//			}
		}
		if(dupcheckwalk[i]==1 || dupcheckwalk[i]==-1){
			curx++;
		}
		else if(dupcheckwalk[i]==2 || dupcheckwalk[i]==-2){
			curx--;
		}
		else if(dupcheckwalk[i]==3 || dupcheckwalk[i]==-3){
			cury++;
		}
		else if(dupcheckwalk[i]==4 || dupcheckwalk[i]==-4){
			cury--;
		}
		else if(dupcheckwalk[i]==5){
			curz++;
		}
		else if(dupcheckwalk[i]==6){
			curz--;
		}
		else if(dupcheckwalk[i]==0){	//reached end of the walk[] (never gets to (x, y, 0))
			i=i+1;	//just to be safe
			break;
		}
		i++;
//		printf("curx=%d, cury=%d, curz=%d\n", curx, cury, curz);
	}


	//now at (x,y,0) (before executing dupcheckwalk[i])
	if(dupcheckwalk[i-1]==6){	//just moved negz, so we have a duplicate
//		printf("duplicate found with last move was negz (x=%d, y=%d)\n", x, y);
		return 1;	//we have a duplicate
	}
	else if(dupcheckwalk[i]==5){	//next move was to move posz (add the +1 because of the neg.)
//		printf("duplicate found with next move of posz(x=%d, y=%d)\n", x, y);
		return 1;	//we have a duplicate
	}
	else{
		dup=dupcheck(x, y, dup, firstsection);
	}


	return dup;
	
}


int test2secs(int sectionnum) {
/* This function takes in a section number and returns 1 if it consists of only 2 edges... 0 otherwise
*/
	int i, j;
	int i1, j1, i2, j2;

	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			ordertemplate[0][i][j]=0;	//initialize ordertemplate to be empty
		}
	}

	//CASE 1: 2edges in section
	for(i1=0; i1<=M; i1++){
		for(j1=0; j1<=L; j1++){
			ordertemplate[0][i1][j1]=1;	//loop through all possible locations for the first edge
			for(i2=0; i2<=M; i2++){
				for(j2=0; j2<=L; j2++){		//loop through all possible locations for the second edge
					if((i1!=i2) || (j1!=j2)){	//that isn't in the location where the first edge is.
						ordertemplate[0][i2][j2]=2;
						if(num_section(0)==sectionnum){
							return 1;
						}
						ordertemplate[0][i2][j2]=0;
					}
				}
			}
			ordertemplate[0][i1][j1]=0;
		}
	}


	
	return 0;
}





#include "../utils/unsgnlong_vecalloc.c" /*unsigned long int unsgnlong_vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with unsigned long integer enteries */

#include "../utils/vecalloc.c"		/*vec_ent * vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with enteries of type vec_ent */

#include "../utils/matalloc.c"		/*mat_ent **matalloc(int rowlow, int rowhigh, int collow, int colhigh)*/
/* used for dynamically allocating memory for matrices with enteries of type mat_ent */



#include "../experiments/search.c"
#include "../utils/addtowalk.c"
#include "recursivegenerate_all_interestPK.c"

	









