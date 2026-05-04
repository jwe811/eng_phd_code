/***************************************************************************/
/* To run in background type:	nohup nice +19 HingeMakerC >& outfile &	   */
/* landscape mode with 1 column type: enscript -1 -r Pmath220 filename.c   */
/* to compile with optimization type: gcc HM97.c -lm -O3 -o nameprogram
/***************************************************************************/

//FULL HINGE VERSION

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

/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDYIED **/
#define	M 1				/* number of vertical edges                                */
#define	L 2				/* number of horizontal edges                              */

#define totalspan 4			//span of generated polygon (max span of 6 as of now. too many polygons with span=7 (>2,129,088,216). need to make changes to data types)
#define force 0				/* fixed force. */

//#define	connectivity_inv 0.388795	/* When the connectivity is known its inverse is defined here.   */
						/* Doing this will save time, as the rtflsp will not need  */
						/* to be executed to determine the connectivity.		 */
						/* The connectivities (inversed) for various cases follow:       */ 
						/* (1-1) 0.547397 , (1-2) 0.437382 (+1) 0.498950, (1-3) 0.388795 */
						/* 0-? values in Klein article */

//#undef	connectivity_inv	/* Use if connectivity is not known, comment out otherwise.	 */
#define x_L 0.2			/* left end point to use in bracketing method to find x_0        */
#define x_R 0.8			/* right end point to use in bracketing method to find x_0       */

#define		SS	/* TS for two-span // SS for section-section or when conective constat is not known!*/
			/* Note: TS version will find connective constant as well, but is slower */
			/* Also need to run TS version to get data on expected number of times */
			/* a given k-span occurs in an n step polygon 			*/

#define max_tspans   963096		/* vector dimension for two-span matrix  		*/
					/* Use the following list to optimize max_tspans 	*/
					/* This optimization only works with TS mode            */
					/* When using SS mode set max_tspans = max_sections     */
					/* M:1 L:4 => 129413546 valid two-spans 		*/
					/* M:1 L:3 =>    963096 valid two-spans */
					/* M:1 L:2 =>      9702 valid two-spans */
					/* M:1 L:1 =>       108 valid two-spans */

					/* M:0 L:2 =>         7 valid two-spans    3 valid sections    3 column states*/
					/* M:0 L:3 =>        36 valid two-spans    9 valid sections    8 column states*/
					/* M:0 L:4 =>       172 valid two-spans   25 valid sections   20 column states*/
					/* M:0 L:5 =>       814 valid two-spans   74 valid sections   50 column states*/
					/* M:0 L:6 =>      3906 valid two-spans  224 valid sections  126 column states*/
					/* M:0 L:7 =>     19150 valid two-spans  711 valid sections  322 column states*/
					/* M:0 L:8 =>     96022 valid two-spans 2319 valid sections  834 column states*/
					/* M:0 L:9 =>    491699 valid two-spans      valid sections 2187 column states*/
					/* M:0 L:10 =>  2565713 valid two-spans      valid sections 5797 column states*/

/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/


#define	vM (M+1)				/* number of vertices in vertical direction   = M + 1      */
#define	vL (L+1)				/* number of vertices in horizontal direction = L + 1      */

#if M==0
#define		CS     /* define CS for column states when M=0 */
#endif

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

#elif vM*vL<11
#define maxOrdNum	699942
#define max_sections	5599530
#define vec_length	5599531

#endif

#define	makeheaderij (1 << (j + i * vL))	/* given i,j yields appropriate hnodebin */

#define		vec_ent double	/*vector entries will be of type defined here*/
#define		mat_ent double	/*matrix entries will be of type defined here*/

#define	newline		printf("\n")	/* Allows the use of newline; in place of printf("\n"); 		  */
#define VALID		1		/* This simply allows statements such as "if(something == VALID) {...}    */
#define NOTVALID 	!VALID		/* This simply allows statements such as "if(something == NOTVALID) {...} */

#define hingeedges (L + 2*M*L + M)

#if (M==1 && L==1)
#define maxendhinges 15
#elif (M==1 && L==2)
#define maxendhinges 127
#elif (M==1 && L==3)
#define maxendhinges 1023
#endif

/***************************************************************************/
/****************** Structures used in this program ************************/
/***************************************************************************/


struct hinge_span {	/* This data structure is used to store information about two-spans */

	/*unsigned long int 	inorder */		/* no longer needed as it is implicitely defined */
	unsigned long int 	outorder	:26;	/* the section number of the righthand section   */	
	unsigned long int 	edgecount	:6;	/* the number of edges in the hinge		 */
	struct hinge_span 	*nexthinge;		/* pointer to another hinge_span data structure  */
	unsigned long int	num_con_endhinges_left;	// the number of connecting endhinges that connect to the left of section "inorder"
	unsigned long int*	con_endhinges_left;		// array that holds the connecting endhinge numbers of the endhinges that connect
										// to the left of section "inorder"
	unsigned long int	num_con_endhinges_right;
	unsigned long int*	con_endhinges_right;

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

unsigned long int	num_outsections[vec_length];   /* For a given section (which has a number section_num assigned to it)*/
						       /* num_outsections[section_num] is the number of two spans with this  */
						       /* section as a first section                                         */

vec_ent			*L_Evector[2]; 			/*An array of pointers to vectors*/
							/* L_Evector[0] denotes left eigenvector 0 */
							/* L_Evector[1] denotes left eigenvector 1 */ 
vec_ent			*R_Evector[2]; 			/*An array of pointers to vectors*/
							/* R_Evector[0] denotes right eigenvector 0 */
							/* R_Evector[1] denotes right eigenvector 1 */
 
unsigned long int	*tspans_outsection[vec_length];	/*An array of pointers to integers */
							/*records the "outsection" of the two-span     */
							/* tspans_outsection[insection][nth_tspan] */
							/* nth_tspan = nth two-span with firstsection = insection */
							

unsigned long int	*tspans_edges[vec_length];	/*An array of pointers to integers*/
							/*records the number of "edges" in the two-span */
							/* example:          */
							/*tspans_edges[1][5] gives the number of edges in */
							/*the first section and hinge of the two-span         */
							/*with first section = 1			      */
							/*and is the 5th such two-span with first section = 1 */

unsigned long int	*tspans_nrr[vec_length];	/*An array of pointers to integers*/
							/*records the tspan "number"      */
							/* defined by the array index */
							/*example:		*/
							/*tspans_nrr[1][1] gives the two-span number */
							/* of the two-span with first section = 1 */
							/* and the second section is the first such */
							/* section which connects with first section = 1 */

unsigned long int	num_connections_left[vec_length];	//number of endhinges that connect to [section]

unsigned long int	*conn_endhinges_left[vec_length];	//conn_endhinges[1][2] gives the endhinge number of 2nd
							//endhinge that connects to the left of section 1

unsigned long int	num_connections_right[vec_length];
unsigned long int	*conn_endhinges_right[vec_length];


#if defined(force)
double			fval=force;			/* fixed force */
#endif
int			total_end_hinges;		/*the total number of possible endhinges	*/
unsigned short int occupied[vM][vL];	//this will indicate which verticies are occupied (i,j)
unsigned long int curindex=0;			//index of which endhinge is being created for ordertemplate (initially zero)
unsigned long int edgesin_temp = 0;	//will count the number of edges in temp[]
unsigned long int sec_polygons = 0;	//counts the number of polygons that start with a certain section.
unsigned long int tot_polygons = 0;	//used to count the total number of polygons that are of span=span
int			leftside;				//leftside=1 if testing endhinges on the left; otherwise (testing endhinges on right) leftside=0.




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

void			record_1_span(void);
unsigned short int	one_span_byte_bit(unsigned long int, unsigned long int);

void            	recordtemplate( int (*)[3]);
unsigned long int	num_section(int);
unsigned short int	paircomp(unsigned short int, unsigned short int);


/**** dominant eigenvalue-eigenvector functions  ********************************************/

double			max_eval_LRvec(double fugacity); /* Uses powermethod */
double			get_Beta(double);

/****  User input functions  ***************************************************************/

unsigned long int	getval(void);
unsigned long int	getval_mnmx(unsigned long int min, unsigned long int max);
double			enterpats(double,double,double);

void			scan_template(int side);

/****  output functions  *******************************************************************/

void			echo_template(int side);

/****  Memory allocation functions  ********************************************************/

struct hinge_span	*newhinge(void);
unsigned long int	*unsgnlong_vecalloc(int low, int high);
vec_ent			*vecalloc(int low, int high);
mat_ent			**matalloc(int rowlow, int highrow, int lowcol, int highcol);

/****  Numerical recipes in C  *************************************************************/

double 		rtflsp(double(*func)(double),double x1,double x2,double xacc);
void		qcksrt(unsigned long int,unsigned long int []);
void		qcksrtII(unsigned long int,unsigned long int [],unsigned long int []);

/*******************************************************************************************/

void			get_kspan(unsigned short int k);
unsigned short int	entsecval(unsigned short int side);
unsigned short int	entseccon(unsigned short int side);
unsigned long int	choose(unsigned long int, unsigned long int);
unsigned long int 	factorial(unsigned long int n);

struct endhinge* newendhinge(void);
void generate_endhinges(void);
void iterate(int k, unsigned short int *edgestoadd, int counter, unsigned short int totedges, unsigned long int* endnum, struct edge* edges, unsigned long int* endhingesadded);
void printendhinge(int k);	//prints endhinge k
void printtempendhinge(int k);	//prints temp[k]
unsigned long int findnumconendhingesleft(void);
unsigned long int findnumconendhingesright(void);
void connect(unsigned long int small, unsigned long int big, unsigned long int cindex);
void goup(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int uindex);
void godown(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int dindex);
void goleft(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int lindex);
void goright(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int rindex);
unsigned short int testmatch(unsigned long int i, unsigned long int j);
void generatepolys(unsigned long int curspan, unsigned long int sec_num);

void printsection(int sectionnum);


/***************************************************************************/
/******************* end of functions used in this program *****************/
/***************************************************************************/


main(void)
{
/**** variable used in main ************************************************/

	/*int	           	(*pointordNum)[3]; *//*pointer to an array of int*/
	int			ordNum[3],side = 0;
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
	void            conv_to_array(void);

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
	printf("\nThis is the Full Hinge version\n");


	clock();	/* Start clock to time program */




	//generate all endhinges
		//give endhinges a number (just use the order that they were generated)
	//perform check if endhinge connects to a sections (loop over all endhinges, and sections)
		//record which endhinges connect to the left side of each section, and similarily for the right side of each section
		//must be done BEFORE conv_array (because hinge_span structures are freed during conv_array)


	generate_endhinges();
	printf("\nALL ENDHINGES GENERATED!!!!! (there were a total of %d endhinges)\n\n", total_end_hinges);

	//prints all endhinges (check)
/*	for(i=0; i<=total_end_hinges-1; i++){
		printendhinge(i);
	}
*/
	ordNum[0] = 1;	/* first edge in section 0 gets numbered 1 (left side of 2-span)	*/
	ordNum[1] = 1;	/* first edge in section 1 gets numbered 1 (right side of 2-span)	*/
	ordNum[2] = 0;	/* to start there are no edges in the 2-span (edges in hinge)		*/

	for (i=1;i<=max_sections;i++){
		current_hinge_span[i] = newhinge();		/* free up space of hinge structure */
		first_hinge_span[i] = current_hinge_span[i];	/* make first_hinge_span[i] point to this first one */
	}
	printf("space created for hinge structures\n");
	fflush(stdout);
	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
//			printf("entering hinge at i=%d, j=%d\n", i, j);
			enterhinge(makeheaderij, makeheaderij, i, j, side, &ordNum);
		}
	}

/*************************************************************************************************/
	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
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


	printf("converting linked list to arrays...\n");
	conv_to_array();
	printf("done converting to arrays.\n\n");


	int previousi=0;
	int previous=0;
	for(i=1;i<=max_sections;i++){
		if(num_outsections[i]>0){
			for(j=1; j<=num_outsections[i]; j++){
				if(previous==tspans_outsection[i][j] && i==previousi ){
					printf("(prev) section %d connects to %d, (tspan #%d)\n", previousi, tspans_outsection[previousi][j-1], tspans_nrr[previousi][j-1]);
					printf("section %d connects to %d, (tspan #%d)\n", i, tspans_outsection[i][j], tspans_nrr[i][j]);
					printsection(previousi);
					printsection(tspans_outsection[previousi][j-1]);
				}
				previous=tspans_outsection[i][j];
				previousi=i;
			}


/*			if(num_connections_left[i]>0){
				printf("num_connections_left[%d]=%d, ", i, num_connections_left[i]);
				for(j=1; j<=num_connections_left[i];j++){
					printf("conn_endhinges_left[%d][%d]=%d, ", i, j, conn_endhinges_left[i][j]);
				}
				printf("\n");
			}
			if(num_connections_right[i]>0){
				printf("num_connections_right[%d]=%d, ", i, num_connections_right[i]);
				for(j=1; j<=num_connections_right[i];j++){
					printf("conn_endhinges_right[%d][%d]=%d, ", i, j, conn_endhinges_right[i][j]);
				}
				printf("\n");
			}
*/		}
		else{
//			printf("section %d has no outsections (not a valid section number)\n", i);
		}
	}

	//count the number of polygons of span=span
	unsigned long int mult_starthinge;	//#of starthinges that connect to the left of the first section

	for(i=1; i<=max_sections; i++){		//for each section
		if(num_outsections[i]>0){
			mult_starthinge = num_connections_left[i];
//			printf("looking for polygons that start with section %d\n", i);
			sec_polygons=0;
			generatepolys(1, i);	//generate all polygons that have section i as the first section.
			tot_polygons += sec_polygons*mult_starthinge;	//takes into account the different starthinges.
			if(mult_starthinge!=0){
//				printf("found %d polygons that start with section %d. Multiply this by %d for the different starthinges, and add this to total. we now have %d total polygons\n", sec_polygons, i, mult_starthinge, tot_polygons);
			}
		}
	}
	printf("\nTOTAL NUMBER OF POLYGONS THAT ARE IN A M=%d, L=%d tube, of span=%d are %d!! (Fullhinge version)\n\n", M, L, totalspan, tot_polygons);

	exit(1);

//#if defined(TS)
//	L_Evector[0] = vecalloc(0,max_tspans); /*global*/
//	L_Evector[1] = vecalloc(0,max_tspans); /*global*/ /*temp*/
//	R_Evector[0] = vecalloc(0,max_tspans); /*global*/
//	R_Evector[1] = vecalloc(0,max_tspans); /*global*/ /*temp*/
//#else
//	L_Evector[0] = vecalloc(0,max_sections); /*global*/
//	L_Evector[1] = vecalloc(0,max_sections); /*global*/ /*temp*/
//	R_Evector[0] = vecalloc(0,max_sections); /*global*/
//	R_Evector[1] = vecalloc(0,max_sections); /*global*/ /*temp*/
//#endif

//#if defined(connectivity_inv)
//	max_eval_LRvec(connectivity_inv); /*rtflsp not necessary once connectivity is known*/
//#else
//	printf("Starting false position method to determine connective constant.\n");
//	#if defined(force)
//		printf("\n=> %f is fugacity which yields 1 as maximum eigenvalue (after multiplying matrix by e^%f)\n",
//					connectivity_inv=rtflsp(&max_eval_LRvec,x_L,x_R, 0.0000001), fval);
//	#else
//		printf("\n=> %f is fugacity which yields 1 as maximum eigenvalue\n",
//					connectivity_inv=rtflsp(&max_eval_LRvec,x_L,x_R, 0.0000001));
//	#endif
//		printf("\n=> %f is the connective constant.\n",1/connectivity_inv);
//#endif

//	LmultR = 0;
//#if defined(TS)
//	for (i=1;i<= max_tspans;i++){
//		LmultR += (L_Evector[0][i])*(R_Evector[0][i]);
//	}
//#else
//	for (i=1;i<= max_sections;i++){
//		LmultR += (L_Evector[0][i])*(R_Evector[0][i]);
//	}
//#endif

	/* LmultR = L_Evector*R_Evector *** use this for normalizing the vectors */
//	Ephi_Ephi = (L_Evector[0][1])*(R_Evector[0][1])/LmultR;
//	Beta = get_Beta(connectivity_inv);

//#if defined(TS)
//	Alpha = connectivity_inv*connectivity_inv*Ephi_Ephi/Beta;
//#else /* defined(SS) */
//	Alpha = connectivity_inv*connectivity_inv*connectivity_inv*connectivity_inv*Ephi_Ephi/Beta;
//#endif

	/*may now free L_Evector[1] this has yet to be implemented */
	/*may now free R_Evector[1] this has yet to be implemented */

//	newline;
//	printf("Alpha = %f = %e\n",Alpha,Alpha);
//	printf("Beta  = %f = %e\n\n",Beta,Beta);
//	printf("L_Evector[phi]*R_Evector[phi] = %f = %e\n\n",Ephi_Ephi,Ephi_Ephi);

/************************************************************************************************/	
/*showmenu*/ /*Needs to be implemented*/

//	kappa = connectivity_inv;
	/*system("clear");*/
//	get_kspan(4);
//	Enum_Enum = enterpats(Beta, LmultR, kappa);
//	printf("Enum_Enum=%f", Enum_Enum);
	system("date");  /* prints the date and time */
	printf("\nProgram Complete.\n\n");


}				/* end of main */
/***************************************************************************/


/***** start of function conv_to_array *************************************/
/* called only by main */
void
conv_to_array(void) /* tspans_outsection[section_num][1..num_outsections[section_num]] */
{
	unsigned long int i;
	unsigned long int section_num,tspan_num=0,nth_outsection;
	unsigned long int arr_ent;	/*array entry*/
	struct hinge_span	*hinge_to_free;

	for (section_num = 1; section_num <= max_sections; section_num++) {

		num_outsections[section_num] = 0;	/* initialize to zero */
		current_hinge_span[section_num] = first_hinge_span[section_num];
		/* this causes current_hinge_span to point to the first element in */
		/* the various linked lists					   */

		while ((*current_hinge_span[section_num]).nexthinge != NULL) {
			current_hinge_span[section_num] = (*current_hinge_span[section_num]).nexthinge;
			num_connections_left[section_num] = current_hinge_span[section_num]->num_con_endhinges_left;
			num_connections_right[section_num] = current_hinge_span[section_num]->num_con_endhinges_right;
			num_outsections[section_num]++; /*count that number of two-spans with section_num as first section*/

		} /* this while loop counts the number of valid two spans with first section = section_num */

		tspans_outsection[section_num] = unsgnlong_vecalloc(1, num_outsections[section_num]);/*allocate memory*/

		tspans_edges[section_num] = unsgnlong_vecalloc(1, num_outsections[section_num]);/*allocate memory e=edges*/
							/*for recording the number of edges in the two-spans */
							/*with first section = section_num    */

		tspans_nrr[section_num] = unsgnlong_vecalloc(1, num_outsections[section_num]);/*allocate memory n=number*/
							/*for recording the two-spans number of a two-span which has */
							/* first section = section_num    */
		conn_endhinges_left[section_num] = unsgnlong_vecalloc(1, num_connections_left[section_num]);
		conn_endhinges_right[section_num] = unsgnlong_vecalloc(1, num_connections_right[section_num]);

	}
	for (section_num = 1; section_num <= max_sections; section_num++) {
		arr_ent = 1;	/* initialize to zero */
		current_hinge_span[section_num] = first_hinge_span[section_num];
		for(i=1; i<=num_connections_left[section_num]; i++){
			conn_endhinges_left[section_num][i] = current_hinge_span[section_num]->nexthinge->con_endhinges_left[i-1];
		}
		for(i=1; i<=num_connections_right[section_num]; i++){
			conn_endhinges_right[section_num][i] = current_hinge_span[section_num]->nexthinge->con_endhinges_right[i-1];
		}
		while ((*current_hinge_span[section_num]).nexthinge != NULL) {
			hinge_to_free = current_hinge_span[section_num];
			current_hinge_span[section_num] = (*current_hinge_span[section_num]).nexthinge;
			free(hinge_to_free);
			tspans_outsection[section_num][arr_ent] = (*current_hinge_span[section_num]).outorder;
			tspans_edges[section_num][arr_ent] = (*current_hinge_span[section_num]).edgecount;
			arr_ent++;
		}
	}
	for (section_num = 1; section_num <= max_sections; section_num++) {
		qcksrtII(num_outsections[section_num], tspans_outsection[section_num],tspans_edges[section_num]);
	}

	for (section_num=1; section_num<=max_sections;section_num++){
		for (nth_outsection=1;nth_outsection<=num_outsections[section_num];nth_outsection++){
			tspans_nrr[section_num][nth_outsection] = ++tspan_num;
		}
	}
			
return;
}

/***************************************************************************/

void
enterhinge(bits_16 occstatebin, bits_16 hnodebin, int i, int j, int side, int (*pointordNum)[3])
{

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
	return;

}				/* end of function enterhinge */

/***************************************************************************/

void
leavehinge(bits_16 occstatebin, bits_16 hnodebin, int i, int j, int side, int (*pointordNum)[3])
/* called by: rowedges,coledges,enterhinge */
{
	int             ii;		/*vertex number in vertical direction*/
	int		jj;		/*vertex number in horizontal direction */
	int		shift = 0;	/*keeps track of how far the headernode is shifted */
					/*from the 0,0 position going left to right first */
					/*then down as each row of vertices is exhausted  */
	int fullhinge=1;
	int loopi;
	int loopj;

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

			for(loopi=0; loopi<=M; loopi++){
				for(loopj=0; loopj<=L; loopj++){
					if((occstatebin >> loopj + loopi * (vL) & 1) == 0){	// if vertex i,j is empty
					/* the expression (occstatebin >> j + i * (vL)) & 1  evaluates to 0 if the i,j vertex is empty */
					/* otherwise it evaluates to 1 */
//						printf("found a 2span, but not all verticies are filled\n");
						fullhinge=0;
					}
				}
			}

			if(fullhinge==1){
				valid_2_spans++; /*This is a valid 2 span so it is counted*/
				record_1_span(); /*The component one spans are therefore valid so need to be recorded*/
				recordtemplate(pointordNum); /*the information contained in ordertemplate is recorded*/
			}
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
				enterhinge(occstatebin + (1 << (shift)), 1 << (shift), ii, jj, side, pointordNum);
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
void
recordtemplate( int (*pointordNum)[3])
/* This function records the pertinent information about a two-span */
/* which has been verified as being valid into a linked list */
{
	unsigned long int i;
	unsigned long int j;
	unsigned long int	inNum = num_section(0);

/*if(inNum==4){
	for (i = 0; i <= M; i++) {
		printf("\n");
		for (j = 0; j <= L; j++)
			printf("%i ", ordertemplate[0][i][j]);
	}
	printf("\n\n");
}
*/

	unsigned short int	Num_edges= (*pointordNum)[0] + (*pointordNum)[2] - 1;
				/*number of edges = edges in first section + edges in hinge */
				/*the correction (-1) is necessary since (*pointordNum)[0] */
				/*is the number that the next edge in the first section will */
				/*get */

	(*current_hinge_span[inNum]).nexthinge = newhinge(); /* add a newhinge to the linked list */
	
	current_hinge_span[inNum] = (*current_hinge_span[inNum]).nexthinge; /* let current_hinge_span[inNum] point to newhinge */
	(*current_hinge_span[inNum]).outorder = num_section(1);	/* record the ordering on the out side */
	(*current_hinge_span[inNum]).edgecount = Num_edges;	/* record the number of edges in first section and in hinge */

	//calculate num_con_endhinges_left
	//find which endhinges connect to left side of ordertemplate[0] (left side)
	//record their "number" in con_edges_left[num_con_endhinges_left]

	//Then do the same for endhinges that connect to the right side of the section inNum

//	printf("recorded template. left side is: ");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
//			printf("ordertemplate[0][%d][%d] = %d, ", i, j, ordertemplate[0][i][j]);
		}
	}

//	printf("(SECTION NUMBER: %d)\n", inNum);

	unsigned long int num_con_endhinges_left = findnumconendhingesleft();	//finds the number of endhinges that ordertemplate[0] connects to
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
//		printf("the endhinges that connect to the left of that ordertemplate are:\n");
		for(i=0; i<=num_con_endhinges_left-1; i++){
//			printendhinge(con_edges_left[i]);
		}
	
		current_hinge_span[inNum]->num_con_endhinges_left = num_con_endhinges_left;	//number of endhinges that connect to the left of section inNum
		current_hinge_span[inNum]->con_endhinges_left = unsgnlong_vecalloc(0, num_con_endhinges_left-1);	//list of endhinge "numbers" that connect to the left of section inNum
		for(i=0; i<=num_con_endhinges_left-1; i++){
			current_hinge_span[inNum]->con_endhinges_left[i] = con_edges_left[i];
		}
	}
	curindex=0;	//reset curindex
	edgesin_temp=0;	//reset edgesin_temp
	//freetemp





//	DO THE SAME THING FOR THE RIGHT SIDE OF SECTION INNUM	/////////////////////////////////////////////




	unsigned long int num_con_endhinges_right = findnumconendhingesright();	//finds the number of endhinges that ordertemplate[0] connects to
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
	
		//so con_edges_right contains the indicies of endhinge[] that connect to the right of ordertemplate[0]
//		printf("the endhinges that connect to the right of that ordertemplate are:\n");
		for(i=0; i<=num_con_endhinges_right-1; i++){
//			printendhinge(con_edges_right[i]);
		}
	
		current_hinge_span[inNum]->num_con_endhinges_right = num_con_endhinges_right;	//number of endhinges that connect to the right of section inNum
		current_hinge_span[inNum]->con_endhinges_right = unsgnlong_vecalloc(0, num_con_endhinges_right-1);	//list of endhinge "numbers" that connect to the right of section inNum
		for(i=0; i<=num_con_endhinges_right-1; i++){
			current_hinge_span[inNum]->con_endhinges_right[i] = con_edges_right[i];
		}
	}
	curindex=0;	//reset curindex
	edgesin_temp=0;	//reset edgesin_temp
	//freetemp



	return;
}

#if defined(TS)
	#if defined(force)
		#include "../transfer_matrix/pw_meth_ts_LRvec_fcheck.c" /* for applying the power method to the 2-span 2-span matrix */
	#else
		#include "../transfer_matrix/pw_meth_ts_LRvec.c" /* for applying the power method to the 2-span 2-span matrix */ 
	#endif
	#include "../analysis/statistics/get_Beta_ts.c"      /* for calculating Beta */
#endif

#if defined(SS)
	#if defined(force)
		#include "../transfer_matrix/pw_meth_ss_LRvec_fcheck.c" /* for applying the power method to the section-section matrix */
	#else
		#include "../transfer_matrix/pw_meth_ss_LRvec.c" /* for applying the power method to the section-section matrix */
	#endif
	#include "../analysis/statistics/get_Beta_ss.c"	/* for calculating Beta */
#endif

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

#include "../sections/printsection.c"

#endif







/***************************************************************************/

unsigned long int
getval(void)
/* gets an integer from the user */
{
	char            input[100];
	return (atoi(gets(input)));
}

unsigned long int
getval_mnmx(unsigned long int min, unsigned long int max)
/* get unsigned long int from user between the values min and max */
{
	char            input[100];
	unsigned long int x;
	gets(input);
	x = atoi(input);
	if ((min <= x) && (x <= max)) {
		return (x);
	} else {
		return (min);	/* if the value isn't in the range, return the min */
	}
}

void
scan_template(int side)
{
/* allows the user to enter the edges in a section of a k-span */
unsigned short int i,j;

	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++)
			scanf("%i",&(ordertemplate[side][i][j]));
	}
}

void
echo_template(int side)
{
/* print the edges in a section to the user */
unsigned short int i,j;

	for (i = 0; i <= M; i++) {
		printf("\n");
		for (j = 0; j <= L; j++)
			printf("%i ", ordertemplate[side][i][j]);
	}

}

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
	return nextnewhinge;
}
/***************************************************************************/
/*THE FOLLOWING WAS ADDED MAY 26, 1998 */

/* actually need to return two-span # */
/* or first section # and last section # */
/* and total number of edges */
void
get_kspan(unsigned short int k)
{

	unsigned short int n = 0;
	unsigned short int valid=0;
	unsigned long int tspan_1_num;
	unsigned long int tspan_k_num;
	unsigned long int sec_1_num;
	unsigned long int sec_k_num;
	unsigned long int totedges=0;
	unsigned short int side=0;

	printf("The number of sections to be entered is: %i\n",k);
	printf("Enter %i row(s) of %i numbers for each section\n\n", vM, vL);

	n++;
	do	/* get first section */
	{	
		printf("\nEnter section #1:\n\n");
		scan_template(0);	/* enter section #n */
		valid = entsecval(0);	/* validity test, but not necessary to test connecting */
	} while (valid == 0);

	sec_1_num = num_section(0);

	while (n < k) {

		n++;

		side ^= 1; /*switch sides*/

		do	/* get number of edges in hinge */
		{
			printf("\nEnter number of edges in hinge:");
			scanf("%i",&(totedges));
			getval_mnmx(1,100);
			/*scanf("%i", &(totedges));*/
			printf("\ntotedges are:%i\n",totedges);
			valid=1;
		} while (valid == 0);

		do	/* get next section */
		{
			printf("\nEnter section #%i.\n\n", n);
			scan_template(side);	/* enter section #n */
			valid = entsecval(side) * entseccon(side^1);	/* valid and connecting */
		} while (valid == 0);
	}

	sec_k_num = num_section(side);
}

unsigned short int
entsecval(unsigned short int side) /* entered section validity test */
{
	unsigned long int section_num = num_section(side);

	if (1 <= section_num && section_num <= max_sections) {	/* ensure in range */
		if (num_outsections[section_num] != 0) {	/* ensure it connects to
								 * other sections */
		return(1); /*entered section is valid*/
		} else {
			printf("This section does not connect to anything.\n");
		}
	} else {
		printf("The entered section is not valid!\n");
	}
	return(0); /*The entered section is not valid.*/
}


unsigned short int		/* entered section connects test */
entseccon(unsigned short int side)
{

	unsigned long int section_1 = num_section(side);
	unsigned long int section_2 = num_section(side ^ 1);
	unsigned short int result = 0;
	unsigned long int nth_tspan;

	for (nth_tspan = 1; nth_tspan <= num_outsections[section_1]; nth_tspan++) {
		if (section_2 == tspans_outsection[section_1][nth_tspan]) {
			printf("\n%6u:%6u : edges=%3u : two span# =%3u\n", \
			       section_1, tspans_outsection[section_1][nth_tspan], \
			       tspans_edges[section_1][nth_tspan], \
			       tspans_nrr[section_1][nth_tspan]);
			result = 1;
		}
	}
	if (result == 0) {
		printf("This section does not connect to the previous section.\n");
	}
	return (result);
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

unsigned long int choose(unsigned long int n, unsigned long int r) {
	if(r>n) {
		return 0;
	}
	return factorial(n)/(factorial(r)*factorial(n-r));
}

void generate_endhinges(void){
	unsigned short int i, j, k;		//looping variables
	unsigned short int index=0;		//initially zero

	struct edge edges[hingeedges];		//array of edges.
	//give each edge in the hinge a number (0 to (hingeedges-1))

	//label all of the horizontal edges first (left to right, move down, left to right, move down, etc....)

	for(j=0; j<=M; j++){	//for each row
		for(i=0; i<=L-1; i++){		//for each edge in the row
			edges[index].f.x = i;
			edges[index].f.y = j;
			edges[index].l.x = i+1;
			edges[index].l.y = j;
			index++;
		}
	}

	//next, label all of the vertical edges (left to right, move down, left to right, move down, etc....)
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
	for(i=0; i<=index-1; i++){
//		printf("EDGE %d: (%d,%d) to (%d,%d)\n", i, edges[i].f.x, edges[i].f.y, edges[i].l.x, edges[i].l.y);
	}


	//use this to easily add edges to endhinge....
	//need an algorithm that goes through all combinations of edges in the end hinge.

	//FULL HINGE VERSION:
	
	unsigned short int totedges=1;		//total edges in the endhinge. loop from (1 to hingeedges)
	unsigned long int endnum=0;		//the number/index for the current endhinge
	unsigned long int endhingestoadd;	//given totedges, endhingestoadd is the number of possible combinations of endhinges that need to be added
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

			//FULL HINGE CHANGE
			//check if all of the edges in edges[edgestoadd[]] fill all of the verticies in the hinge

			unsigned short int filledhinge=1;
			unsigned short int filledvertex=0;
			for(z=0; z<=M; z++){	//loop through all verticies
				for(a=0; a<=L; a++){
					filledvertex=0;
					for(j=0; j<=totedges-1; j++){	//look through all edges all edges
						if((edges[edgestoadd[j]].f.x==a || edges[edgestoadd[j]].l.x==a) && (edges[edgestoadd[j]].f.y==z || edges[edgestoadd[j]].l.y==z)){//if vertex (z,a) is occupied by edges[edgestoadd[j]]
							filledvertex=1;
							break;
						}
					}
					if(filledvertex==0){	//vertex (z,a) isn't filled by any edge in edgestoadd[]
						filledhinge=0;
					}
				}
			}

			if(filledhinge==1){
				//add the totedges to a endhinge structure
				firstendhinge[*endnum] = newendhinge();
				currentendhinge[*endnum] = firstendhinge[*endnum];
				for(j=0; j<=totedges-1; j++){
					//add edgestoadd[j] to endhinge
					currentendhinge[*endnum]->edge.f.x = edges[edgestoadd[j]].f.x;
					currentendhinge[*endnum]->edge.f.y = edges[edgestoadd[j]].f.y;
					currentendhinge[*endnum]->edge.l.x = edges[edgestoadd[j]].l.x;
					currentendhinge[*endnum]->edge.l.y = edges[edgestoadd[j]].l.y;
//					printf("just added edge %d to endhinge #%d\n", edgestoadd[j], *endnum);
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
	// finds and counts how many endhinges connect to the right of ordertemplate[0] ("left" side of ordertemplate)
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

	//adds all endhinges that connect point where small is to point where big is.
//	printf("CONNECT called. connecting edge %d to %d. cindex=%d\n", small, big, cindex);

	unsigned long int i, j;
	struct point start;	//starting point (where small is)
	struct point end;	//end point (where big is)

	//assign the proper x and y values for the start and end points
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			if(ordertemplate[0][i][j] == small){
				start.x = j;
				start.y = i;
			}
			else if(ordertemplate[0][i][j] == big){
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

			//NOW CHECK IF OCCUPIED IS FULL (FULLVERSION)
			int filled=1; //1 if filled, 0 if empty vertex
			int i, j;
			for(i=0; i<=M; i++){
				for(j=0; j<=L; j++){
					if(occupied[i][j]==0){	//NOT FULL
						filled=0;
						break;
					}
				}
			}

			if(filled==1){
//				printf("no more points need to be connected. increase curindex from %d to %d. temp[%d] is now a valid endhinge\n", curindex, curindex+1, curindex);
//				printtempendhinge(curindex);
				curindex++;			//finished temp[curindex], look for other endhinges that connect ordertemplate.
			}
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

			//NOW CHECK IF OCCUPIED IS FULL (FULLVERSION)
			int filled=1; //1 if filled, 0 if empty vertex
			int i, j;
			for(i=0; i<=M; i++){
				for(j=0; j<=L; j++){
					if(occupied[i][j]==0){	//NOT FULL
						filled=0;
						break;
					}
				}
			}

			if(filled==1){
//				printf("no more points need to be connected. increase curindex from %d to %d. temp[%d] is now a valid endhinge\n", curindex, curindex+1, curindex);
//				printtempendhinge(curindex);
				curindex++;			//finished temp[curindex], look for other endhinges that connect ordertemplate.
			}

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

			//NOW CHECK IF OCCUPIED IS FULL (FULLVERSION)
			int filled=1; //1 if filled, 0 if empty vertex
			int i, j;
			for(i=0; i<=M; i++){
				for(j=0; j<=L; j++){
					if(occupied[i][j]==0){	//NOT FULL
						filled=0;
						break;
					}
				}
			}

			if(filled==1){
//				printf("no more points need to be connected. increase curindex from %d to %d. temp[%d] is now a valid endhinge\n", curindex, curindex+1, curindex);
//				printtempendhinge(curindex);
				curindex++;			//finished temp[curindex], look for other endhinges that connect ordertemplate.
			}

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

			//NOW CHECK IF OCCUPIED IS FULL (FULLVERSION)
			int filled=1; //1 if filled, 0 if empty vertex
			int i, j;
			for(i=0; i<=M; i++){
				for(j=0; j<=L; j++){
					if(occupied[i][j]==0){	//NOT FULL
						filled=0;
						break;
					}
				}
			}

			if(filled==1){
//				printf("no more points need to be connected. increase curindex from %d to %d. temp[%d] is now a valid endhinge\n", curindex, curindex+1, curindex);
//				printtempendhinge(curindex);
				curindex++;			//finished temp[curindex], look for other endhinges that connect ordertemplate.
			}

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
unsigned short int testmatch(unsigned long int i, unsigned long int j){
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


void generatepolys(unsigned long int curspan, unsigned long int sec_num){
	unsigned long int i;
//	printf("curspan=%d, sec_num=%d\n", curspan, sec_num);
	if(curspan<totalspan){	//not at the end yet
		for(i=1; i<=num_outsections[sec_num]; i++){
			generatepolys(curspan+1, tspans_outsection[sec_num][i]);
		}
	}
	else{	//reached the appropriate span. add the number of endhinges that connect to the right of that last section.
		//if num_connections_right=0, then sec_num isn't a valid ending section, so we add zero polygons to sec_polygons
		sec_polygons += num_connections_right[sec_num];
		if(num_connections_right[sec_num]>0){
//			printf("%d is the final section. added %d polygons to sec_polygons. sec_polygons now is %d\n", sec_num, num_connections_right[sec_num], sec_polygons);
		}
	}
}


#include "../sections/enterpats.c"

#include "../utils/unsgnlong_vecalloc.c" /*unsigned long int unsgnlong_vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with unsigned long integer enteries */

#include "../utils/vecalloc.c"		/*vec_ent * vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with enteries of type vec_ent */

#include "../utils/matalloc.c"		/*mat_ent **matalloc(int rowlow, int rowhigh, int collow, int colhigh)*/
/* used for dynamically allocating memory for matrices with enteries of type mat_ent */


#include "../analysis/statistics/rtflsp.c"	/*double rtflsp(double(*func)(double),double x1,double x2,double xacc)*/
			/*From Numerical Recipes in C section 9.2*/
			/*This uses the regula falsi, (false position method) */



#undef M	//qcksrt uses a different M

#include "../utils/qcksrt_int.c"	/*void qcksrt(unsigned long int n,unsigned long int arr[])*/
			/*From Numerical Recipes in C section 8.2*/
			/*Sorts an array arr[1..n] into ascending numerical order using */
			/*the Quicksort algorithm. n is input; arr is replaced on output */
			/*by its sorted rearrangement.*/

#include "../utils/qcksrtII_int.c"	/*void qcksrt(unsigned long int n,unsigned long int arr[],double brr[])*/
				/* sorts array arr[] and makes corresponding changes in array brr[] */






	









