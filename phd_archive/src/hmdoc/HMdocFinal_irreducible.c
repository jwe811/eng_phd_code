/*
2-spans from irreducible polygons only


COMMENTS BY JEREMY ENG, APRIL 26, 2013

This is the heavily editted and most current, up to date version of Allan Duffy's HMdoc.c.
Some changes that have been made:
 1.	Force has been implemented into the model. Basically, the transfer matrix is multiplied by e^f.
	This results in finding the x_0 that results in a dominant eigenvalue of e^-f, instead of 1.
	This change was made where the power method is used (pw_meth_ss_LRvec_fcheck.c and
	pw_meth_ts_LRvec_fcheck.c).  For the rest of the calculations, whenever you use x_0, you must
	multiply it by e^f. This has been changed for finding beta (get_Beta_ssCheck.c and get_Beta_tsCheck.c),
	and also for finding alpha (in function main()).
2.	Added a check that if you sum over all 2-spans of the expected number of occurrences (without dividing
	by beta), then you get one. This is done with the function sumtoone(), stored in sumtoone_ss.c and 
	sumtoone_ts.c.
3.	Added a check that checks that the transition probabilites for each 2-span/section add up to one. This
	is done in the function transprobcheck(), stored in transprobcheck_ss.c and transprobcheck_ts.c.
	For the 2-span 2-span version, it looks at each 2-span and then sums the transition probabilites over
	every connecting 2-span. For the section-section version, it looks at each section, and then sums the
	transition probabilities over every connecting section.
4.	The function enterpats() now works properly for both versions. Basically enterpats() allows the user
	to enter specific k-span, and it calculates the expected number of occurrences of that specific k-span.
	The function enterpats() is done in enterpats_ss.c and enterpats_ts.c
5.	Implemented automatic setting of some defined variables.  maxOrdNum, max_sections, vec_length, and CS all
	set themselves automatically, depending on the prism size.  If M==1, max_tspans also sets automatically.

Things that you can change (for different system):
M - number of vertical edges in prism
L - number of horizontal edges in prism
force - force in the model
SS or TS - section-section version or 2-span 2-span version. Both should give the same results

Things you must change manually:
x_L and x_R - used in the root finding. x_0 must lie in between x_L and x_R.
In function main(), can uncomment or comment out the functions:
	sumtoone()
	transprobcheck()
	enterpats()

END OF COMMENTS BY JEREMY ENG
*/


/***************************************************************************/
/* To run in background type:	nohup nice +19 HingeMakerC >& outfile &	   */
/* landscape mode with 1 column type: enscript -1 -r Pmath220 filename.c   */
/* to compile with optimization type: gcc HM97.c -lm -O3 -o nameprogram
/***************************************************************************/

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

/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/
//IF LOOKING AT SLIT, USE M=0
#define	L 2				/* number of horizontal edges                              */
#define	M 1				/* number of vertical edges                                */


#define force 0			/* fixed force in model. */

//#define	connectivity_inv 0.163125	/* When the connectivity is known its inverse is defined here.   */
						/* Doing this will save time, as the rtflsp will not need  */
						/* to be executed to determine the connectivity.		 */
						/* The connectivities (inversed) for various cases follow: (F=0)       */ 
						/* (1-1) 0.547397 , (1-2) 0.437382 (+1) 0.498950, (1-3) 0.388795 */
						/* 0-? values in Klein article */

//#undef	connectivity_inv	/* Use if connectivity is not known, comment out otherwise.	 */
#define x_L 0.4		/* left end point to use in bracketing method to find x_0        */
#define x_R 0.7			/* right end point to use in bracketing method to find x_0       */

#define		TS			/* TS for two-span // SS for section-section or when connective constant is not known!*/
						/* Note: TS version will find connective constant as well, but is slower */
						/* Also need to run TS version to get data on expected number of times */
						/* a given k-span occurs in an n step polygon 			*/
						/* SS version now can also be used to calculate expected number of times */
						/* a given k-span occurs in an n-step polygon	*/

//Only have to set max_tspans if M==0, otherwise it'll be set automatically
//#define max_tspans   9702		/* vector dimension for two-span matrix  		*/
					/* Use the following list to optimize max_tspans 	*/
					/* This optimization only works with TS mode            */
					/* When using SS mode set max_tspans = max_sections     */
					/* M:1 L:4 => 129413546 valid two-spans */
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

//Automatically set max_tspans
#if M==1 && L==1
#define max_tspans	108

#elif M==1 && L==2
#define max_tspans	9702

#elif M==1 && L==3
#define max_tspans	963096

#elif M==1 && L==4
#define max_tspans	129413546

#elif M==2 && L==2
#define max_tspans	12095392
#endif

#define	vM (M+1)				/* number of vertices in vertical direction   = M + 1      */
#define	vL (L+1)				/* number of vertices in horizontal direction = L + 1      */

#if M==0
#define		CS     /* define CS for column states when M=0 */
#endif

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


/***************************************************************************/
/****************** Structures used in this program ************************/
/***************************************************************************/


struct hinge_span {	/* This data structure is used to store information about two-spans */

	/*unsigned long int 	inorder */		/* no longer needed as it is implicitely defined */
	unsigned long int 	outorder	:26;	/* the section number of the righthand section   */	
	unsigned long int 	edgecount	:6;	/* the number of edges in the hinge and first section		 */
	struct hinge_span 	*nexthinge;		/* pointer to another hinge_span data structure  */

	int			walk[vM*vL*8];	/* this will hold the random walks that the 2span creates (walks are seperated by -222 or -333)	*/
	int			startpoint[2];	//contains the (x,y) coordinates of where the walk is starting (we know z=0)

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

int			curwalk[vM*vL*8];	/* this will hold the current walks, to be recorded when the 2span is complete. (-222 or -333 separates walks)*/
int			startingx=0, startingy=0;	//starting coordinates for x and y

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

int			**tspans_walk[vec_length];


#if defined(force)
double			fval=force;			/* fixed force */
#endif

int recordedsections[max_sections];
int failedsections[max_sections];
int twoedgesections=0;
int fouredgesections=0;
int sixedgesections=0;

int twoedgethrown=0;
int fouredgethrown=0;
int sixedgethrown=0;

unsigned long int* ConnectInfo[max_tspans];
unsigned long int num_outtspans[max_tspans];

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
unsigned long int	choose(unsigned short int, unsigned short int);

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
void		qcksrtII(unsigned long int,unsigned long int [],unsigned long int [], int* []);

/*******************************************************************************************/

void			get_kspan(unsigned short int k);
unsigned short int	entsecval(unsigned short int side);
unsigned short int	entseccon(unsigned short int side);

/**** Added Functions  ********************************************************************/
double			sumtoone(double,double,double);
double			transprobcheck(double, double, double);

#if defined(TS)
void				groupings(double, double, double);
double			expoccur(double, double, double, int);
void				group22(double, double, double);
void				group145(double, double, double);
void				group330(double, double, double);
void				graphexpoccurvsedges(double, double, double);
void				graphexpoccurvsnumconnect(double, double, double);
void				graphexpoccurvsconnectto(double, double, double);
void				graphexpoccurvstotalconnect(double, double, double);
void				graphedgesvsnumconnect(double, double, double);
void				graphedgesvsconnectto(double, double, double);
void				graphedgesvstotalconnect(double, double, double);
#endif

void				findsection(double, double, double);
void				printsection(int);
void				probkspan(double, double, double);
void				probkspan2(double, double, double);

int hookedhinge(int hingewalk[]);
int allstraight();
void addtowalk(int printablewalk[], int nextwalk[]); /* This function nextwalk to printablewalk

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
	for(i=0; i<=vM*vL*8-1; i++){
		curwalk[i]=0;
	}

	for(i=0; i<=max_sections; i++){
		recordedsections[i]=0;
		failedsections[i]=0;
	}

	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			printf("entering at i=%d,j=%d\n", i, j);
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

/*
	//print some walks
	i=1;
	if(num_outsections[i]>0){
		for(j=1; j<=num_outsections[i]; j++){
			printf("print 2-span made up of sections %d and %d:", i, tspans_outsection[i][j]);
			k=0;
			while(tspans_walk[i][j][k]!=0){
				printf("%d ", tspans_walk[i][j][k]);
				k++;
			}
			printf("\n");
		}
	}
*/

	int totalsecs=0;
	int total2spans=0;
	for(i=1; i<=max_sections; i++){
		if(num_outsections[i]>0){
			totalsecs++;
			total2spans=total2spans + num_outsections[i];
		}
	}
	printf("Check: Total sections=%d, total 2spans=%d\n", totalsecs, total2spans);


//Full version runs check here, for irreduciblity and "actual" valid 2span/sections counts



//////////////////////////////////////////////////////////////////////////////////////
// Check irreducibility. Can i get from any 2-span to any other 2-span ///////////////
//////////////////////////////////////////////////////////////////////////////////////

	//need some way to store information about 3-spans
	//maybe array of size=total2spans. would contain array of which tspans it can connect to.

//	int* ConnectInfo[total2spans+1];	made this global
	int num_connects;
	int threespanwalk[vM*vL*8];
	int l;
	int invalidtspans[total2spans];	//will record 2spans which do not connect to anything

	for(i=0; i<=total2spans; i++){
		invalidtspans[i]=0;
	}

	for(i=1; i<=max_sections; i++){
		for(j=1; j<=num_outsections[i]; j++){
			num_connects=0;
			for(k=1; k<=num_outsections[tspans_outsection[i][j]]; k++){
				num_connects++;
			}
			ConnectInfo[tspans_nrr[i][j]] = unsgnlong_vecalloc(0,num_connects);
		}
	}	

	int index;
	for(i=1; i<=total2spans; i++){
		num_outtspans[i]=0;
	}
	
	for(i=1; i<=max_sections; i++){
		for(j=1; j<=num_outsections[i]; j++){
			index=1;
			for(k=1; k<=num_outsections[tspans_outsection[i][j]]; k++){
				//combine walks and test if there's a hook
//				printf("Examining tspan %d to tspan %d\n", tspans_nrr[i][j], tspans_nrr[tspans_outsection[i][j]][k]);
				for(l=0; l<=vM*vL*8-1; l++){
					threespanwalk[l]=tspans_walk[i][j][l];
				}
				addtowalk(threespanwalk, tspans_walk[tspans_outsection[i][j]][k]);
				if(hookedhinge(threespanwalk)==0){//if no hooks in hinge
//					printf("no hook. recording\n");
					ConnectInfo[tspans_nrr[i][j]][index] = tspans_nrr[tspans_outsection[i][j]][k];
					num_outtspans[tspans_nrr[i][j]]++;
					index++;
				}
//				else{ printf("hooked. not recording\n");}
			}
			if(index==1){
//				printf("tspan %d cannot connect to any other tspan\n", tspans_nrr[i][j]);
				l=0;
				while(invalidtspans[l]!=0){
					l++;
				}
				invalidtspans[l]=tspans_nrr[i][j];
				
			}
		}
	}

//	printf("The invalid tspans numbers are: ");
	i=0;
	while(invalidtspans[i]!=0){
//		printf("%d ", invalidtspans[i]);
		i++;
	}
	printf("\nThere are %d invalid tspans (absorbing).\n\n", i);

/*
	for(i=1; i<=total2spans; i++){
		printf("2span %d. num_outtspans=%d: ", i, num_outtspans[i]);
		for(j=1; j<=num_outtspans[i]; j++){
			if(ConnectInfo[i][j]!=0){
				printf("%d ", ConnectInfo[i][j]);
			}
		}
		printf("\n");
	}
*/



//delete absorbing states.
	for(i=1; i<=total2spans; i++){
//		printf("2span %d\n", i);
		for(j=1; j<=num_outtspans[i]; j++){//for each outtspan, check if it's in invalid sections.
//			printf("checking the %dth outtspan of %d total. This one is tspan number %d\n", j, num_outtspans[i], ConnectInfo[i][j]);
			for(k=0; k<=total2spans-1; k++){
				if(invalidtspans[k]==ConnectInfo[i][j]){
//					printf("2span %d connects to 2span %d, but it's an invalid 2span. Deleting.\n", i, ConnectInfo[i][j]);
					l=j;
					while(ConnectInfo[i][l+1]!=NULL){
						ConnectInfo[i][l] = ConnectInfo[i][l+1];
						l++;
					}
					num_outtspans[i]--;
					j--;
					break;
				}
			}
		}
	}


	//repeat check for absorbing states until there are none.
	int foundsome=1;
	int recorded;

	while(foundsome==1){
		foundsome=0;
		for(i=1; i<=total2spans; i++){
//			printf("checking tspan %d\n", i);
			if(num_outtspans[i]==0){
//				printf("tspan %d has no outtspans. check if it's already been recorded in invalidtspans[]\n", i);
				//check if it's in invalidtspans[]
				recorded=0;
				j=0;
				while(invalidtspans[j]!=0){
					if(invalidtspans[j]==i){
//						printf("it has already been recorded in spot %d\n", j);
						recorded=1;
					}
					j++;
				}
				if(recorded==0){
//					printf("found a new invalid tspan. must add it to the list. exitting.\n");
					foundsome=1;
					exit(1);
				}
			}
		}
	}


/*
	for(i=1; i<=total2spans; i++){
		printf("2span %d. num_outspans=%d: ", i, num_outtspans[i]);
		for(j=1; j<=num_outtspans[i]; j++){
			printf("%d ", ConnectInfo[i][j]);
		}
		printf("\n");
	}
*/






	printf("Performing check for irreducibility.\n");
	int testsection=1;	//actually the test 2span
	int tstsec=testsection;


	//check for irreducibility. Graph theory way.
	unsigned long int outs[total2spans+1];		//number of out paths from this 2span in the tree (if zero, we know its a leaf)
	unsigned long int ConnSecs[total2spans+1];	
	unsigned long int ConnectingSections[total2spans+1];//will indicate if the 2span has been determined if it can be connected

	printf("The test 2-span for irreducibility is 2span %d.\n", testsection);

	for(i=1; i<=total2spans; i++){
		ConnectingSections[i]=0;
		outs[i]=0;
	}
	ConnectingSections[testsection]=1; //initialized that it's in its own class

	int secs_in_class_one = 1;  //actually 2-spans in class1
	int updated=1;	//updated=1 when there's been an update in the iteration (initialized to 1)

	while(updated==1){
		updated=0;
		for(i=1; i<=total2spans; i++){
//			printf("Looking at 2span %d\n", i);
			if(ConnectingSections[i]==1){
				for(j=1; j<=num_outtspans[i]; j++){
//					printf("2span %d connects to 2span %d\n", i, ConnectInfo[i][j]);
					if(ConnectingSections[ConnectInfo[i][j]]==0){
						outs[i]++;
						updated=1;
						secs_in_class_one++;
						ConnectingSections[ConnectInfo[i][j]]=1;
					}
				}
			}
		}
	}

	printf("2span %d can connect to %d 2spans in class one.\n", tstsec, secs_in_class_one);
/*
	printf("2spans which are 'leaves': ");
	for(i=1; i<=total2spans; i++){
		if(ConnectingSections[i]==1 && outs[i]==0){
			printf("%d ", i);
		}
	}
	printf("\n");

*/

	//check if each 2span that's a leaf can connect to 2span tstsec.

	int can=0;	//indicates if the leaf can connect to tstsec
	int foundinvalid=0;
	for(k=1; k<=total2spans; k++){
		if(ConnectingSections[k]==1 && outs[k]==0){	//k is a tspan that's a leaf.
			for(i=1; i<=total2spans; i++){
				ConnSecs[i]=0;
			}
			can=0;
			ConnSecs[k]=1;
			updated=1;
			while(updated==1){
				updated=0;
				for(i=1; i<=total2spans; i++){
					if(ConnSecs[i]==1){
						for(j=1; j<=num_outtspans[i]; j++){
							if(ConnectInfo[i][j]!=0){
								if(ConnectInfo[i][j]==tstsec){
									can=1;
									break;
								}
								else if(ConnSecs[ConnectInfo[i][j]]==0){
									updated=1;
									ConnSecs[ConnectInfo[i][j]]=1;
								}
							}
						}
					}
				}

			}
			if(can==1){
//				printf("leaf %d can connect to section %d\n", k, tstsec);
			}
			else{
//				printf("leaf %d CANNOT CONNECT to 2span %d\n", k, tstsec);
				//check if tspan k is in invalidtspans
				foundinvalid=0;
				i=0;
				while(invalidtspans[i]!=0){
					if(invalidtspans[i]==k){
//						printf("BUT 2SPAN %d IS INVALID. CONTINUE ON.\n", k);
						foundinvalid=1;
						break;
					}
					i++;
				}
				if(foundinvalid==0){
					printf("Found a leaf that cannot connect back to tspan %d.\n", tstsec);
					printf("Class 1 with tspan %d is not irreducible\n", tstsec);
					//print 2span tstsec and 2span k
					for(i=1;i<=max_sections; i++){
						for(j=1; j<=num_outsections[i]; j++){
							if(tspans_nrr[i][j]==k || tspans_nrr[i][j]==tstsec){
								//print tspan
								printf("Printing tspan %d: ",tspans_nrr[i][j]);
								l=0;
								while(tspans_walk[i][j][l]!=0){
									printf("%d ", tspans_walk[i][j][l]);
									l++;
								}
								printf("\n");
								if(tspans_nrr[i][j]==k){
//									printf("tspan %d can connect to %d other tspans\n", k, num_outtspans[k]);
								}
							}
						}
					}
					printf("Exitting\n");
					exit(1);
				}

			}
		}
	}

/*
	int nonzeroterms=0;
	int zeroterms=0;
	for(i=1; i<=total2spans; i++){
		if(ConnectingSections[i]!=0){
			nonzeroterms++;
		}
		else{
			zeroterms++;
		}
	}
	printf("ConnectingSections[] has %d non-zero entries. Also, %d zero entries.\n", nonzeroterms, zeroterms);
*/


	//find tspans which are not in class one
	int notinclass[total2spans];
	for(i=0; i<=total2spans-1; i++){
		notinclass[i]=0;
	}
	index=0;
	int count=0;
	for(i=1; i<=total2spans; i++){
		if(ConnectingSections[i]==0 && num_outtspans[i]>0){
			foundinvalid=0;
			j=0;
			while(invalidtspans[j]!=0){
				if(invalidtspans[j]==i){
					foundinvalid=1;
					break;
				}
				j++;
			}
			if(foundinvalid==0){			
//				printf("2span %d is not in class one\n", i);
				notinclass[index]=i;
				index++;
				count++;
/*				for(j=1; j<=max_sections; j++){
					for(k=1; k<=num_outsections[j]; k++){
						if(tspans_nrr[j][k]==i){
						printf("Printing tspan %d: ",tspans_nrr[j][k]);
							l=0;
							while(tspans_walk[j][k][l]!=0){
								printf("%d ", tspans_walk[j][k][l]);
								l++;
							}
							printf("\n");
							break;
						}
					}
				}
*/			}
		}
	}

printf("There are %d 2spans not in class one.\n", count);


//finding which tspans can connect to the 2spans not in class one
	int foundavalid=0;
	for(i=1; i<=total2spans; i++){
		for(j=1; j<=num_outtspans[i]; j++){
			k=0;
			while(notinclass[k]!=0){
				if(ConnectInfo[i][j]==notinclass[k]){
					foundavalid=1;
					printf("tspan %d connects to %d.\n", i, notinclass[k]);
				}
				k++;
			}
		}
	}

	if(foundavalid==1){
		printf("We found a tspan which is not in class one and is a valid tspan. There exists more than one class, so we do NOT have irreducibility!");
	}
	else{
		printf("We have irreducibility! Every valid tspan not in class one is invalid (nothing can get into it)\n");
//delete these invalid 2spans.
/*		i=0;
		while(notinclass[i]!=0){
			for(j=1; j<=num_outtspans[notinclass[i]]; j++){
				ConnectInfo[notinclass[i]][j]=0;
			}
			num_outtspans[notinclass[i]]=0;
			i++;
		}
*/
	}




//check for periodicity


	int testvector[total2spans+1];
	for(i=1; i<=total2spans; i++) {
		if(i==testsection){	//initialize testvector be all zero, expect for the entry corresponding to the testsection
			testvector[i]=1;
		}
		else{
			testvector[i]=0;
		}
	}
	int temptestvector[total2spans+1];

	unsigned long int testpower=201;	//power of the matrix that will be multiplied by testvector (steps in markov chain)
	for(k=1; k<=testpower; k++){
		for(i=1; i<=total2spans; i++){
			temptestvector[i] = 0; //set temp=vector of zeros
		}
		for(i=1; i<=total2spans; i++){
			for(j=1; j<=num_outtspans[i]; j++){
				temptestvector[ConnectInfo[i][j]] += testvector[i]; //temp=test*matrix
			}
		}
		for(i=1; i<=total2spans; i++){
			testvector[i]=temptestvector[i];	//set test=temp
		}
//		for(i=1; i<=total2spans; i++){
//			if(testvector[i] != 0){
//				printf("testvector[%d]=%d\n", i, testvector[i]);
//			}
//		}
//		printf("iter=%d\n", k);
	}

	int itsinvalid=0;
	int isperiodic=0;
	for(i=1; i<=total2spans; i++){
//		printf("checking tspan %d\n", i);
		if(testvector[i]==0){ //tspan 'testsection' cannot connect to tspan i in 'testpower' moves
			//check if it's an invalid tspan.
//			printf("tspan %d is zero...\n", i);
			itsinvalid=0;
			j=0;
			while(invalidtspans[j]!=0){
				if(i==invalidtspans[j]){
//					printf("but it's absorbing\n");
					itsinvalid=1;
					break;
				}
				j++;
			}
			if(itsinvalid==0){
				j=0;
				while(notinclass[j]!=0){
					if(i==notinclass[j]){
//						printf("but it's only a starter\n");
						itsinvalid=1;
						break;
					}
					j++;
				}
			}
			//finshed check
			if(itsinvalid==0){
				printf("tspan %d cannot connect to tspan %d in %d moves\n", i, tstsec, testpower);
				printf("PERIODIC!! exitting\n");
				isperiodic=1;
				exit(1);
			}
			
		}
	}
	if(isperiodic==0){
		printf("Not periodic! (testpower=%d)\n", testpower);
	}





#if defined(TS)
	L_Evector[0] = vecalloc(0,max_tspans); /*global*/
	L_Evector[1] = vecalloc(0,max_tspans); /*global*/ /*temp*/
	R_Evector[0] = vecalloc(0,max_tspans); /*global*/
	R_Evector[1] = vecalloc(0,max_tspans); /*global*/ /*temp*/
#else
	L_Evector[0] = vecalloc(0,max_sections); /*global*/
	L_Evector[1] = vecalloc(0,max_sections); /*global*/ /*temp*/
	R_Evector[0] = vecalloc(0,max_sections); /*global*/
	R_Evector[1] = vecalloc(0,max_sections); /*global*/ /*temp*/
#endif

#if defined(connectivity_inv)
	max_eval_LRvec(connectivity_inv); /*rtflsp not necessary once connectivity is known*/
#else
	printf("\nStarting false position method to determine connective constant.\n");
	#if defined(force)
		printf("\n=> %f is fugacity which yields 1 as maximum eigenvalue (after multiplying matrix by e^%f)\n",
					connectivity_inv=rtflsp(&max_eval_LRvec,x_L,x_R, 0.0000001), fval);
	#else
		printf("\n=> %f is fugacity which yields 1 as maximum eigenvalue\n",
					connectivity_inv=rtflsp(&max_eval_LRvec,x_L,x_R, 0.0000001));
	#endif
		printf("\n=> %f is the connective constant.\n",1/connectivity_inv);
#endif

	LmultR = 0;

#if defined(TS)
	for (i=1;i<= max_tspans;i++){
		LmultR += (L_Evector[0][i])*(R_Evector[0][i]);
	}
#else
	for (i=1;i<= max_sections;i++){
		LmultR += (L_Evector[0][i])*(R_Evector[0][i]);
	}
#endif

/*
//CHECK IF LEFT EIGENVECTOR IS A LEFT EIGENVECTOR
	//initialize resultant vectors to zeros
	for(i=1; i<=max_sections; i++){
		L_Evector[1][i] = 0;
		R_Evector[1][i] = 0;
	}
	float pwfugacity;
	//Take left eigenvector, multiply by trans matrix, and see if you get the same number
	for(i=1; i<=max_sections; i++){	//for each section
		for(j=1; j<=num_outsections[i]; j++){	//for each section that i connects to,
			pwfugacity = pow(connectivity_inv, tspans_edges[i][j]);
			L_Evector[1][tspans_outsection[i][j]] += pwfugacity*(L_Evector[0][i])*exp(fval);
			R_Evector[1][i] += pwfugacity*(R_Evector[0][tspans_outsection[i][j]])*exp(fval);
		}
	}

	//check if you have the same numbers
	for(i=1; i<=max_sections; i++){
		if(num_outsections[i]>0){
			if((L_Evector[0][i] - L_Evector[1][i] > 0.0000001) || (L_Evector[0][i] - L_Evector[1][i] < -0.0000001)){
				printf("L_Evector[0][%d] = %f = %e, L_Evector[1][%d] = %f = %e\n", i, L_Evector[0][i], L_Evector[0][i], i, L_Evector[1][i], L_Evector[1][i]);
//				printf("R_Evector[0][%d] = %f = %e, R_Evector[1][%d] = %f = %e\n", i, R_Evector[0][i], R_Evector[0][i], i, R_Evector[1][i], R_Evector[1][i]);
			}
		}
	}
	printf("x_0 = %f\n", connectivity_inv);
*/


	/* LmultR = L_Evector*R_Evector *** use this for normalizing the vectors */
	Ephi_Ephi = (L_Evector[0][1])*(R_Evector[0][1])/LmultR;
	Beta = get_Beta(connectivity_inv);

#if defined(TS)
	Alpha = connectivity_inv*connectivity_inv*exp(2*fval)*Ephi_Ephi/Beta;
#else /* defined(SS) */
	Alpha = connectivity_inv*connectivity_inv*connectivity_inv*connectivity_inv*exp(4*fval)*Ephi_Ephi/Beta;
#endif

	/*may now free L_Evector[1] this has yet to be implemented */
	/*may now free R_Evector[1] this has yet to be implemented */
	printf("LmultR=%f\n", LmultR);
	printf("L_Evector[0][1]=%f\n", L_Evector[0][1]);
	printf("R_Evector[0][1]=%f\n", R_Evector[0][1]);

	newline;
	printf("Alpha = %f = %e\n",Alpha,Alpha);
 	printf("Beta  = %f = %e\n\n",Beta,Beta);
	printf("L_Evector[phi]*R_Evector[phi] = %f = %e\n\n",Ephi_Ephi,Ephi_Ephi);

/************************************************************************************************/	
/*showmenu*/ /*Needs to be implemented*/

	kappa = connectivity_inv;
	/*system("clear");*/

//USE TO CHECK IF EXPECTED NUMBER OF OCCURRENCES ADD UP TO ONE (after not dividing by beta)
//	double sumofprobs = sumtoone(Beta, LmultR, kappa);

//USE TO CHECK IF THE TRANSITION PROBABILITES ADD UP TO ONE FOR EACH SECTION/TSPAN
//	double transprobs = transprobcheck(Beta, LmultR, kappa);

//USE IF YOU WANT TO SEE WHICH TSPANS HAVE THE SAME EXPECTED NUMBER OF OCCURRENCES (GROUPS THEM)
//	groupings(Beta, LmultR, kappa);	//different results depending on the tolerance you set for considering them equal
								//can change this tolerance in the function groupings()
//	group22(Beta, LmultR, kappa);		//prints the number of edges of the 2-spans in group 22
//	group145(Beta, LmultR, kappa);	//prints the number of edges of the 2-spans in group 145
//	group330(Beta, LmultR, kappa);	//prints the number of edges of the 2-spans in group 330

//USE TO FIND THE SECTION WITH THE EXPECTED NUMBER OF OCCURRENCES
//ONLY USE IN SECTION-SECTION MODE
//	findsection(Beta, LmultR, kappa);

//USE IF YOU WANT TO TAKE IN POLYGON FILES AND FIND THE SUM OF THE EXPECTED NUMBER OF OCCURRENCES (ALSO PRINTS PROBABILITY) OF
//THE K-SPAN CREATED BY REMOVING THE ENDHINGES OF THE POLYGON.
//	probkspan(Beta, LmultR, kappa);
//	probkspan2(Beta, LmultR, kappa);

//USE IF YOU WANT TO ENTER PATTERNS AND GET THE EXPECTED NUMBER OF OCCURRENCES OF A K-SPAN
//	Enum_Enum = enterpats(Beta, LmultR, kappa);

#if defined(TS)
//	graphexpoccurvsedges(Beta, LmultR, kappa);		//creates .csv file to graph expected number of occurrences vs number of edges
//	graphexpoccurvsnumconnect(Beta, LmultR, kappa);	//creates .csv file to graph expected number of occurrences vs number of
											//sections that can come after the 2-span
//	graphexpoccurvsconnectto(Beta, LmultR, kappa);	//creates .csv file to graph expected number of occurrences vs number of
											//sections that can come before the 2-span
//	graphexpoccurvstotalconnect(Beta, LmultR, kappa);	//creates .csv file to graph expected number of occurrences vs the number
											//of sections that can come before the 2-span plus the number of sections
											//that can come after the 2-span

//	graphedgesvsnumconnect(Beta, LmultR, kappa);		//graphs edges in a 2-span vs number of sections that can come after the 2-span
//	graphedgesvsconnectto(Beta, LmultR, kappa);		//graphs edges in a 2-span vs number of sections that can come before the 2-span
//	graphedgesvstotalconnect(Beta, LmultR, kappa);	//graphs edges in a 2-span vs number of sections that can come before the 2-span plus
											//the number of sections that can come after the 2-span

#endif

/*
	printf("L_Evector[0][%d]=%f=%e, R_Evector[0][%d]=%f=%e\n", 372, L_Evector[0][372], L_Evector[0][372], 372, R_Evector[0][372], R_Evector[0][372]);
	printf("L_Evector[0][%d]=%f=%e, R_Evector[0][%d]=%f=%e\n", 48, L_Evector[0][48], L_Evector[0][48], 48, R_Evector[0][48], R_Evector[0][48]);
	printf("L_Evector[0][%d]=%f=%e, R_Evector[0][%d]=%f=%e\n", 436, L_Evector[0][436], L_Evector[0][436], 436, R_Evector[0][436], R_Evector[0][436]);
	printf("L_Evector[0][%d]=%f=%e, R_Evector[0][%d]=%f=%e\n", 230, L_Evector[0][230], L_Evector[0][230], 230, R_Evector[0][230], R_Evector[0][230]);

	printf("num_outsections[%d]=%d\n", 372, num_outsections[372]);
	printf("num_outsections[%d]=%d\n", 48, num_outsections[48]);

	int edgesum=0;
	for(i=1; i<=num_outsections[372]; i++){
		edgesum += tspans_edges[372][i];
		printf("tspans_edges[372][%d]=%d\n", i, tspans_edges[372][i]);
	}
	printf("edge sum for 372 is %d\n", edgesum);

	edgesum=0;
	for(i=1; i<=num_outsections[48]; i++){
		edgesum += tspans_edges[48][i];
		printf("tspans_edges[48][%d]=%d\n", i, tspans_edges[48][i]);
	}
	printf("edge sum for 48 is %d\n", edgesum);


	//check how many sections go into section 372 and 48
	int incount372=0;
	int incount48=0;

	for(i=1; i<=max_sections; i++){
		for(j=1; j<=num_outsections[i]; j++){
			if(tspans_outsection[i][j]==372){
//				printsection(i);
				incount372++;
			}
			else if(tspans_outsection[i][j]==48){
				printsection(i);
				incount48++;
			}
		}
	}

	printf("%d sections connect to section 372\n", incount372);
	printf("%d sections connect to section 48\n", incount48);

	printf("num_outsections[74]=%d\n", num_outsections[74]);
	for(i=1; i<=num_outsections[74]; i++){
		printf("outsection of 74 is %d\n", tspans_outsection[74][i]);
	}



	printf("number of two edge sections=%d\n", twoedgesections);
	printf("number of four edge sections=%d\n", fouredgesections);
	printf("number of six edge sections=%d\n", sixedgesections);

	printf("number of two edge sections thrown away=%d\n", twoedgethrown);
	printf("number of four edge sections thrown away=%d\n", fouredgethrown);
	printf("number of six edge sections thrown away=%d\n", sixedgethrown);

	printf("total number of valid sections=%d=%d. They are sections: ", twoedgesections+fouredgesections+sixedgesections, valid_sections);
	i=0;
	while(recordedsections[i]!=0){
		printf("%d\n",recordedsections[i]);
		printsection(recordedsections[i]);
		printf("left eigen: %f, right eigen: %f", L_Evector[0][recordedsections[i]], R_Evector[0][recordedsections[i]]);
		i++;
	}
	printf("\n");
	i=0;
	printf("The thrown away sections are: ");
	while(failedsections[i]!=0){
		printf("%d, ",failedsections[i]);
		i++;
	}
*/

	system("date");  /* prints the date and time */
	printf("\nProgram Complete.\n\n");


}				/* end of main */
/***************************************************************************/


/***** start of function conv_to_array *************************************/
/* called only by main */
void
conv_to_array(void) /* tspans_outsection[section_num][1..num_outsections[section_num]] */
{
	int i;
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

			num_outsections[section_num]++; /*count that number of two-spans with section_num as first section*/

		} /* this while loop counts the number of valid two spans with first section = section_num */

		tspans_outsection[section_num] = unsgnlong_vecalloc(1, num_outsections[section_num]);/*allocate memory*/

		tspans_edges[section_num] = unsgnlong_vecalloc(1, num_outsections[section_num]);/*allocate memory e=edges*/
							/*for recording the number of edges in the two-spans */
							/*with first section = section_num    */

		tspans_nrr[section_num] = unsgnlong_vecalloc(1, num_outsections[section_num]);/*allocate memory n=number*/
							/*for recording the two-spans number of a two-span which has */
							/* first section = section_num    */
		tspans_walk[section_num] = unsgnlong_vecalloc(1, num_outsections[section_num]);
		for(i=1; i<=num_outsections[section_num]; i++){
			tspans_walk[section_num][i] = vecalloc(0,vM*vL*8-1);
		}




	}

	for (section_num = 1; section_num <= max_sections; section_num++) {
		arr_ent = 1;	/* initialize to zero */
		current_hinge_span[section_num] = first_hinge_span[section_num];
		while ((*current_hinge_span[section_num]).nexthinge != NULL) {
			hinge_to_free = current_hinge_span[section_num];
			current_hinge_span[section_num] = (*current_hinge_span[section_num]).nexthinge;
			free(hinge_to_free);
			tspans_outsection[section_num][arr_ent] = (*current_hinge_span[section_num]).outorder;
			tspans_edges[section_num][arr_ent] = (*current_hinge_span[section_num]).edgecount;
			//copy over walk
			for(i=0; i<=vM*vL*8-1; i++){
				tspans_walk[section_num][arr_ent][i] = (*current_hinge_span[section_num]).walk[i];
			}
			arr_ent++;
		}
	}

	for (section_num = 1; section_num <= max_sections; section_num++) {
		qcksrtII(num_outsections[section_num], tspans_outsection[section_num],tspans_edges[section_num], tspans_walk[section_num]);
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


	int             ii;		/*vertex number in vertical direction*/
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

			//check if there are any "hooks" in the 2-span hinge. If no hook, record it.
			if(hookedhinge(curwalk)==0){ // 0 if curwalk[] does not contain a hook
				if(allstraight()==0){	//0 if curwalk does NOT have all straight edges
					valid_2_spans++; /*This is a valid 2 span so it is counted*/
					record_1_span(); /*The component one spans are therefore valid so need to be recorded*/
					recordtemplate(pointordNum); /*the information contained in ordertemplate is recorded*/
				}
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
void
recordtemplate( int (*pointordNum)[3])
/* This function records the pertinent information about a two-span */
/* which has been verified as being valid into a linked list */
{
	int i,j;
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
		i++;
	}

	return;
}

#if defined(TS)
	#if defined(force)
		#include "../transfer_matrix/pw_meth_ts_LRvec_fcheck_irreducible.c" /* for applying the power method to the 2-span 2-span matrix */
		#include "../analysis/statistics/get_Beta_tsCheck.c"      /* for calculating Beta */
	#else
		#include "../transfer_matrix/pw_meth_ts_LRvec.c" /* for applying the power method to the 2-span 2-span matrix */
		#include "../analysis/statistics/get_Beta_ts.c"      /* for calculating Beta */ 
	#endif

#endif

#if defined(SS)
	#if defined(force)
		#include "../transfer_matrix/pw_meth_ss_LRvec_fcheck.c" /* for applying the power method to the section-section matrix */
		#include "../analysis/statistics/get_Beta_ssCheck.c"	/* for calculating Beta */
	#else
		#include "../transfer_matrix/pw_meth_ss_LRvec.c" /* for applying the power method to the section-section matrix */
		#include "../analysis/statistics/get_Beta_ss.c"	/* for calculating Beta */
	#endif
	
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



#endif







/***************************************************************************/
/*
unsigned long int
getval(void)
// gets an integer from the user
{
	char            input[100];
	return (atoi(gets(input)));
}
*/


/*
unsigned long int
getval_mnmx(unsigned long int min, unsigned long int max)
// get unsigned long int from user between the values min and max
{
	char            input[100];
	unsigned long int x;
	gets(input);
	x = atoi(input);
	if ((min <= x) && (x <= max)) {
		return (x);
	} else {
		return (min);	// if the value isn't in the range, return the min
	}
}
*/


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
	int i;
	for(i=0; i<=vM*vL*8-1; i++){
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






/*THE FOLLOWING WAS ADDED MAY 26, 1998 */

/* actually need to return two-span # */
/* or first section # and last section # */
/* and total number of edges */
/*
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
	do	// get first section
	{	
		printf("\nEnter section #1:\n\n");
		scan_template(0);	// enter section #n
		valid = entsecval(0);	// validity test, but not necessary to test connecting
	} while (valid == 0);

	sec_1_num = num_section(0);

	while (n < k) {

		n++;

		side ^= 1; //switch sides

		do	// get number of edges in hinge 
		{
			printf("\nEnter number of edges in hinge:");
			scanf("%i",&(totedges));
			getval_mnmx(1,100);
			//scanf("%i", &(totedges));
			printf("\ntotedges are:%i\n",totedges);
			valid=1;
		} while (valid == 0);

		do	// get next section
		{
			printf("\nEnter section #%i.\n\n", n);
			scan_template(side);	// enter section #n
			valid = entsecval(side) * entseccon(side^1);	 //valid and connecting
		} while (valid == 0);
	}

	sec_k_num = num_section(side);
}
*/
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

int hookedhinge(int hingewalk[]){	//returns 1 if hingewalk's hinge has a hook in it. 0 otherwise
	int i;
//	printf("hookedhinge() is running\n");

/*	printf("hingewalk is:\n");
	i=0;
	while(hingewalk[i]!=0){
		printf("%d\n", hingewalk[i]);
		i++;
	}
*/
	int hooked=0;
	i=0;
	while(hingewalk[i+2]!=0){
		if(hingewalk[i]==5 && hingewalk[i+1]==3 && hingewalk[i+2]==6){	//+z+y-z
			return 1;
		}
		else if(hingewalk[i]==5 && hingewalk[i+1]==4 && hingewalk[i+2]==6){	//+z-y-z
			return 1;
		}
		else if(hingewalk[i]==5 && hingewalk[i+1]==1 && hingewalk[i+2]==6){	//+z+x-z
			return 1;
		}
		else if(hingewalk[i]==5 && hingewalk[i+1]==2 && hingewalk[i+2]==6){	//+z-x-z
			return 1;
		}
	/////
		else if(hingewalk[i]==6 && hingewalk[i+1]==3 && hingewalk[i+2]==5){	//-z+y+z
			return 1;
		}
		else if(hingewalk[i]==6 && hingewalk[i+1]==4 && hingewalk[i+2]==5){	//-z-y+z
			return 1;
		}
		else if(hingewalk[i]==6 && hingewalk[i+1]==1 && hingewalk[i+2]==5){	//-z+x+z
			return 1;
		}
		else if(hingewalk[i]==6 && hingewalk[i+1]==2 && hingewalk[i+2]==5){	//-z-x+z
			return 1;
		}
	/////
		else if(hingewalk[i]==1 && hingewalk[i+1]==3 && hingewalk[i+2]==2){	//+x+y-x
			return 1;
		}
		else if(hingewalk[i]==1 && hingewalk[i+1]==4 && hingewalk[i+2]==2){	//+x-y-x
			return 1;
		}
		else if(hingewalk[i]==1 && hingewalk[i+1]==5 && hingewalk[i+2]==2){	//+x+z-x
			return 1;
		}
		else if(hingewalk[i]==1 && hingewalk[i+1]==6 && hingewalk[i+2]==2){	//+x-z-x
			return 1;
		}
	/////
		else if(hingewalk[i]==2 && hingewalk[i+1]==3 && hingewalk[i+2]==1){	//-x+y+x
			return 1;
		}
		else if(hingewalk[i]==2 && hingewalk[i+1]==4 && hingewalk[i+2]==1){	//-x-y+x
			return 1;
		}
		else if(hingewalk[i]==2 && hingewalk[i+1]==5 && hingewalk[i+2]==1){	//-x+z+x
			return 1;
		}
		else if(hingewalk[i]==2 && hingewalk[i+1]==6 && hingewalk[i+2]==1){	//-x-z+x
			return 1;
		}
	/////
		else if(hingewalk[i]==3 && hingewalk[i+1]==1 && hingewalk[i+2]==4){	//+y+x-y
			return 1;
		}
		else if(hingewalk[i]==3 && hingewalk[i+1]==2 && hingewalk[i+2]==4){	//+y-x-y
			return 1;
		}
		else if(hingewalk[i]==3 && hingewalk[i+1]==5 && hingewalk[i+2]==4){	//+y+z-y
			return 1;
		}
		else if(hingewalk[i]==3 && hingewalk[i+1]==6 && hingewalk[i+2]==4){	//+y-z-y
			return 1;
		}
	/////
		else if(hingewalk[i]==4 && hingewalk[i+1]==1 && hingewalk[i+2]==3){	//-y+x+y
			return 1;
		}
		else if(hingewalk[i]==4 && hingewalk[i+1]==2 && hingewalk[i+2]==3){	//-y-x+y
			return 1;
		}
		else if(hingewalk[i]==4 && hingewalk[i+1]==5 && hingewalk[i+2]==3){	//-y+z+y
			return 1;
		}
		else if(hingewalk[i]==4 && hingewalk[i+1]==6 && hingewalk[i+2]==3){	//-y-z+y
			return 1;
		}
		i++;
	}

//	printf("hooked=%d\n", hooked);
//	exit(1);
	return hooked;

}

int allstraight(){
//	printf("allstraight() is running\n");
	int i=0;

/*	printf("curwalk is:\n");
	i=0;
	while(curwalk[i]!=0){
		printf("%d\n", curwalk[i]);
		i++;
	}
*/

	i=0;
	while(curwalk[i]!=0){
		if(curwalk[i]==5){	//just entered hinge
			if(curwalk[i+1]!=5){	//not all straight
//				printf("not all straight\n");
				return 0;
			}
			i=i+2;	//get past 55 
		}
		else if(curwalk[i]==6){
			if(curwalk[i+1]!=6){	//not all straight
//				printf("not all straight\n");
				return 0;
			}
			i=i+2;	//get past 66
		}
		i++;
	}
	return 1;

}

#include "../sections/findsection.c"
#include "../sections/printsection.c"

#if defined(TS)
#include "../analysis/statistics/sumtoone_ts.c"
//#include "../sections/enterpatstspan_ts.c"
#include "../sections/enterpats_ts.c"
#include "../analysis/statistics/transprobcheck_ts.c"
#include "../analysis/statistics/probkspan_ts.c"
#else
#include "../analysis/statistics/sumtoone_ss.c"
//#include "../sections/enterpatstspan_ss.c"
#include "../sections/enterpats_ss.c"
#include "../analysis/statistics/transprobcheck_ss.c"
#include "../analysis/statistics/probkspan_ss.c"
#include "../analysis/statistics/probkspan_ss2.c"
#endif

#if defined(TS)
#include "../analysis/counters/tspan_groupings_percentequal.c"
#include "../analysis/statistics/expoccur.c"
#include "../analysis/counters/group22.c"
#include "../analysis/counters/group145.c"
#include "../analysis/counters/group330.c"
#include "../analysis/statistics/graphexpoccurvsedges.c"
#include "../analysis/statistics/graphexpoccurvsnumconnect.c"
#include "../analysis/statistics/graphexpoccurvsconnectto.c"
#include "../analysis/statistics/graphexpoccurvstotalconnect.c"
#include "../analysis/statistics/graphedgesvsnumconnect.c"
#include "../analysis/statistics/graphedgesvsconnectto.c"
#include "../analysis/statistics/graphedgesvstotalconnect.c"
#endif

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

#include "../utils/qcksrtII_int_irreducible.c"	/*void qcksrt(unsigned long int n,unsigned long int arr[],double brr[])*/
				/* sorts array arr[] and makes corresponding changes in array brr[] */

#include "../utils/addtowalk.c"


	









