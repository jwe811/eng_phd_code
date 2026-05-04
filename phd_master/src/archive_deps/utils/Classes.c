/*
This program will generate all 2-spans with a full hinge. Then it will find how many "connecting classes" there are (and how many sections/2spans are in each).
Then it will run the power method on the transfer matrices for each class.


COMMENTS BY JEREMY ENG, MAY 21, 2013

This is the heavily editted and most current, up to date version of Allan Duffy's HMdoc.c.

!!!!!!!!! This version is the "full hinge" version, which means that only 2-spans that have all of its hinge vertices
are generated. This means there are fewer valid 2-spans. Also, an additional check is needed to make sure the compact 2-spans can actually occur in a compact polygon.


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
#define	L 4				/* number of horizontal edges                              */
#define	M 0				/* number of vertical edges                                */

#define force 0				/* fixed force. */

//#define	connectivity_inv 0.388795	/* When the connectivity is known its inverse is defined here.   */
						/* Doing this will save time, as the rtflsp will not need  */
						/* to be executed to determine the connectivity.		 */
						/* The connectivities (inversed) for various cases follow:       */ 
						/* (1-1) 0.547397 , (1-2) 0.437382 (+1) 0.498950, (1-3) 0.388795 */
						/* 0-? values in Klein article */

//#undef	connectivity_inv	/* Use if connectivity is not known, comment out otherwise.	 */
#define x_L 0.1		/* left end point to use in bracketing method to find x_0        */
#define x_R 0.9			/* right end point to use in bracketing method to find x_0       */

#define		SS	/* TS for two-span // SS for section-section or when conective constat is not known!*/
			/* Note: TS version will find connective constant as well, but is slower */
			/* Also need to run TS version to get data on expected number of times */
			/* a given k-span occurs in an n step polygon 			*/

//Only have to set max_tspans if M==0, otherwise it'll be set automatically
//#define max_tspans   963096		/* vector dimension for two-span matrix  		*/
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

//Automatically set max_tspans (counts based on full version)
#if M==1 && L==1
#define max_tspans	64

#elif M==1 && L==2
#define max_tspans	5488

#elif M==1 && L==3
#define max_tspans	510672

#elif M==1 && L==4
#define max_tspans	99413546		//NON-FULL HINGE COUNT=129413546

#elif M==2 && L==2
#define max_tspans	10136280		//NON-FULL HINGE COUNT=10136280
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
	unsigned long int 	edgecount	:6;	/* the number of edges in the hinge		 */
	struct hinge_span 	*nexthinge;		/* pointer to another hinge_span data structure  */

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
#if defined(force)
double			fval=force;			/* fixed force */
#endif

/* note: The following are also global: onespans[], first_hinge_span, current_hinge_span */
/* They are discussed in the structures section above*/

//global arrays that need to be defined here for some weird reason...
unsigned long int ConnectingSections[max_sections+1];
unsigned long int Class_one[max_sections];		//will contain sections in "class one"
unsigned long int Class_two[max_sections];		//will contain sections in "class two"
unsigned long int Class_three[max_sections];	//will contain sections in "class three"
unsigned long int Class_four[max_sections];	//will contain sections in "class four"
unsigned long int Class_five[max_sections];	//will contain sections in "class five"



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

unsigned short int	entsecval(unsigned short int side);
unsigned short int	entseccon(unsigned short int side);

/**** Added Functions  ********************************************************************/
double			sumtoone(double,double,double);
double			transprobcheck(double, double, double);

#if defined(TS)
void				groupings(double, double, double);
double			expoccur(double, double, double, int);
#endif

void				findsection(double, double, double);
void				printsection(int);
int				edgesinsection(int);

unsigned short int	checkconnection(int, unsigned long int*, unsigned long int, unsigned long int, unsigned short int);
void				probkspan(double, double, double);

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

	fflush(stdout);

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
			printf("enterhinge at i=%d, j=%d\n", i, j);
			enterhinge(makeheaderij, makeheaderij, i, j, side, &ordNum);
		}
	}

/*************************************************************************************************/
	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	newline; /* results are printed to the terminal */
	printf("There are a total of %9i potential two-span patterns.\n", --num_2_spans);
	printf("There are a total of %9i   'compact' two-span patterns.\n", valid_2_spans);
	for (byteloc = 0; byteloc < maxOrdNum; byteloc++) {
		for (bitloc = 0; bitloc < 8; bitloc++) {
			if (one_span_byte_bit(byteloc,bitloc)==1) {valid_sections++;}
		}
	}
	printf("There are a total of %9i   valid   section orderings (before removing).\n\n", valid_sections);
/************************************************************************************************/

	printf("converting linked list to arrays...\n");
	conv_to_array();
	printf("done converting to arrays.\n\n");



//FIND CLASSES OF SECTIONS

	unsigned long int valsec[valid_sections];	//array that will hold the valid section numbers (section numbers that actually exist)
	unsigned long int cursec=0;				//used in loop for filling valsec[]

	//initialize valsec to all zeros
	for(i=0; i<=valid_sections-1; i++){
		valsec[i]=0;
	}
	int numinvalid=0;

	for(i=1; i<=max_sections; i++){		//for all possible section numbers
		if(num_outsections[i]>0){	//if this section number exists
			valsec[cursec]=i;	//record it in valsec
			cursec++;		//increase index of valsec 
		}
		else{
//			printf("%d,", i);
			numinvalid++;
		}
	}
//	printf("are the not-valid section numbers.\n");
//	printf("max_sections=%d\n", max_sections);
//	printf("There are %d non-valid sections\n", numinvalid);
	printf("There are %d valid sections\n", cursec);

//	unsigned long int Class_one[valid_sections];		//will contain sections in "class one"
//	unsigned long int Class_two[valid_sections];		//will contain sections in "class two"
//	unsigned long int Class_three[valid_sections];	//will contain sections in "class three"
//	unsigned long int Class_four[valid_sections];	//will contain sections in "class four"
//	unsigned long int Class_five[valid_sections];	//will contain sections in "class five"

	//initialize these vectors to zeros
	for(i=0; i<=valid_sections-1; i++){
		Class_one[i] = 0;
		Class_two[i] = 0;
		Class_three[i] = 0;
		Class_four[i] = 0;
		Class_five[i] = 0;
	}


	unsigned long int testsection=1;	//A section in "class one"
	printf("\nLet Class_one contain section=%d\n\n", testsection);

//Graph Theory way:
//	unsigned long int ConnectingSections[max_sections+1];	//will indicate if the section (index) has been determined if it can be connected to from "testsection"
	for(i=1; i<=max_sections; i++){
		ConnectingSections[i]=0;
	}
	ConnectingSections[testsection]=1; //initialized that it's in its own class
	int secs_in_class_one = 1;

	int updated=1;	//updated=1 when there's been an update in the iteration (initialized to 1)
	while(updated==1){
		updated=0;
		for(i=1; i<=max_sections; i++){
			if(ConnectingSections[i]==1){
				for(j=1; j<=num_outsections[i]; j++){
					if(ConnectingSections[tspans_outsection[i][j]]==0){
						updated=1;
						secs_in_class_one++;
						ConnectingSections[tspans_outsection[i][j]]=1;
					}
				}
			}
		}
	}

	printf("There are %d sections in class one.\n", secs_in_class_one);

//Markov Chain way:
/*
	double testvector[max_sections+1];	//this will be the vector corresponding to the testsection.
	for(i=1; i<=max_sections; i++) {
		if(i==testsection){	//initialize testvector be all zero, expect for the entry corresponding to the testsection
			testvector[i]=1;
		}
		else{
			testvector[i]=0;
		}
	}

	double temptestvector[max_sections+1];	//used for the vector-matrix multiplication

	unsigned long int testpower=200;	//power of the matrix that will be multiplied by testvector (steps in markov chain)
	printf("Using a testpower of %d.\n\n", testpower);
	for(k=1; k<=testpower; k++){
		for(i=1; i<=max_sections; i++){
			temptestvector[i] = 0; //set temp=vector of zeros
		}
		for(i=1; i<=max_sections; i++){
			for(j=1; j<=num_outsections[i]; j++){
				temptestvector[tspans_outsection[i][j]] += testvector[i]*tspans_edges[i][j]; //temp=test*matrix
			}
		}
		for(i=1; i<=max_sections; i++){
			temptestvector[i] = temptestvector[i] + testvector[i];  //incase of periodicity, once an entry is non-zero, keep it as non-zero
			testvector[i]=temptestvector[i];	//set test=temp
		}
//		for(i=1; i<=max_sections; i++){
//			if(testvector[i] != 0){
//				printf("testvector[%d]=%d\n", i, testvector[i]);
//			}
//		}
//		printf("iter=%d\n", k);
	}

	unsigned long int num_class_one=0;	//this will count the number of sections in class one
	unsigned long int num_other=0;	//this will count the number of sections NOT in class one

	int indx=0;	//index variable for Class_one[]

	i=0;
	while(valsec[i]>0){	//loop through all valid sections in valsec
//		printf("i=%d, testing section %d\n", i, valsec[i]);
		if(testvector[valsec[i]] == 0){	// check if testvector[] is 0 (if it's zero, then this section doesnt connect to testsection)
//			printf("Section %d is not in class one\n", valsec[i]);
			num_other++;
		}
		else{	//section valsec[i] is in class one
//			printf("Section %d is in class one\n", valsec[i]);
		     num_class_one++;
			Class_one[indx] = valsec[i];	//record that valsec[i] is in Class one
			indx++;
			//remove valsec[i] from valsec[]
			for(j=i; j<=valid_sections-2; j++){
				valsec[j] = valsec[j+1];
			}
			valsec[j]=0;
			i--;
		}
		i++;
	}

	i=0;
	printf("Class one contains %d=%d sections: ", num_class_one, indx);
	while(Class_one[i]>0){
		printf("%d ", Class_one[i]);
		i++;
	}
	printf("\nPrinted %d sections.\n", i);
	printf("\n");

	printf("There are %d valid sections that aren't in class one: ", num_other);
	i=0;
	while(valsec[i]>0){
		printf("%d ", valsec[i]);
		i++;
	}
	printf("\nPrinted %d sections.\n", i);
	printf("\n");


/*

// FINDING CLASS TWO!!!!!!!!!!!!!!!!!!!!!
	testsection=75;
	printf("\nLet Class_two contain section=%d\n\n", testsection);

	for(i=1; i<=max_sections; i++) {
		if(i==testsection){	//initialize testvector be all zero, expect for the entry corresponding to the testsection
			testvector[i]=1;
		}
		else{
			testvector[i]=0;
		}
	}

	for(k=1; k<=testpower; k++){
		for(i=1; i<=max_sections; i++){
			temptestvector[i] = 0; //set temp=vector of zeros
		}
		for(i=1; i<=max_sections; i++){
			for(j=1; j<=num_outsections[i]; j++){
				temptestvector[tspans_outsection[i][j]] += testvector[i]*tspans_edges[i][j]; //temp=test*matrix
			}
		}
		for(i=1; i<=max_sections; i++){
			testvector[i]=temptestvector[i];	//set test=temp
		}
//		for(i=1; i<=max_sections; i++){
//			if(testvector[i] != 0){
//				printf("testvector[%d]=%d\n", i, testvector[i]);
//			}
//		}
//		printf("iter=%d\n", k);
	}

	unsigned long int num_class_two=0;	//this will count the number of sections in class two
	num_other=0;
	indx=0;

	i=0;
	while(valsec[i]>0){	//loop through all valid sections in valsec
//		printf("i=%d, testing section %d\n", i, valsec[i]);
		if(testvector[valsec[i]] == 0){	// check if testvector[] is 0 (if it's zero, then this section doesnt connect to testsection)
//			printf("Section %d is not in class two\n", valsec[i]);
			num_other++;
		}
		else{	//section valsec[i] is in class two
//			printf("Section %d is in class two\n", valsec[i]);
		     num_class_two++;
			Class_two[indx] = valsec[i];	//record that valsec[i] is in Class two
			indx++;
			//remove valsec[i] from valsec[]
			for(j=i; j<=valid_sections-2; j++){
				valsec[j] = valsec[j+1];
			}
			valsec[j]=0;
			i--;
		}
		i++;
	}

	i=0;
	printf("Class two contains %d=%d sections: ", num_class_two, indx);
	while(Class_two[i]>0){
		printf("%d ", Class_two[i]);
		i++;
	}
	printf("\nPrinted %d sections.\n", i);
	printf("\n");

	printf("There are %d valid sections that aren't in class two: ", num_other);
	i=0;
	while(valsec[i]>0){
		printf("%d ", valsec[i]);
		i++;
	}
	printf("\nPrinted %d sections.\n", i);
	printf("\n");


// FINDING CLASS THREE!!!!!!!!!!!!!!!!!!!!!
	testsection=17;
	printf("\nLet Class_three contain section=%d\n\n", testsection);

	for(i=1; i<=max_sections; i++) {
		if(i==testsection){	//initialize testvector be all zero, expect for the entry corresponding to the testsection
			testvector[i]=1;
		}
		else{
			testvector[i]=0;
		}
	}

	for(k=1; k<=testpower; k++){
		for(i=1; i<=max_sections; i++){
			temptestvector[i] = 0; //set temp=vector of zeros
		}
		for(i=1; i<=max_sections; i++){
			for(j=1; j<=num_outsections[i]; j++){
				temptestvector[tspans_outsection[i][j]] += testvector[i]*tspans_edges[i][j]; //temp=test*matrix
			}
		}
		for(i=1; i<=max_sections; i++){
			testvector[i]=temptestvector[i];	//set test=temp
		}
//		for(i=1; i<=max_sections; i++){
//			if(testvector[i] != 0){
//				printf("testvector[%d]=%d\n", i, testvector[i]);
//			}
//		}
//		printf("iter=%d\n", k);
	}

	unsigned long int num_class_three=0;	//this will count the number of sections in class three
	num_other=0;
	indx=0;

	i=0;
	while(valsec[i]>0){	//loop through all valid sections in valsec
//		printf("i=%d, testing section %d\n", i, valsec[i]);
		if(testvector[valsec[i]] == 0){	// check if testvector[] is 0 (if it's zero, then this section doesnt connect to testsection)
//			printf("Section %d is not in class three\n", valsec[i]);
			num_other++;
		}
		else{	//section valsec[i] is in class three
//			printf("Section %d is in class three\n", valsec[i]);
		     num_class_three++;
			Class_three[indx] = valsec[i];	//record that valsec[i] is in Class three
			indx++;
			//remove valsec[i] from valsec[]
			for(j=i; j<=valid_sections-2; j++){
				valsec[j] = valsec[j+1];
			}
			valsec[j]=0;
			i--;
		}
		i++;
	}

	i=0;
	printf("Class three contains %d=%d sections: ", num_class_three, indx);
	while(Class_three[i]>0){
		printf("%d ", Class_three[i]);
		i++;
	}
	printf("\nPrinted %d sections.\n", i);
	printf("\n");

	printf("There are %d valid sections that aren't in class three: ", num_other);
	i=0;
	while(valsec[i]>0){
		printf("%d ", valsec[i]);
		i++;
	}
	printf("\nPrinted %d sections.\n", i);
	printf("\n");


// FINDING CLASS FOUR!!!!!!!!!!!!!!!!!!!!!
	testsection=1816;
	printf("\nLet Class_four contain section=%d\n\n", testsection);

	for(i=1; i<=max_sections; i++) {
		if(i==testsection){	//initialize testvector be all zero, expect for the entry corresponding to the testsection
			testvector[i]=1;
		}
		else{
			testvector[i]=0;
		}
	}

	for(k=1; k<=testpower; k++){
		for(i=1; i<=max_sections; i++){
			temptestvector[i] = 0; //set temp=vector of zeros
		}
		for(i=1; i<=max_sections; i++){
			for(j=1; j<=num_outsections[i]; j++){
				temptestvector[tspans_outsection[i][j]] += testvector[i]*tspans_edges[i][j]; //temp=test*matrix
			}
		}
		for(i=1; i<=max_sections; i++){
			testvector[i]=temptestvector[i];	//set test=temp
		}
//		for(i=1; i<=max_sections; i++){
//			if(testvector[i] != 0){
//				printf("testvector[%d]=%d\n", i, testvector[i]);
//			}
//		}
//		printf("iter=%d\n", k);
	}

	unsigned long int num_class_four=0;	//this will count the number of sections in class four
	num_other=0;
	indx=0;

	i=0;
	while(valsec[i]>0){	//loop through all valid sections in valsec
//		printf("i=%d, testing section %d\n", i, valsec[i]);
		if(testvector[valsec[i]] == 0){	// check if testvector[] is 0 (if it's zero, then this section doesnt connect to testsection)
//			printf("Section %d is not in class four\n", valsec[i]);
			num_other++;
		}
		else{	//section valsec[i] is in class four
//			printf("Section %d is in class four\n", valsec[i]);
		     num_class_four++;
			Class_four[indx] = valsec[i];	//record that valsec[i] is in Class four
			indx++;
			//remove valsec[i] from valsec[]
			for(j=i; j<=valid_sections-2; j++){
				valsec[j] = valsec[j+1];
			}
			valsec[j]=0;
			i--;
		}
		i++;
	}

	i=0;
	printf("Class four contains %d=%d sections: ", num_class_four, indx);
	while(Class_four[i]>0){
		printf("%d ", Class_four[i]);
		i++;
	}
	printf("\nPrinted %d sections.\n", i);
	printf("\n");

	printf("There are %d valid sections that aren't in class four: ", num_other);
	i=0;
	while(valsec[i]>0){
		printf("%d ", valsec[i]);
		i++;
	}
	printf("\nPrinted %d sections.\n", i);
	printf("\n");


// FINDING CLASS FIVE!!!!!!!!!!!!!!!!!!!!!
	testsection=10482;
	printf("\nLet Class_five contain section=%d\n\n", testsection);

	for(i=1; i<=max_sections; i++) {
		if(i==testsection){	//initialize testvector be all zero, expect for the entry corresponding to the testsection
			testvector[i]=1;
		}
		else{
			testvector[i]=0;
		}
	}

	for(k=1; k<=testpower; k++){
		for(i=1; i<=max_sections; i++){
			temptestvector[i] = 0; //set temp=vector of zeros
		}
		for(i=1; i<=max_sections; i++){
			for(j=1; j<=num_outsections[i]; j++){
				temptestvector[tspans_outsection[i][j]] += testvector[i]*tspans_edges[i][j]; //temp=test*matrix
			}
		}
		for(i=1; i<=max_sections; i++){
			testvector[i]=temptestvector[i];	//set test=temp
		}
//		for(i=1; i<=max_sections; i++){
//			if(testvector[i] != 0){
//				printf("testvector[%d]=%d\n", i, testvector[i]);
//			}
//		}
//		printf("iter=%d\n", k);
	}

	unsigned long int num_class_five=0;	//this will count the number of sections in class five
	num_other=0;
	indx=0;

	i=0;
	while(valsec[i]>0){	//loop through all valid sections in valsec
//		printf("i=%d, testing section %d\n", i, valsec[i]);
		if(testvector[valsec[i]] == 0){	// check if testvector[] is 0 (if it's zero, then this section doesnt connect to testsection)
//			printf("Section %d is not in class five\n", valsec[i]);
			num_other++;
		}
		else{	//section valsec[i] is in class five
//			printf("Section %d is in class five\n", valsec[i]);
		     num_class_five++;
			Class_five[indx] = valsec[i];	//record that valsec[i] is in Class five
			indx++;
			//remove valsec[i] from valsec[]
			for(j=i; j<=valid_sections-2; j++){
				valsec[j] = valsec[j+1];
			}
			valsec[j]=0;
			i--;
		}
		i++;
	}

	i=0;
	printf("Class five contains %d=%d sections: ", num_class_five, indx);
	while(Class_five[i]>0){
		printf("%d ", Class_five[i]);
		i++;
	}
	printf("\nPrinted %d sections.\n", i);
	printf("\n");

	printf("There are %d valid sections that aren't in class five: ", num_other);
	i=0;
	while(valsec[i]>0){
		printf("%d ", valsec[i]);
		i++;
	}
	printf("\nPrinted %d sections.\n", i);
	printf("\n");
*/




//NOW RUN POWER METHOD FOR A CLASS O_O

//Class one:
	//take sections not in class one "out" of the matrix
	//need to make the appropriate changes to num_outsections, tspans_edges, tspans_outsection.

	unsigned long int m; //looping variable
	unsigned long int dummy;
	int num_remove = valid_sections - num_class_one;

	printf("Using transfer matrix for class one (class which contains section %d).\n", testsection);
	printf("There are %d sections in this class--must remove %d sections that aren't in this class.\n", num_class_one, num_remove);


	for(k=0; k<=num_remove-1; k++){	//loop through all sections that don't connect
//		printf("Taking out section %d\n", valsec[k]);
//		printsection(valsec[k]);
		for(j=1; j<=num_outsections[valsec[k]]; j++){
			tspans_edges[valsec[k]][j] = 0;						//CHANGED EDGES
			tspans_outsection[valsec[k]][j] = 0;					//CHANGED OUTSECTION
		}
		num_outsections[valsec[k]]=0;                            //CHANGED NUM_OUT
	}

	//num_outsections, tspans_edges, tspans_outsection are now adjusted after taking out sections that aren't connected to the test section.

	printf("Finished taking out invalid sections\n");












BLAH


#if defined(TS)
	L_Evector[0] = vecalloc(0,max_tspans); //global
	L_Evector[1] = vecalloc(0,max_tspans); //global temp
	R_Evector[0] = vecalloc(0,max_tspans); //global
	R_Evector[1] = vecalloc(0,max_tspans); //global temp
#else
	L_Evector[0] = vecalloc(0,max_sections); //global
	L_Evector[1] = vecalloc(0,max_sections); //global temp
	R_Evector[0] = vecalloc(0,max_sections); //global
	R_Evector[1] = vecalloc(0,max_sections); //global temp
#endif

#if defined(connectivity_inv)
	max_eval_LRvec(connectivity_inv); //rtflsp not necessary once connectivity is known
#else
	printf("Starting false position method to determine connective constant. (Using T.M. of Class one)\n");
	#if defined(force)
		printf("\n=> %f is fugacity which yields 1 as maximum eigenvalue (after multiplying matrix by e^%f) (using the T.M. for class one containing section=%d)\n",
					connectivity_inv=rtflsp(&max_eval_LRvec,x_L,x_R, 0.0000001), fval, testsection);
	#else
		printf("\n=> %f is fugacity which yields 1 as maximum eigenvalue\n",
					connectivity_inv=rtflsp(&max_eval_LRvec,x_L,x_R, 0.0000001));
	#endif
		printf("\n=> %f is the connective constant.\n",1/connectivity_inv);
#endif

	printf("REMINDER THAT THIS IS FOR CLASS ONE! (Contains section %d)\n", testsection);
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

	int NonZeroSec=0;
	int NonZeroTspan=0;
	LmultR = 0;
#if defined(TS)
	for (i=1;i<= max_tspans;i++){
		LmultR += (L_Evector[0][i])*(R_Evector[0][i]);
		if(L_Evector[0][i]>0){
		  NonZeroTspan++;
		}
	}
	printf("I BELIEVE THAT THERE ARE %d VALID COMPACT T-SPANS\n", NonZeroTspan);
#else
	for (i=1;i<= max_sections;i++){
		LmultR += (L_Evector[0][i])*(R_Evector[0][i]);
		if(L_Evector[0][i]>0){
		  NonZeroSec++;
		  NonZeroTspan+=num_outsections[i];
		}
	}
	printf("I BELIEVE THAT THERE ARE %d VALID FULL-HINGE SECTIONS, AND %d VALID FULL-HINGE T-SPANS (in class one).\n", NonZeroSec, NonZeroTspan);
#endif





	// LmultR = L_Evector*R_Evector *** use this for normalizing the vectors
	Ephi_Ephi = (L_Evector[0][1])*(R_Evector[0][1])/LmultR;
	Beta = get_Beta(connectivity_inv);

#if defined(TS)
	Alpha = connectivity_inv*connectivity_inv*exp(2*fval)*Ephi_Ephi/Beta;
#else // defined(SS)
	Alpha = connectivity_inv*connectivity_inv*connectivity_inv*connectivity_inv*exp(4*fval)*Ephi_Ephi/Beta;
#endif

	//may now free L_Evector[1] this has yet to be implemented 
	//may now free R_Evector[1] this has yet to be implemented 

	newline;
	printf("Alpha = %f = %e\n",Alpha,Alpha);
 	printf("Beta  = %f = %e\n\n",Beta,Beta);
	printf("L_Evector[phi]*R_Evector[phi] = %f = %e\n\n",Ephi_Ephi,Ephi_Ephi);

	kappa = connectivity_inv;

//USE TO CHECK IF EXPECTED NUMBER OF OCCURRENCES ADD UP TO ONE (after not dividing by beta)
//	double sumofprobs = sumtoone(Beta, LmultR, kappa);

//USE TO CHECK IF THE TRANSITION PROBABILITES ADD UP TO ONE FOR EACH SECTION/TSPAN
//	double transprobs = transprobcheck(Beta, LmultR, kappa);

//USE IF YOU WANT TO SEE WHICH TSPANS HAVE THE SAME EXPECTED NUMBER OF OCCURRENCES (GROUPS THEM)
//	groupings(Beta, LmultR, kappa);


//USE TO FIND THE SECTION WITH THE EXPECTED NUMBER OF OCCURRENCES
//ONLY USE IN SECTION-SECTION MODE
//	findsection(Beta, LmultR, kappa);

//USE IF YOU WANT TO TAKE IN A POLYGON FILE AND FIND THE SUM OF THE EXPECTED NUMBER OF OCCURRENCES (ALSO PRINTS PROBABILITY) OF
//THE K-SPAN CREATED BY REMOVING THE ENDHINGES OF THE POLYGON.
//	probkspan(Beta, LmultR, kappa);

//USE IF YOU WANT TO ENTER PATTERNS AND GET THE EXPECTED NUMBER OF OCCURRENCES OF A K-SPAN
//	Enum_Enum = enterpats(Beta, LmultR, kappa);



















/************************************************************************************************/	
/*showmenu*/ /*Needs to be implemented*/

	kappa = connectivity_inv;
	/*system("clear");*/
//	Enum_Enum = enterpats(Beta, LmultR, kappa);
//	printf("Enum_Enum=%f", Enum_Enum);
/*
	printf("number of two edge sections=%d\n", twoedgesections);
	printf("number of four edge sections=%d\n", fouredgesections);
	printf("number of six edge sections=%d\n", sixedgesections);

	printf("number of two edge sections thrown away=%d\n", twoedgethrown);
	printf("number of four edge sections thrown away=%d\n", fouredgethrown);
	printf("number of six edge sections thrown away=%d\n", sixedgethrown);

	printf("total number of valid sections=%d=%d\n", twoedgesections+fouredgesections+sixedgesections, sectioncount);
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
	int		ii;		/*vertex number in vertical direction*/
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
						break;
					}
				}
			}

			if(fullhinge==1){
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

	return;
}

#if defined(TS)
	#include "../transfer_matrix/pw_meth_ts_LRvec_fcheck.c" /* for applying the power method to the 2-span 2-span matrix */
	#include "../analysis/statistics/get_Beta_ts.c"      /* for calculating Beta */
#endif

#if defined(SS)
	#if L==2 && M==2
		#include "../transfer_matrix/pw_meth_ss_LRvec_fcheck_2x2Full.c" /* for applying the power method to the section-section matrix */
	#else
		#include "../transfer_matrix/pw_meth_ss_LRvec_fcheck.c" /* for applying the power method to the section-section matrix */
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

#include "noncrossing.c"	/* include if M = 0 */
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

unsigned short int checkconnection(int section, unsigned long int consection[], unsigned long int curstep, unsigned long int maxsteps, unsigned short int validity){
	if(validity==1){
		return validity;
	}
	int i, j;
	//check if section "section" connects to anything in consection[]
	for(i=1; i<=num_outsections[section]; i++){	//for each outsection of section
		j=0;
		while(consection[j]!=0){
			if(tspans_outsection[section][i]==consection[j]){	//if it connects
				validity=1;
				return validity;
			}
			j++;
		}
	}
	//if reached here, section doesnt connect to consection. Now check recursively if any of section's outsections can connect to consection.
	if(curstep<maxsteps){
		for(i=1; i<=num_outsections[section]; i++){	//for each outsection of section
			validity = checkconnection(tspans_outsection[section][i], consection, curstep+1, maxsteps, validity);
			if(validity==1){
				return validity;
			}
		}
	}
	return validity;
}


#include "../sections/findsection.c"
#include "../sections/printsection.c"
#include "../sections/edgesinsection.c"

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
#endif

#if defined(TS)
#include "../analysis/counters/tspan_groupings_percentequal.c"
#include "../analysis/statistics/expoccur.c"
#endif

#include "unsgnlong_vecalloc.c" /*unsigned long int unsgnlong_vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with unsigned long integer enteries */

#include "vecalloc.c"		/*vec_ent * vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with enteries of type vec_ent */

#include "matalloc.c"		/*mat_ent **matalloc(int rowlow, int rowhigh, int collow, int colhigh)*/
/* used for dynamically allocating memory for matrices with enteries of type mat_ent */


#include "../analysis/statistics/rtflsp.c"	/*double rtflsp(double(*func)(double),double x1,double x2,double xacc)*/
			/*From Numerical Recipes in C section 9.2*/
			/*This uses the regula falsi, (false position method) */



#undef M	//qcksrt uses a different M

#include "qcksrt_int.c"	/*void qcksrt(unsigned long int n,unsigned long int arr[])*/
			/*From Numerical Recipes in C section 8.2*/
			/*Sorts an array arr[1..n] into ascending numerical order using */
			/*the Quicksort algorithm. n is input; arr is replaced on output */
			/*by its sorted rearrangement.*/

#include "qcksrtII_int.c"	/*void qcksrt(unsigned long int n,unsigned long int arr[],double brr[])*/
				/* sorts array arr[] and makes corresponding changes in array brr[] */






	









