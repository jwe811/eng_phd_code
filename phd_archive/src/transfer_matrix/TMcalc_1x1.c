//This nordering version will group sections that have the same left connectivity together, to create a smaller TM
//Only using SS version now.


/*
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
5.	Implemented automatic setting of some defined variables.  max_sections, vec_length, and CS all
	set themselves automatically, depending on the prism size.

Things that you can change (for different system):
M - number of vertical edges in prism
L - number of horizontal edges in prism
force - force in the model

Things you must change manually:
x_L and x_R - used in the root finding. x_0 must lie in between x_L and x_R.
In function main(), can uncomment or comment out the functions:
	sumtoone()
	transprobcheck()
	enterpats()

END OF COMMENTS BY JEREMY ENG
*/

//can now use function avgtwostrings.c to find avg # two string sections per span.


/***************************************************************************/
/* To run in background type:	nohup nice +19 HingeMakerC >& outfile &		*/
/* landscape mode with 1 column type: enscript -1 -r Pmath220 filename.c   */
/* to compile with optimization type: gcc HM97.c -lm -O3 -o nameprogram	*/
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
#define	L 1				/* number of horizontal edges                              */
#define	M 1				/* number of vertical edges                                */

#define force 0			/* fixed force in model. */

//#define	thex 0.361863		/* When the connectivity is known its inverse is defined here (x_0).   */
						/* Doing this will save time, as the rtflsp will not need  */
						/* to be executed to determine the connectivity.		 */
						/* The connectivities (inversed) for various cases follow: (f=0)       */ 
						/* (1-1) 0.547397 , (1-2) 0.437382 , (1-3) 0.388795 , (1-4) 0.361863 , (2-2) 0.366126 */
						/* 0-? values in Klein article */


#define x_L 0.437381	/* left end point to use in bracketing method to find x_0        */
#define x_R 0.437383	/* right end point to use in bracketing method to find x_0       */

#define		SS			//Only using SS version now


/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/

//Automatically set max_sections, based on previous runs.
#if M==1 && L==1
#define max_sections 8
#define max_tspans 48

#elif M==1 && L==2
#define max_sections 73
#define max_tspans 1829

#elif M==1 && L==3
#define max_sections 742
#define max_tspans 70306

#elif M==1 && L==4
#define max_sections 9309
#define max_tspans 3165653

#elif M==1 && L==5
#define max_sections 138038
#define max_tspans 165637127

#elif M==2 && L==2
#define max_sections 2619
#define max_tspans 513585

#elif M==2 && L==3
#define max_sections 138322
#define max_tspans 201423784
#endif

#define	vec_length	(max_sections+1)
#define	vM	(M+1)				/* number of vertices in vertical direction   = M + 1      */
#define	vL	(L+1)				/* number of vertices in horizontal direction = L + 1      */

#if M==0
#define		CS     /* define CS for column states when M=0 */
#define max_sections 1000
#define max_tspans 1000
#endif

#define	vec_ent	double		/*vector entries will be of type defined here*/
#define	mat_ent	double		/*matrix entries will be of type defined here*/

#define	newline	printf("\n")	/* Allows the use of newline; in place of printf("\n"); 		  */
#define	VALID	1		/* This simply allows statements such as "if(something == VALID) {...}    */
#define	NOTVALID 	!VALID		/* This simply allows statements such as "if(something == NOTVALID) {...} */


/***************************************************************************/
/****************** Structures used in this program ************************/
/***************************************************************************/


struct hinge_span {	/* This data structure is used to store information about two-spans */
	/*unsigned long int 	inorder */		/* no longer needed as it is implicitely defined */
	unsigned long int 	outorder;	/* the section number of the righthand section   */	
	unsigned int 	edgecount;	/* the number of edges in the hinge and first section		 */
	struct hinge_span 	*nexthinge;		/* pointer to another hinge_span data structure  */
	unsigned int hedges[vM*vL];

}	*first_hinge_span[vec_length], *current_hinge_span[vec_length]; /* these are global */

	/* first_hinge_span is an array of pointers to hinge_span data structures where each pointer in */
	/* in the array points to the first element of a linked list of hinge_span data structures	 */
	/* The array index implicely defines the "inorder" (or the lefthand section of the two_span)     */

	/* current_hinge_span is a dummy pointer used for traversing the linked lists mentioned above    */
	/* an array of dummy pointers is needed because there is more than one linked list to consider   */

			    


/***************************************************************************/
/****************** end of structures **************************************/
/***************************************************************************/


/***************************************************************************/
/*************** Global variables used in this program *********************/
/***************************************************************************/

unsigned int		ordertemplate[2][vM][vL];	/* The first index is for side (left = 0 ; right = 1) */
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

//if x_0 is already known/defined, set it to connectivity_inv. If not, declare a variable connectivity_inv, and it will be found via rtflsp+power method
#if defined(thex)
double			connectivity_inv=thex;
#else
double			connectivity_inv;
#endif

unsigned int reordertemplate[2][vM][vL];	//holds the re-ordered ordertemplate to reduce number of unique sections.

unsigned short int hingestatus[vM][vL];		//keeps track of which vertices in the hinge are occupied. 0=empty, 1=filled.
unsigned short int alreadyentered[vM][vL];	//keeps track of where we've already initially entered. 0=no, 1=yes.

unsigned short int colhingeedges[M][vL]; //keeps track of which column (i/M) edges in the hinge are occupied. 0=empty, 1=filled.
unsigned short int rowhingeedges[vM][vL]; //keeps track of which row (j/L) edges in the hinge are occupied. 0=empty, 1=filled.

unsigned long int sectionkey[vec_length];	//will hold valid (re-ordered) section nums. index is the key num
									//note: sectionkey[0] is garbage. Key numbers start at 1.

unsigned long int dupcounter=0;	//will count the number of duplicate 2-spans created and not stored.



/* note: The following are also global: first_hinge_span, current_hinge_span */
/* They are discussed in the structures section above */

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

void            	enterhinge(int, int, int, int (*)[3]);
void            	leavehinge(int, int, int, int (*)[3]);
void            	rowedges(int, int, int (*)[3]);
void            	coledges(int, int, int (*)[3]);

/****  two_span validity verifier functions  ******************* *****************************/

unsigned short int	LFlag( int (*)[3]);
//unsigned short int	RFlag( int (*)[3]);
unsigned short int	noncrossing(int, int, int, int, int);	/* used in 0_L cases only   */

/****  one_span/section recorder functions  *************************************************/

void            	recordtemplate( int (*)[3]);
unsigned long int	num_section(int);
unsigned long int	num_section_norder(int);
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
void			qcksrt(unsigned long int,unsigned long int []);
void			qcksrtII(unsigned long int,unsigned long int [],unsigned long int []);

/*******************************************************************************************/

void				get_kspan(unsigned short int k);
unsigned short int	entsecval(unsigned short int side);
unsigned short int	entseccon(unsigned short int side);
//double			avgtwostrings(double, double, double,double);

/**** Added Functions  ********************************************************************/
double			sumtoone(double,double,double);
double			transprobcheck(double, double, double);
double			prob_check(double,double,double);

void				findsection(double, double, double);
void				printsection(int);
void				probkspan(double, double, double);

void 			fillreordertemplate(int, int);	//fills reordertemplate appropriately based on ordertemplate
void				printordtemp();
void				printreordtemp();



/***************************************************************************/
/******************* end of functions used in this program *****************/
/***************************************************************************/


int main(void) {
/**** variable used in main ************************************************/

	/*int	           	(*pointordNum)[3]; *//*pointer to an array of int*/
	int				ordNum[3],side = 0;
	int             	i, j;	   /*used in for loops */

	double			kappa;
	double			Alpha;		/* Alan's Alpha as defined in the notes provided by C. Soteros */
	double			Beta;		/* Beta as defined in the notes provided by C. Soteros */
	double			Ephi_Ephi;	//used in calculating Alan's alpha.
	double			LmultR;		/* inner product of left and right eigenvectors */

/**** functions called by main *********************************************/


	clock_t         clock(void);
	void            conv_to_array(void);

/**************************************************************************/
/**** start of main program ***********************************************/
/**************************************************************************/

	system("clear"); /* clears the screen */
	system("date");  /* prints the date and time */

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
		sectionkey[i]=0;
	}
	printf("space created for hinge structures and sectionkey initialized\n");

	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			hingestatus[i][j]=0;
			alreadyentered[i][j]=0;
		}
	}

	for(i=0; i<=M-1; i++){
		for(j=0; j<=L; j++){
			colhingeedges[i][j]=0;
		}
	}

	for(i=0; i<=M; i++){
		for(j=0; j<=L-1; j++){
			rowhingeedges[i][j]=0;
		}
	}



	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if( !(i==M && j==L) ){	//don't need to enter at last spot
				printf("entering at i=%d,j=%d\n", i, j);
				enterhinge(i, j, side, &ordNum);
			}
			alreadyentered[i][j]=1;	//don't need to explore any 2-spans that enter/exit here.
		}
	}

/*************************************************************************************************/
	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	newline; /* results are printed to the terminal */
	printf("There are a total of %9lu potential two-span patterns.\n", --num_2_spans);
	printf("There are a total of %9lu   valid   two-span patterns.\n", valid_2_spans);
/************************************************************************************************/

	printf("converting linked list to arrays...\n");
	conv_to_array();
	printf("done converting to arrays.\n\n");

	unsigned long int totalsecs=0;
	unsigned long int total2spans=0;
	for(i=1; i<=max_sections; i++){
		if(num_outsections[i]>0){
/*			printf("section %d is valid, it connects to %d other sections: \n", i, num_outsections[i]);
			printf("OGsecnum of %d is %d\n", i, sectionkey[i]);
			for(j=1; j<=num_outsections[i]; j++){
				printf("%d ", tspans_outsection[i][j]);
			}
			printf("\n");
			printsection(i);
*/			totalsecs++;
			total2spans=total2spans + num_outsections[i];
		}
	}
	printf("Check: Total sections=%lu, total 2spans=%lu\n", totalsecs, total2spans);

	//Ham version runs check here, for irreduciblity and "actual" valid 2span/sections counts

	L_Evector[0] = vecalloc(0,max_sections); /*global*/
	L_Evector[1] = vecalloc(0,max_sections); /*global*/ /*temp*/
	R_Evector[0] = vecalloc(0,max_sections); /*global*/
	R_Evector[1] = vecalloc(0,max_sections); /*global*/ /*temp*/

#if defined(thex)
	max_eval_LRvec(connectivity_inv); /*rtflsp not necessary once connectivity is known*/
#else
	printf("Starting false position method to determine connective constant.\n");
	#if defined(force)
		printf("\n=> %f is fugacity which yields 1 as maximum eigenvalue (after multiplying matrix by e^%f)\n",
					connectivity_inv=rtflsp(&max_eval_LRvec,x_L,x_R, 0.00000001), fval);
		printf("fugacity in sci. notation = %e\n", connectivity_inv);
	#else
		printf("\n=> %f is fugacity which yields 1 as maximum eigenvalue\n",
					connectivity_inv=rtflsp(&max_eval_LRvec,x_L,x_R, 0.00000001));
	#endif
		printf("\n=> %f is the connective constant.\n",1/connectivity_inv);
#endif

	//calculate LmultR
	LmultR = 0;
	for (i=1;i<= max_sections;i++){
		LmultR += (L_Evector[0][i])*(R_Evector[0][i]);
	}

/*
//CHECK IF LEFT EIGENVECTOR IS A LEFT EIGENVECTOR
	//initialize resultant vectors to zeros
	for(i=1; i<=max_sections; i++){
		L_Evector[1][i] = 0;
		R_Evector[1][i] = 0;
	}
	double pwfugacity;
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

	//Alan's version...not correct
	Alpha = connectivity_inv*connectivity_inv*connectivity_inv*connectivity_inv*exp(4*fval)*Ephi_Ephi/Beta;

	/*may now free L_Evector[1] this has yet to be implemented */
	/*may now free R_Evector[1] this has yet to be implemented */
//	printf("LmultR=%f\n", LmultR);
//	printf("L_Evector[0][1]=%f\n", L_Evector[0][1]);
//	printf("R_Evector[0][1]=%f\n", R_Evector[0][1]);

	newline;
	printf("Alan's Alpha = %f = %e\n",Alpha,Alpha);
 	printf("Beta  = %f = %e\n\n",Beta,Beta);

	printf("Final: Total sections=%lu, total 2spans=%lu\n", totalsecs, total2spans);
	printf("max_sections was set = %d\n", max_sections);
	printf("\n");

	printf("There are a total of %lu potential two-span patterns.\n", --num_2_spans);
	printf("There are a total of %lu   valid   two-span patterns.\n", valid_2_spans);
	printf("Final number of saved 2-spans=%lu\n", total2spans);
	printf("Number of duplicate 2-spans generated that were not recorded=%lu\n", dupcounter);
//	printf("L_Evector[phi]*R_Evector[phi] = %f = %e\n\n",Ephi_Ephi,Ephi_Ephi);

/************************************************************************************************/	
/*showmenu*/ /*Needs to be implemented*/

	kappa = connectivity_inv;
	/*system("clear");*/

//USE TO CHECK IF EXPECTED NUMBER OF OCCURRENCES ADD UP TO ONE (after not dividing by beta)
//	double sumofprobs = sumtoone(Beta, LmultR, kappa);
//	printf("sum of probs is %f\n", sumofprobs);

//USE TO CHECK IF PROB OF ALL SECTIONS (SS) ADD UP TO ONE.
//	double totprobs = prob_check(Beta, LmultR, kappa);
//	printf("totprobs is %f\n", totprobs);

//USE TO CHECK IF THE TRANSITION PROBABILITES ADD UP TO ONE FOR EACH SECTION
//	double transprobs = transprobcheck(Beta, LmultR, kappa);
//	printf("transprobs is %f\n", transprobs);

//USE TO FIND THE SECTION WITH THE EXPECTED NUMBER OF OCCURRENCES
//ONLY USE IN SECTION-SECTION MODE
//	findsection(Beta, LmultR, kappa);

//USE IF YOU WANT TO TAKE IN POLYGON FILES AND FIND THE SUM OF THE EXPECTED NUMBER OF OCCURRENCES (ALSO PRINTS PROBABILITY) OF
//THE K-SPAN CREATED BY REMOVING THE ENDHINGES OF THE POLYGON.
//	probkspan(Beta, LmultR, kappa);

//USE IF YOU WANT TO ENTER PATTERNS AND GET THE EXPECTED NUMBER OF OCCURRENCES OF A K-SPAN
//	Enum_Enum = enterpats(Beta, LmultR, kappa);

//USE IF YOU WANT TO CALCULATE THE EXP NUM OCCUR OF ALL SECTIONS WITH ONLY 2 EDGES (MUST BE IN SS MODE)
//	double expoccurtwostrings = avgtwostrings(Beta, LmultR, kappa, fval);
//	printf("exp num of occur of two string sections = %f = %e\n", expoccurtwostrings, expoccurtwostrings);





/*
	printf("sectionkey is this:\n");
	i=1;
	while(i<=max_sections){
		if(sectionkey[i]!=0){
			printf("sectionkey[%d]=%d\n", i, sectionkey[i]);
		}
		i++;
	}
*/

	system("date");  /* prints the date and time */
	printf("\nProgram Complete.\n\n");

	return 0;

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
enterhinge(int i, int j, int side, int (*pointordNum)[3])
{
	if( !(alreadyentered[i][j]==1 && side==0) ){
	//	printf("entered hinge (i=%d, j=%d, side=%d\n", i, j, side);
		hingestatus[i][j]=1;
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

		leavehinge(i, j, side ^ 1, pointordNum);
			/* the SAW leaves the hinge out the opposite side, hence (side ^ 1) is passed */
			/* if side = 1 then side ^ 1 = 0 */
			/* if side = 0 then side ^ 1 = 1 */
		rowedges(i, j, pointordNum);
			/* horizontal edges are explored */
		coledges(i, j, pointordNum);
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
		hingestatus[i][j]=0;
	//	printf("finished enterhinge\n");
	}
	return;

}				/* end of function enterhinge */

/***************************************************************************/

void
leavehinge(int i, int j, int side, int (*pointordNum)[3])
/* called by: rowedges,coledges,enterhinge */
{
	if( !(alreadyentered[i][j]==1 && side==0) ){
	//	printf("left hinge (i=%d, j=%d, side=%d\n", i, j, side);
		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		ordertemplate[side][i][j] = (*pointordNum)[side];	/* record the section
									 * edge (i.e. the
									 * leaving edge) */
		(*pointordNum)[side]++;
			/* the edge number is incremented so that the next edge gets a higher number */

		if (side == 0 && (*pointordNum)[1] > 1) {
	 	/* if the leaving edge is on the left and there is at least one edge on the right section */
		/* this may be a valid two span */
			if (LFlag(pointordNum)/* && RFlag(pointordNum)*/){
			/* if the two span connects to phi on the left and on the right then do the following*/

	//			printordtemp();
				fillreordertemplate((*pointordNum)[0]-1, (*pointordNum)[1]-1);
	//			printreordtemp();
	//			printf("\n");

				valid_2_spans++; /*This is a valid 2 span so it is counted. Includes duplicates. */
				recordtemplate(pointordNum); /*the information contained in ordertemplate is recorded*/
				//duplicate check is performed in recordtemplate.
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
				if (hingestatus[ii][jj] == 0) {
				/* if the vertex is not occupied then explore entering the hinge at this vertex */
				/* this maintains the self avoiding nature of the walk */
					enterhinge(ii, jj, side, pointordNum);
				}
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
	//		printf("finished leavehinge\n");
	}
	return;
}
/***************************************************************************/

void
rowedges(int i, int j, int (*pointordNum)[3])
{
//	printf("rowedges called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(hingestatus[i][j-1]==0){
			hingestatus[i][j-1]=1;
			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			rowhingeedges[i][j-1]=1;
			leavehinge(i, j-1, 0, pointordNum);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(i, j-1, 1, pointordNum);
				/* exit the hinge on side 1 with j moved*/
			rowedges(i, j-1, pointordNum);
				/* explore horizontal moves within the hinge */
			coledges(i, j-1, pointordNum);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i][j-1]=0;
			rowhingeedges[i][j-1]=0;
		}
	}
	if(j<L){
		if(hingestatus[i][j+1]==0){
			hingestatus[i][j+1]=1;
			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			rowhingeedges[i][j]=1;
			leavehinge(i, j+1, 0, pointordNum);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(i, j+1, 1, pointordNum);
				/* exit the hinge on side 1 with j moved*/
			rowedges(i, j+1, pointordNum);
				/* explore horizontal moves within the hinge */
			coledges(i, j+1, pointordNum);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i][j+1]=0;
			rowhingeedges[i][j]=0;
		}
	}
//	printf("finished rowedges\n");
	return;
}				/* end of function rowedges */

/***************************************************************************/

void
coledges(int i, int j, int (*pointordNum)[3])
/* called by rowedges,coledges */
{
//	printf("coledges called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(hingestatus[i-1][j]==0){
			hingestatus[i-1][j]=1;
			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			colhingeedges[i-1][j]=1;
			leavehinge(i-1, j, 0, pointordNum);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(i-1, j, 1, pointordNum);
				/* exit the hinge on side 1 with j moved*/
			rowedges(i-1, j, pointordNum);
				/* explore horizontal moves within the hinge */
			coledges(i-1, j, pointordNum);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i-1][j]=0;
			colhingeedges[i-1][j]=0;
		}
	}
	if(i<M){
		if(hingestatus[i+1][j]==0){
			hingestatus[i+1][j]=1;
			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			colhingeedges[i][j]=1;
			leavehinge(i+1, j, 0, pointordNum);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(i+1, j, 1, pointordNum);
				/* exit the hinge on side 1 with j moved*/
			rowedges(i+1, j, pointordNum);
				/* explore horizontal moves within the hinge */
			coledges(i+1, j, pointordNum);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i+1][j]=0;
			colhingeedges[i][j]=0;
		}
	}
//	printf("finished coledges\n");
	return;
}				/* end of function coledges */

/***************************************************************************/



/***************************************************************************/
void
recordtemplate( int (*pointordNum)[3])
/* This function records the pertinent information about a two-span */
/* which has been verified as being valid into a linked list */
{
	int i,j;

	unsigned long int	inNum = num_section_norder(0);
//	printf("OGinNum=%d\n", inNum);
	//Change inNum to it's key number.
	i=1;
	int found=0;
	while(sectionkey[i]!=0){
		if(sectionkey[i]==inNum){
			found=1;
			inNum=i;
			break;
		}
		i++;
	}
	if(found==0){	//need to record inNum in sectionkey
		sectionkey[i]=inNum;
		inNum=i;
	}
//	printf("NEWinNum=%d\n", inNum);

	unsigned long int	outNum = num_section_norder(1);
//	printf("OGoutNum=%d\n", inNum);
	//reset outNum to it's key number.
	i=1;
	found=0;
	while(sectionkey[i]!=0){
		if(sectionkey[i]==outNum){
			found=1;
			outNum=i;
			break;
		}
		i++;
	}
	if(found==0){	//need to record inNum in sectionkey
		sectionkey[i]=outNum;
		outNum=i;
	}
//	printf("NEWoutNum=%d\n", inNum);

	unsigned short int	Num_edges= (*pointordNum)[0] + (*pointordNum)[2] - 1;
				/*number of edges = edges in first section + edges in hinge */
				/*the correction (-1) is necessary since (*pointordNum)[0] */
				/*is the number that the next edge in the first section will */
				/*get */
	int temp_hedges[vM*vL];
	int index=0;
	int edgenum=1;


	for(i=0; i<=vM*vL-1; i++){
		temp_hedges[i]=0;
	}

	for(i=0; i<=vM-1; i++){
		for(j=0; j<=L-1; j++){
			if(rowhingeedges[i][j]==1){
				temp_hedges[index]=edgenum;
				index++;
			}
			edgenum++;
		}
	}
	for(i=0; i<=M-1; i++){
		for(j=0; j<=vL-1; j++){
			if(colhingeedges[i][j]==1){
				temp_hedges[index]=edgenum;
				index++;
			}
			edgenum++;
		}
	}
//	printf("there were %d edgenums\n", edgenum-1);
/*	printf("temp_hedges is: ");
	for(i=0; i<=vM*vL-1; i++){
		if(temp_hedges[i]==0){
			break;
		}
		printf("%d ", temp_hedges[i]);
	}
	printf("\n");
*/


	//check if this 2-span configuration has already been recorded
	int duplicate=0;
	int match;

	struct hinge_span* test = first_hinge_span[inNum];
	while(test->nexthinge != NULL && duplicate==0){
//		printf("looking through linked list\n");
		test = test->nexthinge;
		if(test->outorder==outNum){	//if same 2 sections
//			printf("found one that needs to test\n");
			match=1;	//test if the hinge edges match
			for(i=0; i<=vM*vL-1; i++){
				if(temp_hedges[i]!=test->hedges[i]){ //if an edge doesn't match
					match=0;
					break;
				}
			}
			if(match==1){
				duplicate=1;
				dupcounter++;
				break;
			}
		}
	}
//	printf("duplicate=%d\n", duplicate);
	if(duplicate==0){
//		printf("actually recording\n");
		(*current_hinge_span[inNum]).nexthinge = newhinge(); /* add a newhinge to the linked list */

		current_hinge_span[inNum] = (*current_hinge_span[inNum]).nexthinge; /* let current_hinge_span[inNum] point to newhinge */
		(*current_hinge_span[inNum]).outorder = outNum;	/* record the ordering on the out side */
		(*current_hinge_span[inNum]).edgecount = Num_edges;	/* record the number of edges in first section and in hinge */
		for(i=0; i<=vM*vL-1; i++){
			(*current_hinge_span[inNum]).hedges[i]=temp_hedges[i];
		}

/*		printf("Recorded 2-span. Section %d to %d, with hedges: ", inNum, outNum);
		for(i=0; i<=vM*vL-1; i++){
			printf("%d ", (*current_hinge_span[inNum]).hedges[i]);
		}
		printf("\n");
*/	}
//	printf("finished recording\n");
	return;
}



#include "pw_meth_ss_LRvec_fcheck.c" /* for applying the power method to the section-section matrix */
#include "../analysis/statistics/get_Beta_ssCheck.c"	/* for calculating Beta */


#if defined(CS)
#include "../topology/LFlag_0.c"		/* include if M = 0  */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the left and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

//#include "../topology/RFlag_0.c"		/* include if M = 0  */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the right and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#include "../utils/noncrossing.c"	/* include if M = 0 */
/* This function takes arguements (int side, int a, int b, int c, int d) */
/* and returns 1 if the walk connecting a and b does not cross */
/* the walk connecting c and d it returns 0 otherwise*/
#include "../sections/cstatenum.c"
#include "../sections/cstatenum_norder.c"
#else
#include "../topology/LFlag_norder.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the left and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

//#include "../topology/RFlag.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the right and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#if vM*vL<7
#include "../sections/Num_section_6V.c"  /*unsigned long int num_section(int side)*/
/* This function takes the arguement (int side) and assigns a unique number */
/* to the section on side = side */
#include "../sections/Num_section_6V_nonordered.c"
#elif vM*vL<9
#include "../sections/Num_section_8V.c"
#include "../sections/Num_section_8V_nonordered.c"
#elif vM*vL<11
#include "../sections/Num_section_10V.c"
#include "../sections/Num_section_10V_nonordered.c"
#elif vM*vL<13
#include "../sections/Num_section_12V.c"
#include "../sections/Num_section_12V_nonordered.c"
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
	for(i=0; i<=vM*vL-1; i++){
		nextnewhinge->hedges[i]=0;
	}
	return nextnewhinge;
}
/***************************************************************************/
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
			printf("\n%6lu:%6lu : edges=%3lu : two span# =%3lu\n", section_1, tspans_outsection[section_1][nth_tspan], tspans_edges[section_1][nth_tspan], tspans_nrr[section_1][nth_tspan]);
			result = 1;
		}
	}
	if (result == 0) {
		printf("This section does not connect to the previous section.\n");
	}
	return (result);
}

void fillreordertemplate(int ledges, int redges){ //fills reordertemplate appropriately based on ordertemplate
	int i, j, ii, jj, s;
	int firstentry, secondentry;
	int connectingedge;

//	printf("fillreordertemplate called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			reordertemplate[0][i][j]=0;
			reordertemplate[1][i][j]=0;
		}
	}

//	printf("initialized\n");
//	printf("ledges=%d, redges=%d\n", ledges, redges);
	for(s=0; s<=1; s++){	//for each side
		firstentry=1;
		for(i=0; i<=M; i++){	
			for(j=0; j<=L; j++){
				if(ordertemplate[s][i][j] % 2 != 0 && reordertemplate[s][i][j]==0){	//if odd
//					printf("odd edge found in (%d,%d)\n", j, i);
					reordertemplate[s][i][j] = firstentry;
//					printf("recorded edge %d in reorder\n", firstentry);
					if(ordertemplate[s][i][j]==1){	//if 1, then have to find the location of ledges or redges.
						if(s==0){
							connectingedge=ledges;
						}
						else{
							connectingedge=redges;
						}
					}
					else{	//if not 1, then find location the edgenum-1
						connectingedge=ordertemplate[s][i][j]-1;
					}

					if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
						if(s==0){
							secondentry=ledges;
						}
						else{
							secondentry=redges;
						}
					}
					else{
						secondentry=firstentry+1;
					}
					for(ii=0; ii<=M; ii++){
						for(jj=0; jj<=L; jj++){
							if(ordertemplate[s][ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
								reordertemplate[s][ii][jj] = secondentry;
//								printf("putting %d into that location\n", secondentry);
								if(firstentry==1){
									firstentry++;
								}
								else{
									firstentry += 2;
								}
								break;
							}
						}
					}
				}
				else if(ordertemplate[s][i][j]>0 && reordertemplate[s][i][j]==0){	// if even
//					printf("even edge found in (%d,%d)\n", j, i);
					reordertemplate[s][i][j] = firstentry;
//					printf("recorded edge %d in reorder\n", firstentry);
					if(s==0 && ordertemplate[s][i][j]==ledges){	//if max, then have to find the location of edge 1.
						connectingedge=1;
					}
					else if(s==1 && ordertemplate[s][i][j]==redges){
						connectingedge=1;
					}
					else{	//if not max, then find location the edgenum+1
						connectingedge=ordertemplate[s][i][j]+1;
					}

					if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
						if(s==0){
							secondentry=ledges;
						}
						else{
							secondentry=redges;
						}
					}
					else{
						secondentry=firstentry+1;
					}
					for(ii=0; ii<=M; ii++){
						for(jj=0; jj<=L; jj++){
							if(ordertemplate[s][ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
								reordertemplate[s][ii][jj] = secondentry;
								if(firstentry==1){
									firstentry++;
								}
								else{
									firstentry += 2;
								}
								break;
							}
						}
					}
				}
			}
		}	
	}
}

void printordtemp(){
	int i, j;
	printf("printing ordertemplate. Side 0:\n");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("%d ", ordertemplate[0][i][j]);
		}
		printf("\n");
	}

	printf("ordertemplate side1:\n");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("%d ", ordertemplate[1][i][j]);
		}
		printf("\n");
	}
}

void printreordtemp(){
	int i, j;
	printf("printing reordertemplate. Side 0:\n");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("%d ", reordertemplate[0][i][j]);
		}
		printf("\n");
	}

	printf("reordertemplate side1:\n");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("%d ", reordertemplate[1][i][j]);
		}
		printf("\n");
	}
}







#include "../sections/findsection.c"
#include "../sections/printsection.c"


#include "../analysis/statistics/sumtoone_ss.c"
#include "../analysis/statistics/prob_check_ss.c"
//#include "../sections/enterpatstspan_ss.c"
#include "../sections/enterpats_ss.c"
#include "../analysis/statistics/transprobcheck_ss.c"
#include "../analysis/statistics/probkspan_ss_sectionkey.c"
//#include "../utils/avgtwostrings.c"

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






	









