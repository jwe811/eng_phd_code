//HAMILTONIAN VERSION

//This program will sample Hamiltonian SAPs of a certain span uniformly at random (by span).
//uses 2SAP_TMcalcHam_PrintEvectors.c, which prints out eigenvectors to:
//2SAP_L_EvectorHam_TS_L%dM%d.txt
//2SAP_R_EvectorHam_TS_L%dM%d.txt
//These files must be in same folder.

//compliled with: gcc -lm -Wall -O3 -o 2SAP_MCsample_Ham.o 2SAP_MCsample_Ham.c
//ran using ./2SAP_MCsample_Ham.o

//sample 2SAPs will go to files: MC2SAPsHamL%dM%dspan%dnum%d.txt

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include "../../include/marsaglia.h"


/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/
//for now, need L,M>0, totalspan>0,

#define	L 3				// number of horizontal edges
#define	M 1				// number of vertical edges
#define	totalspan 250		// span of sampled SAPs
#define	samplesize 1000000	// number of samples desired
#define	runnum 1
#define 	seednum 931250

#define	maxpolys 10000	// number of 2SAPs per file.

/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/

//Automatically set max_sections, based on previous runs.
#if M==1 && L==1
#define max_sections 8

#elif M==1 && L==2
#define max_sections 73
#define max_keynum 152
#define max_tspans 490
#define dom_evalue 5.534148126030995

#elif M==1 && L==3
#define max_sections 742
#define max_keynum 4048
#define max_tspans 37454
#define dom_evalue 24.378235811209002

#elif M==1 && L==4
#define max_sections 9309
#define max_keynum 108386
#define max_tspans 2598620
#define dom_evalue 97.848097677835298

#elif M==1 && L==5
#define max_sections 138038

#elif M==2 && L==2
#define max_sections 2619
#define max_keynum 21540
#define max_tspans 495792
#define dom_evalue 53.686533141903084

#elif M==2 && L==3
#define max_sections 138322
#endif

#define	vec_length	(max_sections+1)
#define	vM	(M+1)				/* number of vertices in vertical direction   = M + 1      */
#define	vL	(L+1)				/* number of vertices in horizontal direction = L + 1      */

#if M==0
#define		CS     /* define CS for column states when M=0 */
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
	unsigned long int	inorder2;		// section num of lefthand section of poly2
	unsigned long int 	outorder;		// section num of righthand section of poly1
	unsigned long int	outorder2;	// section num of righthand section of poly2

	unsigned int		hedges[vM*vL];	//contains info about which hinge-edges are filled

	int		*start[3];	//start[0][i] contains the starting x-coordinate of the i'th walk.
						//start[1][i] contains the starting y-coordinate of the i'th walk.
						//start[2][i] contains the starting z-coordinate of the i'th walk.
	int		*end[3];		//end[0][i] contains the ending x-coordinate of the i'th walk.
						//end[1][i] contains the ending y-coordinate of the i'th walk.
						//end[2][i] contains the ending z-coordinate of the i'th walk.
	int		*walks[vM*vL+1];	//walks[j][i] contains the j'th step in the i'th walk.	rest zero's.
	unsigned int	num_of_walks;	//number of walks in this 2-span.

	int		*start2[3];	//same but for 2nd poly
	int		*end2[3];		//same but for 2nd poly
	int		*walks2[vM*vL+1];	//same but for 2nd poly
	unsigned int	num_of_walks2;	//same but for 2nd poly

	struct hinge_span	*nexthinge;	// pointer to another hinge_span data structure

}	*first_hinge_span[vec_length], *current_hinge_span[vec_length]; /* these are global */

	/* first_hinge_span is an array of pointers to hinge_span data structures where each pointer	*/
	/* in the array points to the first element of a linked list of hinge_span data structures	*/
	/* The array index implicely defines the "inorder" (or the lefthand section of the two_span)  	*/

	/* current_hinge_span is a dummy pointer used for traversing the linked lists mentioned above	*/
	/* an array of dummy pointers is needed because there is more than one linked list to consider	*/


struct endhinge {
	unsigned long int	sec2;			//section number of poly2
	unsigned int		hedges[vM*vL];	//contains info about which hinge-edges are filled
	unsigned short int	side;		//0 if a left endhinge, 1 if a right endhinge

	int		*start[3];	//start[0][i] contains the starting x-coordinate of the i'th walk.
						//start[1][i] contains the starting y-coordinate of the i'th walk.
						//start[2][i] contains the starting z-coordinate of the i'th walk.
	int		*end[3];		//end[0][i] contains the ending x-coordinate of the i'th walk.
						//end[1][i] contains the ending y-coordinate of the i'th walk.
						//end[2][i] contains the ending z-coordinate of the i'th walk.
	int		*walks[vM*vL+1];	//walks[j][i] contains the j'th step in the i'th walk.	rest zero's.
	unsigned int	num_of_walks;	//number of walks in this endhinge

	int		*start2[3];		//same but for 2nd poly
	int		*end2[3];			//same but for 2nd poly
	int		*walks2[vM*vL+1];	//same but for 2nd poly
	unsigned int	num_of_walks2;	//same but for 2nd poly

	struct endhinge* nextendhinge;// pointer to the next endhinge structure in the linked list
};

struct endhinge* firstendhinge[vec_length];		/* Array of pointers to endhinge structures, where each	*/
										/* pointer in the array points to the first element in		*/
										/* the linked list of endhinge structures				*/
struct endhinge* currentendhinge[vec_length];	/* Array of dummy pointers that are used to traverse		*/
										/* the endhinge linked lists.							*/



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
unsigned int		ordertemplate2[2][vM][vL]; //same but for 2nd poly

unsigned long int 	num_2_spans = 1;	/* This variable keeps track of the total number of potential two spans   */

unsigned long int	valid_2_spans = 0;	/* This is used to count the number of two spans which are actually valid */

////////////////////// arrays which hold the tspan info ////////////////////////////////

unsigned long int	num_outsections[max_keynum+1];	/* For a given section (which has a number section_num assigned to it)*/
										/* num_outsections[section_num] is the number of two spans with this  */
										/* section as a first section                                         */
 
unsigned long int	*t_outsection[max_keynum+1];	/* An array of pointers to integers 					*/
										/* records the "outsection" of the two-span				*/
										/* t_outsection[insection][nth_tspan]				*/
										/* nth_tspan = nth two-span with firstsection = insection	*/

unsigned long int	*t_nrr[max_keynum+1];		/*An array of pointers to integers*/
										/*records the tspan "number"      */
										/* defined by the array index */
										/*example:		*/
										/*t_nrr[1][1] gives the two-span number */
										/* of the two-span with first section = 1 */
										/* and the second section is the first such */
										/* section which connects with first section = 1 */

unsigned int		*t_num_walks[max_keynum+1];	//t_num_walks[insection][nth_tspan] holds the number of walks in this tspan

int				***t_start[max_keynum+1];	//t_start[insection][nth_tspan][i][nth_walk] holds the i-coordinate (0=x, 1=y, 2=z) of the start point of the nth_walk

int				***t_end[max_keynum+1];		//t_end[insection][nth_tspan][i][nth_walk] holds the i-coordinate (0=x, 1=y, 2=z) of the end point of the nth_walk

int				***t_walks[max_keynum+1];	//t_end[insection][nth_tspan][direc][nth_walk] holds the "direc'th" direction the nth_walk

unsigned int		*t_num_walks2[max_keynum+1];	//same but for poly2
int				***t_start2[max_keynum+1];	//same but for poly2
int				***t_end2[max_keynum+1];		//same but for poly2
int				***t_walks2[max_keynum+1];	//same but for poly2



////////////////////// arrays which hold the endhinge info ////////////////////////////////////

unsigned long int	num_left_endhinges[max_keynum+1];
unsigned int		*Lend_num_walks[max_keynum+1];
int				***Lend_start[max_keynum+1];
int				***Lend_end[max_keynum+1];
int				***Lend_walks[max_keynum+1];

unsigned int		*Lend_num_walks2[max_keynum+1];
int				***Lend_start2[max_keynum+1];
int				***Lend_end2[max_keynum+1];
int				***Lend_walks2[max_keynum+1];



unsigned long int	num_right_endhinges[max_keynum+1];
unsigned int		*Rend_num_walks[max_keynum+1];
int				***Rend_start[max_keynum+1];
int				***Rend_end[max_keynum+1];
int				***Rend_walks[max_keynum+1];

unsigned int		*Rend_num_walks2[max_keynum+1];
int				***Rend_start2[max_keynum+1];
int				***Rend_end2[max_keynum+1];
int				***Rend_walks2[max_keynum+1];



////////////////

unsigned int reordertemplate[2][vM][vL];	//holds the re-ordered ordertemplate to reduce number of unique sections.
unsigned int reordertemplate2[2][vM][vL];	//same but for poly2

unsigned short int hingestatus[vM][vL];		//keeps track of which vertices in the hinge are occupied. 0=empty, 1=filled.
unsigned short int alreadyentered[vM][vL];	//keeps track of where we've already initially entered. 0=no, 1=yes.
unsigned short int alreadyentered2[vM][vL];	//same but for poly2

unsigned short int colhingeedges[M][vL]; //keeps track of which column (i/M) edges in the hinge are occupied. 0=empty, 1=filled.
unsigned short int rowhingeedges[vM][vL]; //keeps track of which row (j/L) edges in the hinge are occupied. 0=empty, 1=filled.

unsigned long int sectionkey[vec_length];	//will hold valid (re-ordered) section nums. index is the key num
									//note: sectionkey[0] is garbage. Key numbers start at 1.

unsigned long int sectionkey2SAP[max_keynum+1][2];

unsigned long int dupcounter=0;	//will count the number of duplicate 2-spans created and not stored.

int curstart[3][vM*vL];		//will hold the start-points of the walks being built
int curend[3][vM*vL];		//will hold the end-points of the walks being built
int curwalks[vM*vL+1][vM*vL];	//will hold the walks as they are being built.
						//curwalks[j][i] contains the j'th step in the i'th walk.	rest zero's.

int curstart2[3][vM*vL];			//same but for poly2
int curend2[3][vM*vL];			//same but for poly2
int curwalks2[vM*vL+1][vM*vL];	//same but for poly2

unsigned int num_walks=0;	//number of walks in the 2-span
unsigned int num_walks2=0;	//same but for poly2

//endhinge stuff
int EndOrdNum[2];		//contains info on # of edges in left section and endhinge
int EndOrdNum2[2];		//contains info on # of edges in left section and endhinge

int endtemplate[vM][vL];		//contains the section associated with the end hinge
int rendtemplate[vM][vL];	//holds the re-ordered endtemplate to reduce number of unique sections.
int endtemplate2[vM][vL];
int rendtemplate2[vM][vL];

int endtemplate3[vM][vL];	//same but for 2nd poly
int rendtemplate3[vM][vL];
int endtemplate4[vM][vL];
int rendtemplate4[vM][vL];

unsigned long int tot_right_endhinges=0;
unsigned long int tot_left_endhinges=0;
unsigned long int num_duplicate_right_endhinges=0;
unsigned long int num_duplicate_left_endhinges=0;
unsigned long int num_tspans=0;
unsigned long int num_duplicate_tspans=0;

unsigned long int num_potential_endhinges=0;
unsigned long int num_valid_endhinges=0;

unsigned long int num_nonlexi_endhinges=0;


//integer division rounds down in case of vM*vL is odd
int built_walks_start[vM*vL/2][3];
int built_walks_end[vM*vL/2][3];
int built_walks_direcs[vM*vL/2][vM*vL*(totalspan+1)];

int num_built_walks;

int built_walks_start2[vM*vL/2][3];	//2nd poly
int built_walks_end2[vM*vL/2][3];	//2nd poly
int built_walks_direcs2[vM*vL/2][vM*vL*(totalspan+1)];	//2nd poly

int num_built_walks2;	//2nd poly

char 			filename[100];	//this will hold the filename of the text file that we want to write to
FILE*			fp;			//file pointer. will point to the file that is to be written to.
unsigned long int filetotal=0;	//number of polys in the file 
unsigned long int filenum=1;		//number of polygon files created.

unsigned short int testvector[max_keynum+1];	//this will be the vector corresponding to the testsection.
unsigned short int temptestvector[max_keynum+1];	//used for the vector-matrix multiplication
unsigned short int testvector2[max_keynum+1];	//this will be the vector corresponding to the testsection.
unsigned short int temptestvector2[max_keynum+1];	//used for the vector-matrix multiplication





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
/* conv_endhinges_to_array	is only called by main and is declared in main */

/****  two_span builder functions  **********************************************************/

void            	enterhinge(int, int, int, int (*)[6], int);
void            	leavehinge(int, int, int, int (*)[6], int);
void            	rowedges(int, int, int (*)[6], int);
void            	coledges(int, int, int (*)[6], int);

void            	enterhinge2(int, int, int, int (*)[6], int);
void            	leavehinge2(int, int, int, int (*)[6], int);
void            	rowedges2(int, int, int (*)[6], int);
void            	coledges2(int, int, int (*)[6], int);

/****  two_span validity verifier functions  ******************* *****************************/

unsigned short int	LFlag2( int (*)[6]);
//unsigned short int	RFlag( int (*)[3]);
unsigned short int	noncrossing(int, int, int, int, int);	/* used in 0_L cases only   */

/****  one_span/section recorder functions  *************************************************/

void            	recordtemplate( int (*)[6]);
unsigned long int	num_section(int);
unsigned long int	num_section_norder(int);
unsigned long int	num_section_norder2(int);
unsigned short int	paircomp(unsigned short int, unsigned short int);
unsigned long int	choose(unsigned short int, unsigned short int);

/****  Memory allocation functions  ********************************************************/

struct hinge_span	*newhinge(unsigned int, unsigned int);
struct endhinge	*newendhinge(unsigned int, unsigned int);
int				*int_vecalloc(int low, int high);
unsigned int		*unsgn_vecalloc(int low, int high);
unsigned long int	*unsgnlong_vecalloc(int low, int high);
vec_ent			*vecalloc(int low, int high);
mat_ent			**matalloc(int rowlow, int highrow, int lowcol, int highcol);


/**** Added Functions  ********************************************************************/
void				findsection(double, double, double);
void				printsection(int);

void 			fillreordertemplate(int, int);	//fills reordertemplate appropriately based on ordertemplate
void 			fillreordertemplate2(int, int);	//fills reordertemplate appropriately based on ordertemplate
void				printordtemp();
void				printreordtemp();

void enterendhinge(int i, int j, int curlength);
void leaveendhinge(int i, int j, int curlength);
void endrowedges(int i, int j, int curlength);
void endcoledges(int i, int j, int curlength);

void enterendhinge2(int i, int j, int curlength2);
void leaveendhinge2(int i, int j, int curlength2);
void endrowedges2(int i, int j, int curlength2);
void endcoledges2(int i, int j, int curlength2);

unsigned short int LFlag_endhinge2();
void recordendtemplate();

void fillrendtemplate();
void fillendtemplate2();
void fillrendtemplate2();

void fillrendtemplate3();
void fillendtemplate4();
void fillrendtemplate4();


unsigned long int num_section_endhinge_norder(void);
unsigned long int num_section_endhinge_norder2(void);
unsigned long int num_section_endhinge_norder3(void);
unsigned long int num_section_endhinge_norder4(void);

void printbuiltwalks();
void printbuiltwalks_all();
void printbuiltwalks2_all();

void add_to_built_walks(unsigned long int secnum, int nth_tspan);
void add_right_endhinge(unsigned long int secnum, int nth_endhinge);

int reverse_direc(int direc);

void printtofile();

/***************************************************************************/
/******************* end of functions used in this program *****************/
/***************************************************************************/


int main(void) {

	unsigned int seed=seednum;
	initran_(&seed);

//	time_t t;
//	srand((unsigned) time(&t));

/**** variable used in main ************************************************/

	/*int	           	(*pointordNum)[3]; *//*pointer to an array of int*/
	int				ordNum[6],side = 0;
	int             	i,j,k,l;	   /*used in for loops */

/**** functions called by main *********************************************/

	clock_t         clock(void);
	void            conv_to_array(void);
	void            conv_endhinges_to_array(void);

/**************************************************************************/
/**** start of main program ***********************************************/
/**************************************************************************/

	system("clear"); /* clears the screen */
	system("date");  /* prints the date and time */

	printf("section-section version\n");

	printf("\nThe variables are as follows M: %i L: %i\n", M, L);


	clock();	/* Start clock to time program */

	ordNum[0] = 1;	/* first edge in section 0 gets numbered 1 (left side of 2-span)	*/
	ordNum[1] = 1;	/* first edge in section 1 gets numbered 1 (right side of 2-span)	*/
	ordNum[2] = 0;	/* to start there are no edges in the 2-span (edges in hinge)		*/
	ordNum[3] = 1;  // same as ordNum[0] but for the second polygon
	ordNum[4] = 1;  // same as ordNum[1] but for the second polygon
	ordNum[5] = 0;  // same as ordNum[2] but for the second polygon

	for (i=1;i<=max_sections; i++){
		current_hinge_span[i] = newhinge(1,1);		/* free up space of hinge structure */
		first_hinge_span[i] = current_hinge_span[i];	/* make first_hinge_span[i] point to this first one */
		currentendhinge[i] = newendhinge(1, 1);		/* free up space of endhinge structure */
		firstendhinge[i] = currentendhinge[i];		/* make firstendhinge[i] point to this first one */
		sectionkey[i]=0;
	}
	printf("space created for hinge structures, endhinge structures, and sectionkey initialized\n");

	//initializations//////////////////////////////////////////////////
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			hingestatus[i][j]=0;
			alreadyentered[i][j]=0;
			alreadyentered2[i][j]=0;
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

	num_walks=0;
	for(i=0; i<=vM*vL-1; i++){
		for(j=0; j<=2; j++){
			curstart[j][i]=-1;
			curstart2[j][i]=-1;
			curend[j][i]=-1;
			curend2[j][i]=-1;
		}
		for(j=0; j<=vM*vL; j++){
			curwalks[j][i]=0;
			curwalks2[j][i]=0;
		}
	}
	//initializations end///////////////////////////////////////////////	


	//create 2-spans
	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if( !(i==M && j==L) ){	//don't need to enter at last spot
				printf("2-span: entering at i=%d,j=%d\n", i, j);
				enterhinge(i, j, side, &ordNum, 0);
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
	printf("%lu recorded tspans\n", num_tspans);
	printf("%lu duplicate tspans that were 'valid', but not recorded\n", dupcounter);

	printf("\n");


////////////////////////////////////////////////////////////////////////


	//generate the endhinges
	//we are going to create "right endhinges" (entering from the left, can't leave to the right)
	//When we record these right endhinges, then we will make appropriate changes to also record it's "mirror" left endhinge

	//initialize stuff to zero.
	num_walks=0;
	num_walks2=0;
	for(i=0; i<=vM*vL-1; i++){
		for(j=0; j<=2; j++){
			curstart[j][i]=-1;
			curend[j][i]=-1;
			curstart2[j][i]=-1;
			curend2[j][i]=-1;
		}
		for(j=0; j<=vM*vL; j++){
			curwalks[j][i]=0;
			curwalks2[j][i]=0;
		}
	}

	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			hingestatus[i][j]=0;
			alreadyentered[i][j]=0;
			alreadyentered2[i][j]=0;
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

	EndOrdNum[0]=1;	//first edge in section will get the label 1
	EndOrdNum[1]=0;	//no edges currently in the endhinge
	EndOrdNum2[0]=1;	//first edge in section will get the label 1
	EndOrdNum2[1]=0;	//no edges currently in the endhinge

	//create endhinges
	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if( !(i==M && j==L) ){	//don't need to enter at last spot
				printf("Endhinge: entering at i=%d,j=%d\n", i, j);
				enterendhinge(i, j, 0);
			}
			alreadyentered[i][j]=1;	//don't need to explore any 2-spans that enter/exit here.
		}
	}

	printf("number of potential endhinges: %lu\n", num_potential_endhinges);
	printf("number of valid endhinges: %lu\n", num_valid_endhinges);

	printf("%lu right endhinges recorded\n", tot_right_endhinges);
	printf("%lu left endhinges recorded\n", tot_left_endhinges);
	printf("%lu non-lexi left endhinges\n", tot_left_endhinges);
	printf("%lu duplicate right endhinges\n", num_duplicate_right_endhinges);
	printf("%lu duplicate left endhinges\n", num_duplicate_left_endhinges);




//tests
/*
	for(i=1; i<=max_sections; i++){
		currentendhinge[i]=firstendhinge[i];
		while(currentendhinge[i]->nextendhinge!=NULL){
			currentendhinge[i]=currentendhinge[i]->nextendhinge;
			if(currentendhinge[i]->side==1){
				printf("Endhinge with sectionkey =%d,%lu. side=%d, num_walks=%d, num_walks2=%d\n", i, currentendhinge[i]->sec2, currentendhinge[i]->side, currentendhinge[i]->num_of_walks, currentendhinge[i]->num_of_walks2);
				for(j=0; j<=currentendhinge[i]->num_of_walks-1; j++){
					printf("walk1: %d. Start=(%d,%d,%d), Directions: ", j, currentendhinge[i]->start[0][j], currentendhinge[i]->start[1][j], currentendhinge[i]->start[2][j]);
					for(k=0; k<=vM*vL; k++){
						printf("%d ", currentendhinge[i]->walks[k][j]);
					}
					printf(". End=(%d,%d,%d)\n", currentendhinge[i]->end[0][j], currentendhinge[i]->end[1][j], currentendhinge[i]->end[2][j]);
				}
				for(j=0; j<=currentendhinge[i]->num_of_walks2-1; j++){
					printf("walk2: %d. Start=(%d,%d,%d), Directions: ", j, currentendhinge[i]->start2[0][j], currentendhinge[i]->start2[1][j], currentendhinge[i]->start2[2][j]);
					for(k=0; k<=vM*vL; k++){
						printf("%d ", currentendhinge[i]->walks2[k][j]);
					}
					printf(". End=(%d,%d,%d)\n", currentendhinge[i]->end2[0][j], currentendhinge[i]->end2[1][j], currentendhinge[i]->end2[2][j]);
				}
				printf("hedges = ");
				for(j=0; j<=vM*vL-1; j++){
					printf("%d ", currentendhinge[i]->hedges[j]);
				}
				printf("\n");
			}
		}
	}
*/




/*************************************************************************************************/
	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
/************************************************************************************************/




	printf("converting hinge_span linked lists to arrays...\n");
	conv_to_array();
	printf("done converting hinge_spans to arrays.\n\n");

	printf("converting endhinge linked lists to arrays...\n");
	conv_endhinges_to_array();
	printf("done converting endhinges to arrays.\n\n");



//tests

/*
	int l;
	for(i=1; i<=max_keynum; i++){
		for(j=1; j<=num_outsections[i]; j++){
			printf("\nsectionkey %d (%lu,%lu) connects to sectionkey %lu (%lu, %lu). It consists of %d walks in poly1, and %d walks in poly2.\n", i, sectionkey2SAP[i][0], sectionkey2SAP[i][1], t_outsection[i][j], sectionkey2SAP[t_outsection[i][j]][0], sectionkey2SAP[t_outsection[i][j]][1], t_num_walks[i][j], t_num_walks2[i][j]);
			printf("poly1:\n");
			for(k=0; k<=t_num_walks[i][j]-1; k++){
				printf("Walk %d. Start (%d,%d,%d). Directions: ", k, t_start[i][j][0][k], t_start[i][j][1][k], t_start[i][j][2][k]);
				for(l=0; l<=vM*vL; l++){
					printf("%d ", t_walks[i][j][l][k]);
				}
				printf(". End (%d,%d,%d)\n", t_end[i][j][0][k], t_end[i][j][1][k], t_end[i][j][2][k]);
			}
			printf("poly2:\n");
			for(k=0; k<=t_num_walks2[i][j]-1; k++){
				printf("Walk %d. Start (%d,%d,%d). Directions: ", k, t_start2[i][j][0][k], t_start2[i][j][1][k], t_start2[i][j][2][k]);
				for(l=0; l<=vM*vL; l++){
					printf("%d ", t_walks2[i][j][l][k]);
				}
				printf(". End (%d,%d,%d)\n", t_end2[i][j][0][k], t_end2[i][j][1][k], t_end2[i][j][2][k]);
			}

		}
	}


	for(i=1; i<=max_keynum; i++){
		printf("sectionkey %d has %lu left endhinges attached\n", i, num_left_endhinges[i]);
		for(j=1; j<=num_left_endhinges[i]; j++){
			printf("left endhinge %d has %d walk(s) in poly1 and %d walk(s) in poly2\n", j, Lend_num_walks[i][j], Lend_num_walks2[i][j]);
			printf("poly1:\n");
			for(k=0; k<=Lend_num_walks[i][j]-1; k++){
				printf("Walk %d. Start (%d,%d,%d). Directions: ", k, Lend_start[i][j][0][k], Lend_start[i][j][1][k], Lend_start[i][j][2][k]);
				for(l=0; l<=vM*vL; l++){
					printf("%d ", Lend_walks[i][j][l][k]);
				}
				printf(". End (%d,%d,%d)\n", Lend_end[i][j][0][k], Lend_end[i][j][1][k], Lend_end[i][j][2][k]);
			}
			printf("poly2:\n");
			for(k=0; k<=Lend_num_walks2[i][j]-1; k++){
				printf("Walk %d. Start (%d,%d,%d). Directions: ", k, Lend_start2[i][j][0][k], Lend_start2[i][j][1][k], Lend_start2[i][j][2][k]);
				for(l=0; l<=vM*vL; l++){
					printf("%d ", Lend_walks2[i][j][l][k]);
				}
				printf(". End (%d,%d,%d)\n", Lend_end2[i][j][0][k], Lend_end2[i][j][1][k], Lend_end2[i][j][2][k]);
			}
		}
	}

*/

	unsigned long int totalsecs=0;
	unsigned long int total2spans=0;
	for(i=1; i<=max_keynum; i++){
		if(num_outsections[i]>0){
			totalsecs++;
			total2spans=total2spans + num_outsections[i];
		}
//		printf("\n");
	}
	printf("Check: Total 2SAP sections=%lu, total 2SAP 2spans=%lu\n", totalsecs, total2spans);


//exit(1);




	//Ham version runs check here, for irreduciblity and "actual" valid 2span/sections counts

/////////CHECK IF THE TRANSMATRIX IS IRREDUCIBLE////////////////////////////
/////////ALSO GET ACTUAL NUMBER OF VALID 2-SPANS AND SECTIONS///////////////
/////////GET RID OF SECTIONS/2-SPANS WHICH CANT OCCUR IN A COMPACT POLY/////

	printf("Finding ACTUAL valid 2-spans and sections\n");


	unsigned long int testsection=1;	//section that is going to be "tested" to see if it eventually can connect to all other sections
//	unsigned long int testvector[vec_length];	//this will be the vector corresponding to the testsection. had to make this global for some reason.

	printf("\nSection that is being tested to see if it connects to all other sections is section=%lu\n\n", testsection);
	printf("It consists of sections %lu and %lu\n", sectionkey2SAP[testsection][0], sectionkey2SAP[testsection][1]);


	for(i=1; i<=max_keynum; i++) {
		if(i==testsection){	//initialize testvector be all zero, expect for the entry corresponding to the testsection
			testvector[i]=1;
		}
		else{
			testvector[i]=0;
		}
	}



//	unsigned long int temptestvector[vec_length];	//used for the vector-matrix multiplication. had to make this global for some reason




	//BE CAREFUL HERE BECAUSE OF PERIODICITY (in 2x2 case). NEED testpower TO BE BOTH EVEN AND ODD.
	//use an even (odd) power here
	unsigned long int testpower=100;	//power of the matrix that will be multiplied by testvector (steps in markov chain)
	printf("Performing check with test power (%lu) \n", testpower);
	for(k=1; k<=testpower; k++){
		for(i=1; i<=max_keynum; i++){
			temptestvector[i] = 0; //set temp=vector of zeros
		}

		for(i=1; i<=max_keynum; i++){
			for(j=1; j<=num_outsections[i]; j++){
				temptestvector[t_outsection[i][j]] += testvector[i]*t_num_walks[i][j]; //temp=test*matrix
			}
		}
		for(i=1; i<=max_keynum; i++){
			testvector[i]=temptestvector[i];	//set test=temp
		}
		for(i=1; i<=max_keynum; i++){
			if(testvector[i] != 0){
				testvector[i]=1;
//				printf("testvector[%d]=%d, \n", i, testvector[i]);
			}
		}
//		printf("iter=%d\n", k);
	}

	if(M==2 && L==2){
		//REPEAT WITH the the opposite of even (odd) TESTPOWER NOW
		testpower=testpower+1;
		printf("Now performing with odd test power (%lu)\n", testpower);

		for(i=1; i<=max_keynum; i++) {
			if(i==testsection){	//initialize testvector2 be all zero, expect for the entry corresponding to the testsection
				testvector2[i]=1;
			}
			else{
				testvector2[i]=0;
			}
		}



		for(k=1; k<=testpower; k++){
			for(i=1; i<=max_keynum; i++){
				temptestvector2[i] = 0; //set temp=vector of zeros
			}

			for(i=1; i<=max_keynum; i++){
				for(j=1; j<=num_outsections[i]; j++){
					temptestvector2[t_outsection[i][j]] += testvector2[i]*t_num_walks[i][j]; //temp=test*matrix
				}
			}
			for(i=1; i<=max_keynum; i++){
				testvector2[i]=temptestvector2[i];	//set test=temp
			}
			for(i=1; i<=max_keynum; i++){
				if(testvector2[i] != 0){
					testvector2[i]=1;
	//				printf("testvector2[%d]=%d\n", i, testvector2[i]);
				}
			}
	//		printf("iter=%d\n", k);
		}
	}


	unsigned long int numnoconnect=0;	//this will count the number of non-connected sections (to the testsection)
	unsigned long int numconnect=0;	//this will count the number of connected sections (to the testsection)
	static unsigned long int nonconnectingsections[max_keynum];	//this will hold the section numbers which DONT connect to testsection
	static unsigned long int connectingsections[max_keynum];		//this will hold the section numbers which DO connect to testsection
	//initialize these vectors to all zeros
	for(i=0; i<=max_keynum-1; i++){
		nonconnectingsections[i] = 0;
		connectingsections[i] = 0;
	}


//	printf("Section %d doesn't connect to sections: \n", testsection);
	int indx=0;	//index variable for nonconnectingsections[]
	int indx2=0;	//index variable for connectingsections[]

	int oddeven=0;

	for(i=1; i<=max_keynum; i++){	//loop through all valid sections
		if(num_outsections[i]>0) {
			if(M==2 && L==2){
				if(testvector[i] == 0 && testvector2[i] == 0){	// check if testvector[] and testvector2[] are both 0 (if they're zero, then this section doesnt connect to testsection (in even or odd number of steps.
		//			printf("%d,", i);
					numnoconnect++;
					nonconnectingsections[indx] = i;	//record that non-connecting section number in nonconnectingsections[]
					indx++;
				}
				else{
					numconnect++;
					connectingsections[indx2] = i;	//record that connecting section number in connectingsections[]
					indx2++;
					if(testvector[i] != 0 && testvector2[i] != 0){
		//				printf("TEST SECTION %d CAN GET TO SECTION %d IN BOTH AN EVEN AND ODD NUMBER OF STEPS!!!!!!!\n", testsection, i);
						oddeven++;
					}
				}
			}
			else{
				if(testvector[i] == 0){	// check if testvector[] is 0 (if it's zero, then this section doesnt connect to testsection
		//			printf("%d,", i);
					numnoconnect++;
					nonconnectingsections[indx] = i;	//record that non-connecting section number in nonconnectingsections[]
					indx++;
				}
				else{
					numconnect++;
					connectingsections[indx2] = i;	//record that connecting section number in connectingsections[]
					indx2++;
				}
			}
		}
	}
//	printf("oddeven=%d\n", oddeven);
	printf("\n\n");
	printf("ACTUAL 'VALID' SECTIONS FOR THE L=%d, M=%d COMPRESSED CASE IS: %d-%lu=%lu=%lu\n", L, M, max_keynum, numnoconnect, max_keynum-numnoconnect,numconnect);
	
	//loop through all 2-spans. if either (should just be able to check one) of the sections in the 2-span are in nonconnectingsections[], subtract one from the "total valid 2-spans"
	unsigned long int actual2spantotal=0;
	i=0;
	while(connectingsections[i] != 0){
		actual2spantotal = actual2spantotal + num_outsections[connectingsections[i]];
		i++;
	}

	printf("ACTUAL 'VALID' 2SPANS FOR THE L=%d, M=%d COMPRESSED CASE IS: %lu\n", L, M, actual2spantotal);
	printf("\n\n");






	//take sections that aren't connected "out" of the matrix
	//need to make the appropriate changes to num_outsections, tspans_edges, tspans_outsection.

	printf("must remove %lu invalid sections\n", numnoconnect);
	for(k=0; k<=numnoconnect-1; k++){	//loop through all sections that don't connect
//		printf("Taking out section %lu, k=%d of %lu\n", nonconnectingsections[k], k, numnoconnect-1);
		//TSPAN REMOVAL
		for(j=1; j<=num_outsections[nonconnectingsections[k]]; j++){
			t_outsection[nonconnectingsections[k]][j] = 0;					//CHANGED OUTSECTION
			for(i=0; i<=t_num_walks[nonconnectingsections[k]][j]-1; i++){
				for(l=0; l<=2; l++){
					t_start[nonconnectingsections[k]][j][l][i]=-1;
					t_end[nonconnectingsections[k]][j][l][i]=-1;
				}
				for(l=0; l<=vM*vL; l++){
					t_walks[nonconnectingsections[k]][j][l][i]=0;
				}
			}
			t_num_walks[nonconnectingsections[k]][j] = 0;						//CHANGED NUM_WALKS

			for(i=0; i<=t_num_walks2[nonconnectingsections[k]][j]-1; i++){
				for(l=0; l<=2; l++){
					t_start2[nonconnectingsections[k]][j][l][i]=-1;
					t_end2[nonconnectingsections[k]][j][l][i]=-1;
				}
				for(l=0; l<=vM*vL; l++){
					t_walks2[nonconnectingsections[k]][j][l][i]=0;
				}
			}
			t_num_walks2[nonconnectingsections[k]][j] = 0;						//CHANGED NUM_WALKS2
		}
		num_outsections[nonconnectingsections[k]]=0;                            //CHANGED NUM_OUT

		//ENDHINGE REMOVAL
		for(j=1; j<=num_left_endhinges[nonconnectingsections[k]]; j++){
			for(i=0; i<=Lend_num_walks[nonconnectingsections[k]][j]-1; i++){
				for(l=0; l<=2; l++){
					Lend_start[nonconnectingsections[k]][j][l][i]=-1;
					Lend_end[nonconnectingsections[k]][j][l][i]=-1;
				}
				for(l=0; l<=vM*vL; l++){
					Lend_walks[nonconnectingsections[k]][j][l][i]=0;
				}
			}
			Lend_num_walks[nonconnectingsections[k]][j] = 0;						//CHANGED NUM_WALKS

			for(i=0; i<=Lend_num_walks2[nonconnectingsections[k]][j]-1; i++){
				for(l=0; l<=2; l++){
					Lend_start2[nonconnectingsections[k]][j][l][i]=-1;
					Lend_end2[nonconnectingsections[k]][j][l][i]=-1;
				}
				for(l=0; l<=vM*vL; l++){
					Lend_walks2[nonconnectingsections[k]][j][l][i]=0;
				}
			}
			Lend_num_walks2[nonconnectingsections[k]][j] = 0;						//CHANGED NUM_WALKS2
		}
		num_left_endhinges[nonconnectingsections[k]]=0;                            //CHANGED NUM_OUT

		for(j=1; j<=num_right_endhinges[nonconnectingsections[k]]; j++){
			for(i=0; i<Rend_num_walks[nonconnectingsections[k]][j]-1; i++){
				for(l=0; l<=2; l++){
					Rend_start[nonconnectingsections[k]][j][l][i]=-1;
					Rend_end[nonconnectingsections[k]][j][l][i]=-1;
				}
				for(l=0; l<=vM*vL; l++){
					Rend_walks[nonconnectingsections[k]][j][l][i]=0;
				}
			}
			Rend_num_walks[nonconnectingsections[k]][j] = 0;						//CHANGED NUM_WALKS

			for(i=0; i<Rend_num_walks2[nonconnectingsections[k]][j]-1; i++){
				for(l=0; l<=2; l++){
					Rend_start2[nonconnectingsections[k]][j][l][i]=-1;
					Rend_end2[nonconnectingsections[k]][j][l][i]=-1;
				}
				for(l=0; l<=vM*vL; l++){
					Rend_walks2[nonconnectingsections[k]][j][l][i]=0;
				}
			}
			Rend_num_walks2[nonconnectingsections[k]][j] = 0;						//CHANGED NUM_WALKS2
		}
		num_right_endhinges[nonconnectingsections[k]]=0;                            //CHANGED NUM_OUT
	}

/*
	char filename3[100];
	sprintf(filename3, "TMinfoL%dM%d_MC.txt", L, M);
	FILE* fp3 = fopen(filename3, "w");
	for(i=1; i<=max_keynum; i++){
		for(j=1; j<=num_outsections[i]; j++){
			fprintf(fp3, "%d %d %lu, (%lu,%lu)->(%lu,%lu) OG:(%lu,%lu)->(%lu,%lu)\n", i, j, t_outsection[i][j], sectionkey2SAP[i][0], sectionkey2SAP[i][1], sectionkey2SAP[t_outsection[i][j]][0], sectionkey2SAP[t_outsection[i][j]][1], sectionkey[sectionkey2SAP[i][0]], sectionkey[sectionkey2SAP[i][1]], sectionkey[sectionkey2SAP[t_outsection[i][j]][0]], sectionkey[sectionkey2SAP[t_outsection[i][j]][1]]);
		}
	}
*/


	//set t_nrr properly now
	for(i=1; i<=max_keynum; i++){
		t_nrr[i] = unsgnlong_vecalloc(1, num_outsections[i]);
	}

	unsigned long int ii, jj;
	unsigned long int tspan_num=0;
	for(ii=1; ii<=max_keynum; ii++){
		for (jj=1; jj<=num_outsections[ii];jj++){
			t_nrr[ii][jj] = ++tspan_num;
//			printf("ii=%lu, jj=%lu, entry=%lu\n", ii, jj, t_nrr[ii][jj]);
		}
	}
	


	printf("Finished taking out invalid sections\n");

	printf("CHECK:\n");
	unsigned long int seccheck=0;
	for(i=1; i<=max_keynum; i++){
		seccheck += num_outsections[i];
	}
	printf("Total tspans check = %lu\n", seccheck);



	printf("Finished generating all 2-spans and endhinges.\n");

	if(totalspan<2){
		printf("totalspan=%d. This program only works for span>=2. Exitting\n", totalspan);
		exit(1);
	}

	printf("\nNOW SAMPLING 2SAPs: %d samples from L=%d, M=%d, span=%d\n", samplesize, L, M, totalspan);


	sprintf(filename, "MC2SAPsHamL%dM%dspan%drun%dnum%lu.txt", L, M, totalspan, runnum, filenum);
	fp = fopen(filename, "w");	//create or overwrite the file "filename"

	if(fp != NULL){
		fprintf(fp, "UofS\n");	//first line in file is always "UofS"
		printf("printed UofS in file '%s'\n", filename);
	}
	else{
		printf("file pointer is pointing to NULL\n");
		exit(1);
	}



	unsigned long int secnum, sec1, sec2;
	int nth_endhinge, nth_walk;
	int curspan;
	int nth_tspan;

	int chosenLEH;	//chosen left endhinge
	int curLEH;	//current left endhinge
	int chosenREH;

	unsigned long int curSample=0;

	double prob;
	double sumofprobs;
	unsigned long int cur_tspan_num;

	//Create and fill left and right eigenvectors
//	double* L_Evector;
//	L_Evector = (double*)malloc(sizeof(double)*(max_tspans+1));
//	if(L_Evector==NULL){
//		fprintf(stderr, "Out of memory");
//		exit(0);
//	}

//	static double R_Evector[max_tspans+1];

	double* R_Evector;
	R_Evector = (double*)malloc(sizeof(double)*(max_tspans+1));
	if(R_Evector==NULL){
		fprintf(stderr, "Out of memory");
		exit(0);
	}


	char filename2[100];		//this will hold the filename of the text file that we will read from
	sprintf(filename2, "2SAP_R_EvectorHam_TS_L%dM%d.txt", L, M);
	FILE* fp2 = fopen(filename2, "r");
	for(i=1; i<= max_tspans; i++){
		fscanf(fp2, "%lf\n", &R_Evector[i]);
	}

/*
	sumofprobs=0.0;
	//check transition probabilities
	for(i=1; i<=max_keynum; i++){
		for(j=1; j<=num_outsections[i]; j++){	//tspan ij
			if(R_Evector[t_nrr[i][j]]>0.00000001){
				printf("\ntspan %lu connects sections %d->%lu\n", t_nrr[i][j], i, t_outsection[i][j]);
				printf("tspan %lu connects to %lu other tspans\n", t_nrr[i][j], num_outsections[t_outsection[i][j]]);
				sumofprobs=0.0;
				for(k=1; k<=num_outsections[t_outsection[i][j]]; k++){
					printf("connection to tspan %lu. Adding (1/%f)(R_Evector[%lu]/R_Evector[%lu]).\n", t_nrr[t_outsection[i][j]][k], dom_evalue, t_nrr[t_outsection[i][j]][k], t_nrr[i][j]);
					printf("= (1/%f)(%f/%f) = %f\n", dom_evalue, R_Evector[t_nrr[t_outsection[i][j]][k]], R_Evector[t_nrr[i][j]], 1/dom_evalue/R_Evector[t_nrr[i][j]]*R_Evector[t_nrr[t_outsection[i][j]][k]]);
					sumofprobs += 1/dom_evalue/R_Evector[t_nrr[i][j]]*R_Evector[t_nrr[t_outsection[i][j]][k]];
					printf("sumofprobs is now %f\n", sumofprobs);
				}
				if(sumofprobs>1.1 || sumofprobs <0.9){
					printf("transition probs don't add up to one.\n");
					exit(1);
				}
			}
			else{
				printf("skipping tspan %lu since it's invalid\n", t_nrr[i][j]);
			}
		}
	}

exit(1);
*/




/*	for(i=1; i<= max_tspans; i++){
//		printf("L_Evector[%d] = %.15f\n", i, L_Evector[i]);
		printf("R_Evector[%d] = %.15f\n", i, R_Evector[i]);
	}
*/


	double maxt_one=0.0;	//used in rejection
	double t_one[tot_left_endhinges];	//each left end hinge has a t_one value.
	curLEH=0;
	for(i=1; i<=max_keynum; i++){
		for(j=1; j<=num_left_endhinges[i]; j++){
			t_one[curLEH]=0;
			for(k=1; k<=num_outsections[i]; k++){
				t_one[curLEH] = t_one[curLEH] + R_Evector[t_nrr[i][k]];
			}
			if(t_one[curLEH]>maxt_one){
				maxt_one = t_one[curLEH];
			}
//			printf("filled t_one[%d]=%f. maxt_one=%f\n", curLEH, t_one[curLEH], maxt_one);
			curLEH++;
		}
	}
	unsigned long int reject_one=0;


	double maxt_two=0.0;	//used in rejection
//	static double t_two[max_tspans+1];

	double* t_two;
	t_two = (double*)malloc(sizeof(double)*(max_tspans+1));
	if(t_two==NULL){
		fprintf(stderr, "Out of memory");
		exit(0);
	}


	for(i=1; i<=max_keynum; i++){
		for(j=1; j<=num_outsections[i]; j++){
//			printf("keynum=%d, outsection_ent=%d, num_right_endhinges[%lu]=%lu, R_Evector[%lu]=%f\n", i, j, t_outsection[i][j], num_right_endhinges[t_outsection[i][j]], t_nrr[i][j], R_Evector[t_nrr[i][j]]);
			if(R_Evector[t_nrr[i][j]]>0.000000000000001){
				t_two[t_nrr[i][j]] = num_right_endhinges[t_outsection[i][j]] / R_Evector[t_nrr[i][j]];
			}
			else{
				t_two[t_nrr[i][j]]=0;
			}
			if(t_two[t_nrr[i][j]] > maxt_two){
				maxt_two = t_two[t_nrr[i][j]];
			}
//			printf("filled t_two[%lu]=%f. maxt_two=%f\n", t_nrr[i][j], t_two[t_nrr[i][j]], maxt_two);
		}
	}
	unsigned long int reject_two=0;




//	double chance=0.0;
	while(curSample<=samplesize-1){

		//initialize built stuff to zeros or -1
		for(i=0; i<=vM*vL/2-1; i++){
			for(j=0; j<=2; j++){
				built_walks_start[i][j]=-1;
				built_walks_end[i][j]=-1;
				built_walks_start2[i][j]=-1;
				built_walks_end2[i][j]=-1;
			}
			for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
				built_walks_direcs[i][j]=0;
				built_walks_direcs2[i][j]=0;
			}
		}

		chosenLEH = floor(ran1real_()*tot_left_endhinges); //random number [0,tot_left_endhinges-1]

//	printf("\ntot_left_endhinges=%lu\n", tot_left_endhinges);
//	chance = (double)1/tot_left_endhinges;
//	printf("chose left endhinge.\n");
//	printf("numerator=%d. denominator=%d\n", 1, tot_left_endhinges);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);

//		printf("\nSample=%lu. chosenLEH=%d, from range [%d,%lu].\n", curSample, chosenLEH, 0, tot_left_endhinges-1);
		curLEH=0;
		for(secnum=1; secnum<=max_keynum; secnum++){
			for(nth_endhinge=1; nth_endhinge<=num_left_endhinges[secnum]; nth_endhinge++){
//				printf("curLEH=%d. secnum=%lu, nth_endhinge=%d\n", curLEH, secnum, nth_endhinge);
				if(curLEH==chosenLEH){
//					printf("chosenLEH=%d was found. Section=%lu, nth_endhinge=%d\n", chosenLEH, secnum, nth_endhinge);

					//CHECK FOR REJECTION
//					printf("accept prob_one=%f\n", t_one[curLEH] / maxt_one);
					if( ran1real_() < t_one[curLEH] / maxt_one ){	//do not reject
//	chance = chance * t_one[curLEH] / maxt_one;
//	printf("passed first rejection with prob=%f.\n", t_one[curLEH] / maxt_one);
//	printf("numerator=%f. denominator=%f\n", t_one[curLEH], maxt_one);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);

						//fill built_walks with left endhinge
						num_built_walks=Lend_num_walks[secnum][nth_endhinge];
						num_built_walks2=Lend_num_walks2[secnum][nth_endhinge];
						for(nth_walk=0; nth_walk<=num_built_walks-1; nth_walk++){
							for(i=0; i<=2; i++){
								built_walks_start[nth_walk][i] = Lend_start[secnum][nth_endhinge][i][nth_walk];
								built_walks_end[nth_walk][i] = Lend_end[secnum][nth_endhinge][i][nth_walk];
							}
							for(i=0; i<=vM*vL; i++){
								built_walks_direcs[nth_walk][i] = Lend_walks[secnum][nth_endhinge][i][nth_walk];
							}
							for(i=vM*vL+1; i<=vM*vL*(totalspan+1)-1; i++){
								built_walks_direcs[nth_walk][i]=0;
							}
						}
						for(nth_walk=0; nth_walk<=num_built_walks2-1; nth_walk++){
							for(i=0; i<=2; i++){
								built_walks_start2[nth_walk][i] = Lend_start2[secnum][nth_endhinge][i][nth_walk];
								built_walks_end2[nth_walk][i] = Lend_end2[secnum][nth_endhinge][i][nth_walk];
							}
							for(i=0; i<=vM*vL; i++){
								built_walks_direcs2[nth_walk][i] = Lend_walks2[secnum][nth_endhinge][i][nth_walk];
							}
							for(i=vM*vL+1; i<=vM*vL*(totalspan+1)-1; i++){
								built_walks_direcs2[nth_walk][i]=0;
							}
						}
						//built_walks is filled with a left endhinge.
						curspan=1;
//						printf("\nLeft endhinge has been added to built_walks:\n");
//						printbuiltwalks_all();
//						printbuiltwalks2_all();
//						printf("\n");

						sec1=secnum;

						//choose first tspan with probability proportional to it's right eigenvector value (from those t-spans who have first section = sec1)
						prob = ran1real_();
//						printf("prob=%f\n", prob);
						nth_tspan=1;
						sumofprobs = R_Evector[t_nrr[sec1][nth_tspan]] / t_one[curLEH];
						while(sumofprobs<prob){
							nth_tspan++;
							sumofprobs = sumofprobs + R_Evector[t_nrr[sec1][nth_tspan]] / t_one[curLEH];
						}
//	chance = chance * R_Evector[t_nrr[sec1][nth_tspan]] / t_one[curLEH];
//	printf("first block chosen with prob=%f. It was the %d-th choice (out of %lu choices). It consists of sections %lu->%lu (tspan number %lu).\n", R_Evector[t_nrr[sec1][nth_tspan]] / t_one[curLEH], nth_tspan, num_outsections[sec1], sec1, t_outsection[sec1][nth_tspan], t_nrr[sec1][nth_tspan]);
//	printf("numerator=%f. denominator=%f\n", R_Evector[t_nrr[sec1][nth_tspan]], t_one[curLEH]);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);

						add_to_built_walks(sec1, nth_tspan);	//add nth_tspan from sec1 to built_walks

//							printf("first block has been added, here is the result:\n");
//							printbuiltwalks_all();
//							printbuiltwalks2_all();
//							printf("\n");


						sec2 = t_outsection[sec1][nth_tspan];	//set sec2
						cur_tspan_num = t_nrr[sec1][nth_tspan];
						if(R_Evector[cur_tspan_num] > 0.000000000000001){
							curspan++;

							while(curspan<totalspan){
//								printf("adding another tspan to section %lu. Probabilities:\n", sec2);
								sumofprobs=0.0;
								for(i=1; i<=num_outsections[sec2]; i++){
									sumofprobs = sumofprobs + R_Evector[t_nrr[sec2][i]] / R_Evector[cur_tspan_num] / dom_evalue;
//										printf("sec %lu -> sec %lu. tspan_num=%lu. R_Evector[%lu]=%f. probability=%f. cumulative prob=%f\n", sec2, t_outsection[sec2][i], t_nrr[sec2][i], t_nrr[sec2][i], R_Evector[t_nrr[sec2][i]], R_Evector[t_nrr[sec2][i]] / R_Evector[cur_tspan_num] / dom_evalue, sumofprobs);
								}
//								printf("sumofprobs=%f\n", sumofprobs);
								if(sumofprobs<0.9 || sumofprobs>1.1){
									printf("second sumofprobs=%f\n", sumofprobs);
									printf("problem with second sumofprobs when leaving tspan number %lu, (made up of %lu and %lu)\n", cur_tspan_num, sec1, sec2);
									printf("This was on sample %lu. Exitting\n", curSample);
									exit(1);
								}
	

								//choose next tspan according to correct probability
								prob = ran1real_();
								prob = prob * sumofprobs; //in case sum of probs doesn't add up to one.
	//							printf("prob=%f\n", prob);
								nth_tspan=1;
								sumofprobs = R_Evector[t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue;
								while(sumofprobs<prob){
									nth_tspan++;
									sumofprobs = sumofprobs + R_Evector[t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue;
								}
	//	chance = chance * R_Evector[t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue;
//		printf("next block chosen with prob=%f. It is sections %lu->%lu (tspan number %lu)\n", R_Evector[t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue, sec2, t_outsection[sec2][nth_tspan], t_nrr[sec2][nth_tspan]);
	//	printf("numerator=%f. denominator=%f. dom_evalue=%f.\n", R_Evector[t_nrr[sec2][nth_tspan]], R_Evector[cur_tspan_num], dom_evalue);
	//	printf("chance = %f\n", chance);
	//	printf("1/chance = %f\n\n", 1/chance);
	//							printf("sec %lu -> sec %lu. tspan_num=%lu. probability=%f. cumulative prob=%f\n", sec2, t_outsection[sec2][nth_tspan], t_nrr[sec2][nth_tspan], R_Evector[t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue, sumofprobs);
								sec1=sec2;
								add_to_built_walks(sec1, nth_tspan);	//add nth_tspan from sec2 to built_walks

//								printf("tspan has added, here is the result:\n");
//								printbuiltwalks_all();
//								printbuiltwalks2_all();
//								printf("\n");
								sec2 = t_outsection[sec1][nth_tspan];
								cur_tspan_num = t_nrr[sec1][nth_tspan];
								if(R_Evector[cur_tspan_num] < 0.000000000000001){
									printf("tspan chosen (during intermediate tspans) that will not result in a poly. Shouldn't ever happen. Exitting\n");
									printf("cur_tspan_num=%lu, sections %lu->%lu (%d-th outsec).\n", cur_tspan_num, sec1, t_outsection[sec1][nth_tspan], nth_tspan);
									exit(1);
								}
								curspan++;
							}
							//CHECK FOR REJECTION
//									printf("accept prob_two=%f\n", t_two[cur_tspan_num] / maxt_two);
							if( ran1real_() < t_two[cur_tspan_num] / maxt_two ){	//do not reject
	//	chance = chance * t_two[cur_tspan_num] / maxt_two;
//				printf("passed 2nd rejection with prob=%f.\n", t_two[cur_tspan_num] / maxt_two);
	//	printf("numerator=%f. denominator=%f\n", t_two[cur_tspan_num], maxt_two);
	//	printf("chance = %f\n", chance);
	//	printf("1/chance = %f\n\n", 1/chance);
							chosenREH = 1 + floor(ran1real_()*num_right_endhinges[sec2]); //rando number [1, num_right_endhinges[sec2]]
	//	chance = chance * (double)1/num_right_endhinges[sec2];
//				printf("chose final endhinge with prob=%f.\n", (double)1/num_right_endhinges[sec2]);
	//	printf("chance = %f\n", chance);
	//	printf("1/chance = %f\n\n", 1/chance);
//										printf("closing off SAP with endhinge: section=%lu, nth_endhinge=%d (possible [%d,%lu]).\n", sec2, chosenREH, 1, num_right_endhinges[sec2]);
								add_right_endhinge(sec2, chosenREH);

		//						printf("FINAL POLY:\n");
		//						printbuiltwalks();

								//sample created
	//	printf("FINAL CHANCE OF CHOOSING:\n");
	//	printf("chance = %f\n", chance);
	//	printf("1/chance = %f\n", 1/chance);
								curSample++;

								//PRINT TO FILE.
								printtofile();
							}
							else{	//rejected at closure
								reject_two++;
							}
						}
						else{
							printf("tspan chosen that will not result in a poly. Shouldn't ever happen. Exitting\n");
							exit(1);
						}
					}
					else{	//rejected at first step
						reject_one++;
					}
				}
				else if(curLEH>chosenLEH){
					break;
				}
				curLEH++;
			}
		}
	}


	fprintf(fp, "-999\n");
	fclose(fp);
	printf("File %lu contains %lu polys.\n", filenum, filetotal);

	printf("\nSampling Complete.\n");

	printf("\nRESULTS: %lu samples created from L=%d, M=%d, span=%d\n", curSample, L, M, totalspan);
	printf("Rejected at first step %lu times\n", reject_one);
	printf("Rejected at last step %lu times\n", reject_two);


	printf("\n");	



	system("date");  /* prints the date and time */
	printf("\nProgram Complete.\n\n");

	return 0;

}				/* end of main */
/***************************************************************************/


/***** start of function conv_to_array *************************************/
/* called only by main */


void conv_to_array(void)
{
	//Fill sectionkey2SAP[key][2]
		//sectionkey2SAP[key][0] = secnum of poly1
		//sectionkey2SAP[key][1] = secnum of poly2
		//key=key number assigned to this ordered pairing.
	//Also fill num_outsections at the same time.
	int i, j;
	int found;

	for(i=1; i<=max_keynum; i++){
		sectionkey2SAP[i][0]=0;
		sectionkey2SAP[i][1]=0;
		num_outsections[i]=0;
	}

	int keycounter=0;

	for(i=1; i<=max_sections; i++){
		current_hinge_span[i] = first_hinge_span[i];
		while(current_hinge_span[i]->nexthinge != NULL) {
			current_hinge_span[i] = current_hinge_span[i]->nexthinge;
//			printf("looking at section composed of sections %d and %lu (on the left)\n", i, current_hinge_span[i]->inorder2);
			j=1;
			found=0;
			while(sectionkey2SAP[j][0]!=0){
				if(sectionkey2SAP[j][0]==i && sectionkey2SAP[j][1]==current_hinge_span[i]->inorder2){
//					printf("combo found in key already\n");
					found=1;
					num_outsections[j]++;
					break;
				}
				j++;
			}
			if(found==0){	//need to record sections into sectionkey2SAP
//				printf("combo not in key yet. adding: j=%d, sec1=%d, sec2=%lu\n", j, i, current_hinge_span[i]->inorder2);
				keycounter++;
				sectionkey2SAP[j][0] = i;
				sectionkey2SAP[j][1] = current_hinge_span[i]->inorder2;
				num_outsections[j]++;
			}	
		}
	}

	printf("Total of %d 2SAP sections\n", keycounter);

/*	for(i=1; i<=max_keynum; i++){
		printf("keynum=%d, sections: %lu,%lu (OG: %lu,%lu)\n", i, sectionkey2SAP[i][0], sectionkey2SAP[i][1], sectionkey[sectionkey2SAP[i][0]], sectionkey[sectionkey2SAP[i][1]]);
		printsection(sectionkey[sectionkey2SAP[i][0]]);
		printsection(sectionkey[sectionkey2SAP[i][1]]);
		printf("\n");
	}
*/

	for(i=1; i<=max_keynum; i++){
		t_outsection[i] = unsgnlong_vecalloc(1, num_outsections[i]); //allocate memory
		t_num_walks[i] = unsgn_vecalloc(1, num_outsections[i]);
		t_num_walks2[i] = unsgn_vecalloc(1, num_outsections[i]);
	}



	unsigned long int section_num;
	unsigned long int sec1key;
	unsigned long int sec2key;
	unsigned int curArrayEnt[max_keynum+1];
	for(i=1; i<=max_keynum; i++){
		curArrayEnt[i]=1;
	}

	struct hinge_span *hinge_to_free;


	//fill t_outsection and t_num_walks and t_num_walks2
	for(section_num = 1; section_num <= max_sections; section_num++) {
		current_hinge_span[section_num] = first_hinge_span[section_num];
		while ((*current_hinge_span[section_num]).nexthinge != NULL) {
			current_hinge_span[section_num] = (*current_hinge_span[section_num]).nexthinge;

//			printf("this hinge is made up of sections %lu,%lu -> %lu,%lu\n", section_num, (*current_hinge_span[section_num]).inorder2, (*current_hinge_span[section_num]).outorder, (*current_hinge_span[section_num]).outorder2);

			//get keynums
			sec1key=0;
			for(i=1; i<=max_keynum; i++){
//				printf("keynum=%d, match: %lu,%lu ?\n", i, sectionkey2SAP[i][0], sectionkey2SAP[i][1]);
				if(sectionkey2SAP[i][0]==section_num && sectionkey2SAP[i][1]==(*current_hinge_span[section_num]).inorder2){
//					printf("sec1key=%d\n", i);
					sec1key=i;
					break;
				}
			}
			if(sec1key==0){
				printf("error: didn't find sec1key. Exitting\n");
				exit(1);
			}
			sec2key=0;
			for(i=1; i<=max_keynum; i++){
				if(sectionkey2SAP[i][0]==(*current_hinge_span[section_num]).outorder && sectionkey2SAP[i][1]==(*current_hinge_span[section_num]).outorder2){
//					printf("sec2key=%d\n", i);
					sec2key=i;
					break;
				}
			}
			if(sec2key==0){
				printf("error: didn't find sec2key. Sections were %lu,%lu. Exitting\n", (*current_hinge_span[section_num]).outorder, (*current_hinge_span[section_num]).outorder2);
				exit(1);
			}
//			printf("sec1key=%lu, sec2key=%lu\n", sec1key, sec2key);
			t_outsection[sec1key][curArrayEnt[sec1key]] = sec2key;
			t_num_walks[sec1key][curArrayEnt[sec1key]] = (*current_hinge_span[section_num]).num_of_walks;
			t_num_walks2[sec1key][curArrayEnt[sec1key]] = (*current_hinge_span[section_num]).num_of_walks2;
			curArrayEnt[sec1key]++;
		}
	}



	//allocate space for t_start, t_end, and t_walks -- AND -- t_start2, t_end2, and t_walks2
	printf("allocating space for t_start, t_end, and t_walks AND t_start2, t_end2, and t_walks2\n");

	for(section_num=1; section_num<=max_keynum; section_num++){
		t_start[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_start[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
		t_start2[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_start2[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		t_end[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_end[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
		t_end2[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_end2[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		t_walks[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_walks[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
		t_walks2[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_walks2[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		for (i = 1; i <= num_outsections[section_num]; i++){
			t_start[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (t_start[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			t_start2[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (t_start2[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			t_end[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (t_end[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			t_end2[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (t_end2[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			t_walks[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (t_walks[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			t_walks2[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (t_walks2[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			for (j = 0; j < 3; j++){
				t_start[section_num][i][j] = (int*)malloc(t_num_walks[section_num][i] * sizeof(int));
		   		if (t_start[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				t_start2[section_num][i][j] = (int*)malloc(t_num_walks2[section_num][i] * sizeof(int));
		   		if (t_start2[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}

				t_end[section_num][i][j] = (int*)malloc(t_num_walks[section_num][i] * sizeof(int));
		   		if (t_end[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				t_end2[section_num][i][j] = (int*)malloc(t_num_walks2[section_num][i] * sizeof(int));
		   		if (t_end2[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
			}
			for (j = 0; j <= vM*vL; j++){
				t_walks[section_num][i][j] = (int*)malloc(t_num_walks[section_num][i] * sizeof(int));
		   		if (t_walks[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				t_walks2[section_num][i][j] = (int*)malloc(t_num_walks2[section_num][i] * sizeof(int));
		   		if (t_walks2[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
			}
		}
	}
	printf("Done allocating space for t_start, t_end, and t_walks AND t_start2, t_end2, and t_walks2\n");



	///Fill t_start, t_end, t_walks -- AND -- t_start2, t_end2, t_walks2
	printf("Filling t_start, t_end, and t_walks AND t_start2, t_end2, and t_walks2\n");
	for(i=1; i<=max_keynum; i++){
		curArrayEnt[i]=1;
	}
	for (section_num = 1; section_num <= max_sections; section_num++) {
//		printf("section_num=%lu\n", section_num);
		current_hinge_span[section_num] = first_hinge_span[section_num];
		while ((*current_hinge_span[section_num]).nexthinge != NULL) {
			hinge_to_free = current_hinge_span[section_num];
			current_hinge_span[section_num] = (*current_hinge_span[section_num]).nexthinge;
			free(hinge_to_free);

			//get keynums
			sec1key=0;
			for(i=1; i<=max_keynum; i++){
//				printf("keynum=%d, match: %lu,%lu ?\n", i, sectionkey2SAP[i][0], sectionkey2SAP[i][1]);
				if(sectionkey2SAP[i][0]==section_num && sectionkey2SAP[i][1]==(*current_hinge_span[section_num]).inorder2){
//					printf("sec1key=%d\n", i);
					sec1key=i;
					break;
				}
			}
			if(sec1key==0){
				printf("error: didn't find sec1key. Exitting\n");
				exit(1);
			}
			sec2key=0;
			for(i=1; i<=max_keynum; i++){
				if(sectionkey2SAP[i][0]==(*current_hinge_span[section_num]).outorder && sectionkey2SAP[i][1]==(*current_hinge_span[section_num]).outorder2){
//					printf("sec2key=%d\n", i);
					sec2key=i;
					break;
				}
			}
			if(sec2key==0){
				printf("error: didn't find sec2key. Sections were %lu,%lu. Exitting\n", (*current_hinge_span[section_num]).outorder, (*current_hinge_span[section_num]).outorder2);
				exit(1);
			}
//			printf("\nsec1key=%lu, sec2key=%lu\n", sec1key, sec2key);



			for(i=0; i<=(current_hinge_span[section_num]->num_of_walks)-1; i++){
				for(j=0; j<=2; j++){
//					printf("recording %d-th coordinate of %d-th walk\n", j, i);
//					printf("keynum=%lu, arr_ent=%lu, j=%d, i=%d\n", sec1key, curArrayEnt[sec1key], j, i);
					t_start[sec1key][curArrayEnt[sec1key]][j][i] = current_hinge_span[section_num]->start[j][i];
//					printf("t_start[%lu][%lu][%d][%d] = %d\n", sec1key, curArrayEnt[sec1key], j, i, t_start[sec1key][curArrayEnt[sec1key]][j][i]);
					t_end[sec1key][curArrayEnt[sec1key]][j][i] = current_hinge_span[section_num]->end[j][i];
//					printf("t_end[%lu][%lu][%d][%d] = %d\n", sec1key, curArrayEnt[sec1key], j, i, t_end[sec1key][curArrayEnt[sec1key]][j][i]);
				}
				for(j=0; j<=vM*vL; j++){
					t_walks[sec1key][curArrayEnt[sec1key]][j][i] = current_hinge_span[section_num]->walks[j][i];
				}
			}

			for(i=0; i<=(current_hinge_span[section_num]->num_of_walks2)-1; i++){
				for(j=0; j<=2; j++){
//					printf("recording %d-th coordinate\n", j);
//					printf("keynum=%lu, arr_ent=%lu, j=%d, i=%d\n", sec1key, curArrayEnt[sec1key], j, i);
					t_start2[sec1key][curArrayEnt[sec1key]][j][i] = current_hinge_span[section_num]->start2[j][i];
//					printf("t_start2[%lu][%lu][%d][%d] = %d\n", sec1key, curArrayEnt[sec1key], j, i, t_start2[sec1key][curArrayEnt[sec1key]][j][i]);
					t_end2[sec1key][curArrayEnt[sec1key]][j][i] = current_hinge_span[section_num]->end2[j][i];
//					printf("t_end2[%lu][%lu][%d][%d] = %d\n", sec1key, curArrayEnt[sec1key], j, i, t_end2[sec1key][curArrayEnt[sec1key]][j][i]);
				}
				for(j=0; j<=vM*vL; j++){
					t_walks2[sec1key][curArrayEnt[sec1key]][j][i] = current_hinge_span[section_num]->walks2[j][i];
				}
			}
			curArrayEnt[sec1key]++;

		}
	}
	printf("Done filling t_start, t_end, and t_walks\n");


/*
	int counter=1;
	for (section_num=1; section_num<=max_keynum;section_num++){
		for (nth_outsection=1;nth_outsection<=num_outsections[section_num];nth_outsection++){
			printf("This tspan (%d=%lu) consists of keys %lu->%lu, which are pairs of section keys (%lu,%lu)->(%lu,%lu), which have OG section_nums (%lu, %lu)->(%lu, %lu)\n", counter, tspans_nrr[section_num][nth_outsection], section_num, tspans_outsection[section_num][nth_outsection], sectionkey2SAP[0][section_num], sectionkey2SAP[1][section_num], sectionkey2SAP[0][tspans_outsection[section_num][nth_outsection]], sectionkey2SAP[1][tspans_outsection[section_num][nth_outsection]], sectionkey[sectionkey2SAP[0][section_num]], sectionkey[sectionkey2SAP[1][section_num]], sectionkey[sectionkey2SAP[0][tspans_outsection[section_num][nth_outsection]]], sectionkey[sectionkey2SAP[1][tspans_outsection[section_num][nth_outsection]]]);
			counter++;
		}
	}
	exit(1);
*/

	





	return;
}




/***** start of function conv_to_array *************************************/
/* called only by main */


void conv_endhinges_to_array(void)
{



	int i, j;

	unsigned long int section_num;
	unsigned long int sec1key;
	struct endhinge *hinge_to_free;

	for(i=1; i<=max_keynum; i++){
		num_left_endhinges[i]=0;
		num_right_endhinges[i]=0;
	}

	for (section_num = 1; section_num <= max_sections; section_num++) {
		currentendhinge[section_num] = firstendhinge[section_num];
		while ((*currentendhinge[section_num]).nextendhinge != NULL) {
			currentendhinge[section_num] = (*currentendhinge[section_num]).nextendhinge;

			//get keynum
			sec1key=0;
			for(i=1; i<=max_keynum; i++){
//				printf("keynum=%d, match: %lu,%lu ?\n", i, sectionkey2SAP[i][0], sectionkey2SAP[i][1]);
				if(sectionkey2SAP[i][0]==section_num && sectionkey2SAP[i][1]==(*currentendhinge[section_num]).sec2){
//					printf("sec1key=%d\n", i);
					sec1key=i;
					break;
				}
			}
			if(sec1key==0){
				printf("error: didn't find sec1key. Exitting\n");
				exit(1);
			}

			if(currentendhinge[section_num]->side==0){
				num_left_endhinges[sec1key]++;
			}
			else{
				num_right_endhinges[sec1key]++;
			}
		}
	}

	for(section_num=1; section_num<=max_keynum; section_num++){
//		printf("num_left_endhinges[%lu]=%lu\n", section_num, num_left_endhinges[section_num]);
//		printf("num_right_endhinges[%lu]=%lu\n", section_num, num_right_endhinges[section_num]);
		Lend_num_walks[section_num] = unsgn_vecalloc(1, num_left_endhinges[section_num]);
		Rend_num_walks[section_num] = unsgn_vecalloc(1, num_right_endhinges[section_num]);
		Lend_num_walks2[section_num] = unsgn_vecalloc(1, num_left_endhinges[section_num]);
		Rend_num_walks2[section_num] = unsgn_vecalloc(1, num_right_endhinges[section_num]);
	}





	//allocate space for Lend_start, Lend_end, and Lend_walks AND Lend_start2, Lend_end2, and Lend_walks2 (and Rend stuff)

	for(section_num=1; section_num<=max_keynum; section_num++){
/////////
		Lend_start[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_start[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
		Lend_start2[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_start2[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		Lend_end[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_end[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
		Lend_end2[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_end2[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		Lend_walks[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_walks[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
		Lend_walks2[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_walks2[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
/////////
		Rend_start[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_start[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
		Rend_start2[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_start2[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		Rend_end[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_end[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
		Rend_end2[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_end2[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		Rend_walks[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_walks[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
		Rend_walks2[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_walks2[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
/////////
		for (i = 1; i <= num_left_endhinges[section_num]; i++){
			Lend_start[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Lend_start[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			Lend_start2[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Lend_start2[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			Lend_end[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Lend_end[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			Lend_end2[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Lend_end2[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			Lend_walks[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (Lend_walks[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			Lend_walks2[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (Lend_walks2[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			for (j = 0; j < 3; j++){
				Lend_start[section_num][i][j] = (int*)malloc(Lend_num_walks[section_num][i] * sizeof(int));
		   		if (Lend_start[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				Lend_start2[section_num][i][j] = (int*)malloc(Lend_num_walks2[section_num][i] * sizeof(int));
		   		if (Lend_start2[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}

				Lend_end[section_num][i][j] = (int*)malloc(Lend_num_walks[section_num][i] * sizeof(int));
		   		if (Lend_end[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				Lend_end2[section_num][i][j] = (int*)malloc(Lend_num_walks2[section_num][i] * sizeof(int));
		   		if (Lend_end2[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
			}
			for (j = 0; j <= vM*vL; j++){
				Lend_walks[section_num][i][j] = (int*)malloc(Lend_num_walks[section_num][i] * sizeof(int));
		   		if (Lend_walks[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
			}
			for (j = 0; j <= vM*vL; j++){
				Lend_walks2[section_num][i][j] = (int*)malloc(Lend_num_walks2[section_num][i] * sizeof(int));
		   		if (Lend_walks2[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
			}
		}
/////////
		for (i = 1; i <= num_right_endhinges[section_num]; i++){
			Rend_start[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Rend_start[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			Rend_start2[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Rend_start2[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			Rend_end[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Rend_end[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			Rend_end2[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Rend_end2[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			Rend_walks[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (Rend_walks[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			Rend_walks2[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (Rend_walks2[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			for (j = 0; j < 3; j++){
				Rend_start[section_num][i][j] = (int*)malloc(Rend_num_walks[section_num][i] * sizeof(int));
		   		if (Rend_start[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				Rend_start2[section_num][i][j] = (int*)malloc(Rend_num_walks2[section_num][i] * sizeof(int));
		   		if (Rend_start2[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}

				Rend_end[section_num][i][j] = (int*)malloc(Rend_num_walks[section_num][i] * sizeof(int));
		   		if (Rend_end[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				Rend_end2[section_num][i][j] = (int*)malloc(Rend_num_walks2[section_num][i] * sizeof(int));
		   		if (Rend_end2[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
			}
			for (j = 0; j <= vM*vL; j++){
				Rend_walks[section_num][i][j] = (int*)malloc(Rend_num_walks[section_num][i] * sizeof(int));
		   		if (Rend_walks[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				Rend_walks2[section_num][i][j] = (int*)malloc(Rend_num_walks2[section_num][i] * sizeof(int));
		   		if (Rend_walks2[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
			}
		}
	}




	///Fill L/Rend_start, L/Rend_end, L/Rend_walks
	unsigned int L_curArrayEnt[max_keynum+1];
	unsigned int R_curArrayEnt[max_keynum+1];
	for(i=1; i<=max_keynum; i++){
		L_curArrayEnt[i]=1;
		R_curArrayEnt[i]=1;
	}
	for (section_num = 1; section_num <= max_sections; section_num++) {
		currentendhinge[section_num] = firstendhinge[section_num];
		while ((*currentendhinge[section_num]).nextendhinge != NULL) {
			hinge_to_free = currentendhinge[section_num];
			currentendhinge[section_num] = (*currentendhinge[section_num]).nextendhinge;
			free(hinge_to_free);

			//get keynum
			sec1key=0;
			for(i=1; i<=max_keynum; i++){
//				printf("keynum=%d, match: %lu,%lu ?\n", i, sectionkey2SAP[i][0], sectionkey2SAP[i][1]);
				if(sectionkey2SAP[i][0]==section_num && sectionkey2SAP[i][1]==(*currentendhinge[section_num]).sec2){
//					printf("sec1key=%d\n", i);
					sec1key=i;
					break;
				}
			}
			if(sec1key==0){
				printf("error: didn't find sec1key. Exitting\n");
				exit(1);
			}

			if(currentendhinge[section_num]->side==0){	//leftendhinge
				Lend_num_walks[sec1key][L_curArrayEnt[sec1key]] = currentendhinge[section_num]->num_of_walks;
				Lend_num_walks2[sec1key][L_curArrayEnt[sec1key]] = currentendhinge[section_num]->num_of_walks2;
				for(i=0; i<=(currentendhinge[section_num]->num_of_walks)-1; i++){
					for(j=0; j<=2; j++){
						Lend_start[sec1key][L_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->start[j][i];
						Lend_end[sec1key][L_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->end[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						Lend_walks[sec1key][L_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->walks[j][i];
					}
				}
				for(i=0; i<=(currentendhinge[section_num]->num_of_walks2)-1; i++){
					for(j=0; j<=2; j++){
						Lend_start2[sec1key][L_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->start2[j][i];
						Lend_end2[sec1key][L_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->end2[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						Lend_walks2[sec1key][L_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->walks2[j][i];
					}
				}
				L_curArrayEnt[sec1key]++;
			}
			else{ //rightendhinge
				Rend_num_walks[sec1key][R_curArrayEnt[sec1key]] = currentendhinge[section_num]->num_of_walks;
				Rend_num_walks2[sec1key][R_curArrayEnt[sec1key]] = currentendhinge[section_num]->num_of_walks2;
				for(i=0; i<=(currentendhinge[section_num]->num_of_walks)-1; i++){
					for(j=0; j<=2; j++){
						Rend_start[sec1key][R_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->start[j][i];
						Rend_end[sec1key][R_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->end[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						Rend_walks[sec1key][R_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->walks[j][i];
					}
				}
				for(i=0; i<=(currentendhinge[section_num]->num_of_walks2)-1; i++){
					for(j=0; j<=2; j++){
						Rend_start2[sec1key][R_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->start2[j][i];
						Rend_end2[sec1key][R_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->end2[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						Rend_walks2[sec1key][R_curArrayEnt[sec1key]][j][i] = currentendhinge[section_num]->walks2[j][i];
					}
				}
				R_curArrayEnt[sec1key]++;
			}
		}
	}
	
	return;
}



/***************************************************************************/

void
enterhinge(int i, int j, int side, int (*pointordNum)[6], int curlength)
{
	if( !(alreadyentered[i][j]==1 && side==0) ){
	//	printf("entered hinge (i=%d, j=%d, side=%d\n", i, j, side);

		curlength++;
		num_walks++;
		if(side==0){
			curstart[0][num_walks-1]=0;
			curwalks[curlength-1][num_walks-1] = 1;
		}
		else{ //side==1
			curstart[0][num_walks-1]=2;
			curwalks[curlength-1][num_walks-1] = 2;
		}
		curstart[1][num_walks-1]=j;
		curstart[2][num_walks-1]=i;


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

		leavehinge(i, j, side ^ 1, pointordNum, curlength);
			/* the SAW leaves the hinge out the opposite side, hence (side ^ 1) is passed */
			/* if side = 1 then side ^ 1 = 0 */
			/* if side = 0 then side ^ 1 = 1 */
		rowedges(i, j, pointordNum, curlength);
			/* horizontal edges are explored */
		coledges(i, j, pointordNum, curlength);
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

		curstart[0][num_walks-1]=-1;
		curstart[1][num_walks-1]=-1;
		curstart[2][num_walks-1]=-1;
		curwalks[curlength-1][num_walks-1] = 0;

		curlength--;
		num_walks--;
	//	printf("finished enterhinge\n");
	}
	return;

}				/* end of function enterhinge */

/***************************************************************************/

void
leavehinge(int i, int j, int side, int (*pointordNum)[6], int curlength)
/* called by: rowedges,coledges,enterhinge */
{
	if( !(alreadyentered[i][j]==1 && side==0) ){
	//	printf("left hinge (i=%d, j=%d, side=%d\n", i, j, side);

		curlength++;
		if(side==0){
			curend[0][num_walks-1]=0;
			curwalks[curlength-1][num_walks-1] = 2;
		}
		else{ //side==1
			curend[0][num_walks-1]=2;
			curwalks[curlength-1][num_walks-1] = 1;
		}
		curend[1][num_walks-1]=j;
		curend[2][num_walks-1]=i;


		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		ordertemplate[side][i][j] = (*pointordNum)[side];	/* record the section
									 * edge (i.e. the
									 * leaving edge) */
		(*pointordNum)[side]++;
			/* the edge number is incremented so that the next edge gets a higher number */

		if (side == 0 && (*pointordNum)[1] > 1) {
	 	/* if the leaving edge is on the left and there is at least one edge on the right section */
			if (*pointordNum[0]-1 <= vM * vL - 2) { //make sure there is at least 2 vertices open
			// so that a second polygon can be considered
				for (ii = 0; ii <= M; ii++) { // these are the "vertical" locations of the vertices being considered
					for (jj = 0; jj <= L; jj++) { // these are the horizontal locations of the vertices being considered
						if (hingestatus[ii][jj] == 0) {
							// if the vertex is not occupied then explore entering the hinge
							// at this vertex with a second polygon
							enterhinge2(ii, jj, side, pointordNum, 0);
						}
						alreadyentered2[ii][jj]=1;
					}
				}
				for (ii = 0; ii <= M; ii++) {
					for (jj = 0; jj <= L; jj++) {
						alreadyentered2[ii][jj]=0;
					}
				}
			}
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
					enterhinge(ii, jj, side, pointordNum,0);
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


		curend[0][num_walks-1]=-1;
		curend[1][num_walks-1]=-1;
		curend[2][num_walks-1]=-1;
		curwalks[curlength-1][num_walks-1] = 0;
		curlength--;
	}
	return;
}
/***************************************************************************/

void
rowedges(int i, int j, int (*pointordNum)[6], int curlength)  //y-direction
{
//	printf("rowedges called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(hingestatus[i][j-1]==0){
			hingestatus[i][j-1]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 4;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			rowhingeedges[i][j-1]=1;
			leavehinge(i, j-1, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(i, j-1, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			rowedges(i, j-1, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			coledges(i, j-1, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i][j-1]=0;
			rowhingeedges[i][j-1]=0;

			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
	if(j<L){
		if(hingestatus[i][j+1]==0){
			hingestatus[i][j+1]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 3;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			rowhingeedges[i][j]=1;
			leavehinge(i, j+1, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(i, j+1, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			rowedges(i, j+1, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			coledges(i, j+1, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i][j+1]=0;
			rowhingeedges[i][j]=0;

			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
//	printf("finished rowedges\n");
	return;
}				/* end of function rowedges */

/***************************************************************************/

void coledges(int i, int j, int (*pointordNum)[6], int curlength)  //z-direction
/* called by rowedges,coledges */
{
//	printf("coledges called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(hingestatus[i-1][j]==0){
			hingestatus[i-1][j]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 6;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			colhingeedges[i-1][j]=1;
			leavehinge(i-1, j, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(i-1, j, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			rowedges(i-1, j, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			coledges(i-1, j, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i-1][j]=0;
			colhingeedges[i-1][j]=0;

			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
	if(i<M){
		if(hingestatus[i+1][j]==0){
			hingestatus[i+1][j]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 5;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			colhingeedges[i][j]=1;
			leavehinge(i+1, j, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(i+1, j, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			rowedges(i+1, j, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			coledges(i+1, j, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i+1][j]=0;
			colhingeedges[i][j]=0;

			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
//	printf("finished coledges\n");
	return;
}				/* end of function coledges */

/***************************************************************************/







/***************************************************************************/

void
enterhinge2(int i, int j, int side, int (*pointordNum)[6], int curlength2)
{
	if( !(alreadyentered2[i][j]==1 && side==0) ){
	//	printf("entered hinge2 (i=%d, j=%d, side=%d\n", i, j, side);

		curlength2++;
		num_walks2++;
		if(side==0){
			curstart2[0][num_walks2-1]=0;
			curwalks2[curlength2-1][num_walks2-1] = 1;
		}
		else{ //side==1
			curstart2[0][num_walks2-1]=2;
			curwalks2[curlength2-1][num_walks2-1] = 2;
		}
		curstart2[1][num_walks2-1]=j;
		curstart2[2][num_walks2-1]=i;


		hingestatus[i][j]=1;
		ordertemplate2[side][i][j] = (*pointordNum)[side+3];
			/* section edge through which the SAW enters the hinge is given a number */
			/* this number represents the order in which it occurs in the section    */

		(*pointordNum)[side+3]++;
			/* the edge number through which the hinge is entered */
			/* is incremented so that the next edge gets a higher number */

		/* the SAW can now do one of three things */
		/* it may leave the hinge in which case leavehinge(...) is called */
		/* it may explore edges in the horizontal direction => rowedges (...) is called */
		/* it may explore edges in the vertical direction => coleges (...) is called    */
		/* in each case the currrent state of the hinge is passed via the parameters    */

		leavehinge2(i, j, side ^ 1, pointordNum, curlength2);
			/* the SAW leaves the hinge out the opposite side, hence (side ^ 1) is passed */
			/* if side = 1 then side ^ 1 = 0 */
			/* if side = 0 then side ^ 1 = 1 */
		rowedges2(i, j, pointordNum, curlength2);
			/* horizontal edges are explored */
		coledges2(i, j, pointordNum, curlength2);
			/* vertical edges are explored */

		(*pointordNum)[side+3]--;
			/* the edge number is decremented since, control is being passed back to  */
			/* the function that called enterhinge. That is to say, the program will  */
			/* now explore alternate ways of entering the hinge (i.e. from different  */
			/* locations). In essense the intial incremention above is nullified so	  */
			/* when enterhinge is called again the entering edge is numbered properly */
		
		ordertemplate2[side][i][j] = 0;
			/* the edge that was recorded must now be erased as this edge may not exist */
			/* as part of the SAW when alternate possibilities are explored. In the case */
			/* that it does not exist it needs to be numbered as "0" (zero).             */
		hingestatus[i][j]=0;

		curstart2[0][num_walks2-1]=-1;
		curstart2[1][num_walks2-1]=-1;
		curstart2[2][num_walks2-1]=-1;
		curwalks2[curlength2-1][num_walks2-1] = 0;

		curlength2--;
		num_walks2--;
	//	printf("finished enterhinge\n");
	}
	return;

}				/* end of function enterhinge */

/***************************************************************************/

void
leavehinge2(int i, int j, int side, int (*pointordNum)[6], int curlength2)
/* called by: rowedges,coledges,enterhinge */
{
	if( !(alreadyentered2[i][j]==1 && side==0) ){
	//	printf("left hinge (i=%d, j=%d, side=%d\n", i, j, side);

		curlength2++;
		if(side==0){
			curend2[0][num_walks2-1]=0;
			curwalks2[curlength2-1][num_walks2-1] = 2;
		}
		else{ //side==1
			curend2[0][num_walks2-1]=2;
			curwalks2[curlength2-1][num_walks2-1] = 1;
		}
		curend2[1][num_walks2-1]=j;
		curend2[2][num_walks2-1]=i;


		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		ordertemplate2[side][i][j] = (*pointordNum)[side+3];	/* record the section
									 * edge (i.e. the
									 * leaving edge) */
		(*pointordNum)[side+3]++;
			/* the edge number is incremented so that the next edge gets a higher number */

		if (side == 0 && (*pointordNum)[4] > 1) {
	 	/* if the leaving edge is on the left and there is at least one edge on the right section */
		/* this may be a valid two span */
			num_2_spans++; /*Keep track of all 2 spans, valid and non-valid, just for interest sake*/
			if (LFlag2(pointordNum)){
			/* if the two span connects to phi on the left and on the right then do the following*/

				//Hamiltonian Check
				int isHam=1;
				for (ii = 0; ii <= M; ii++) {
					for (jj = 0; jj <= L; jj++) {
						if (hingestatus[ii][jj] == 0) {	//not ham
							isHam=0;
							break;
						}
					}
				}
				if(isHam==1){
					valid_2_spans++;                //This is a valid 2 span so it is counted
					//printf("Pattern number: %d \n", valid_2_spans);
					fillreordertemplate((*pointordNum)[0]-1, (*pointordNum)[1]-1);
					fillreordertemplate2((*pointordNum)[3]-1, (*pointordNum)[4]-1);
					recordtemplate(pointordNum);   //the information contained in ordertemplate is recorded
				}
			}

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
					enterhinge2(ii, jj, side, pointordNum,0);
				}
			}
		}
		(*pointordNum)[side+3]--;
			/* the edge number is decremented since, control is being passed back to  */
			/* the function that called leavehinge. That is to say, the program will  */
			/* now explore alternate ways of leaving the hinge (i.e. via different    */
			/* paths. this may be as simple as going out the other side). In essense  */
			/* the intial incremention above is nullified so when leavehinge is 	  */
			/* called again the leaving edge is numbered properly  			  */
		
		ordertemplate2[side][i][j] = 0;
			/* the edge that was recorded must now be erased as this edge may not exist  */
			/* as part of the SAW when alternate possibilities are explored. In the case */
			/* that it does not exist it needs to be numbered as "0" (zero).             */
	//		printf("finished leavehinge\n");


		curend2[0][num_walks2-1]=-1;
		curend2[1][num_walks2-1]=-1;
		curend2[2][num_walks2-1]=-1;
		curwalks2[curlength2-1][num_walks2-1] = 0;
		curlength2--;
	}
	return;
}
/***************************************************************************/

void rowedges2(int i, int j, int (*pointordNum)[6], int curlength2)  //y-direction
{
//	printf("rowedges2 called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(hingestatus[i][j-1]==0){
			hingestatus[i][j-1]=1;

			curlength2++;
			curwalks2[curlength2-1][num_walks2-1] = 4;

			(*pointordNum)[5]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			rowhingeedges[i][j-1]=1;
			leavehinge2(i, j-1, 0, pointordNum, curlength2);
				/* exit the hinge on side 0 with j moved*/
			leavehinge2(i, j-1, 1, pointordNum, curlength2);
				/* exit the hinge on side 1 with j moved*/
			rowedges2(i, j-1, pointordNum, curlength2);
				/* explore horizontal moves within the hinge */
			coledges2(i, j-1, pointordNum, curlength2);
				/* explore vertical moves within the hinge */
			(*pointordNum)[5]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i][j-1]=0;
			rowhingeedges[i][j-1]=0;

			curwalks2[curlength2-1][num_walks2-1] = 0;
			curlength2--;
		}
	}
	if(j<L){
		if(hingestatus[i][j+1]==0){
			hingestatus[i][j+1]=1;

			curlength2++;
			curwalks2[curlength2-1][num_walks2-1] = 3;

			(*pointordNum)[5]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			rowhingeedges[i][j]=1;
			leavehinge2(i, j+1, 0, pointordNum, curlength2);
				/* exit the hinge on side 0 with j moved*/
			leavehinge2(i, j+1, 1, pointordNum, curlength2);
				/* exit the hinge on side 1 with j moved*/
			rowedges2(i, j+1, pointordNum, curlength2);
				/* explore horizontal moves within the hinge */
			coledges2(i, j+1, pointordNum, curlength2);
				/* explore vertical moves within the hinge */
			(*pointordNum)[5]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i][j+1]=0;
			rowhingeedges[i][j]=0;

			curwalks2[curlength2-1][num_walks2-1] = 0;
			curlength2--;
		}
	}
//	printf("finished rowedges2\n");
	return;
}				/* end of function rowedges */

/***************************************************************************/

void coledges2(int i, int j, int (*pointordNum)[6], int curlength2)  //z-direction
/* called by rowedges,coledges */
{
//	printf("coledges2 called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(hingestatus[i-1][j]==0){
			hingestatus[i-1][j]=1;

			curlength2++;
			curwalks2[curlength2-1][num_walks2-1] = 6;

			(*pointordNum)[5]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			colhingeedges[i-1][j]=1;
			leavehinge2(i-1, j, 0, pointordNum, curlength2);
				/* exit the hinge on side 0 with j moved*/
			leavehinge2(i-1, j, 1, pointordNum, curlength2);
				/* exit the hinge on side 1 with j moved*/
			rowedges2(i-1, j, pointordNum, curlength2);
				/* explore horizontal moves within the hinge */
			coledges2(i-1, j, pointordNum, curlength2);
				/* explore vertical moves within the hinge */
			(*pointordNum)[5]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i-1][j]=0;
			colhingeedges[i-1][j]=0;

			curwalks2[curlength2-1][num_walks2-1] = 0;
			curlength2--;
		}
	}
	if(i<M){
		if(hingestatus[i+1][j]==0){
			hingestatus[i+1][j]=1;

			curlength2++;
			curwalks2[curlength2-1][num_walks2-1] = 5;

			(*pointordNum)[5]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			colhingeedges[i][j]=1;
			leavehinge2(i+1, j, 0, pointordNum, curlength2);
				/* exit the hinge on side 0 with j moved*/
			leavehinge2(i+1, j, 1, pointordNum, curlength2);
				/* exit the hinge on side 1 with j moved*/
			rowedges2(i+1, j, pointordNum, curlength2);
				/* explore horizontal moves within the hinge */
			coledges2(i+1, j, pointordNum, curlength2);
				/* explore vertical moves within the hinge */
			(*pointordNum)[5]--;
				/* decrement the number of edges in the hinge */
			hingestatus[i+1][j]=0;
			colhingeedges[i][j]=0;

			curwalks2[curlength2-1][num_walks2-1] = 0;
			curlength2--;
		}
	}
//	printf("finished coledges2\n");
	return;
}				/* end of function coledges */

/***************************************************************************/



















/***************************************************************************/
void
recordtemplate( int (*pointordNum)[6])
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
//	printf("OGoutNum=%d\n", outNum);
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
//	printf("NEWoutNum=%d\n", outNum);


	unsigned long int	inNum2 = num_section_norder2(0);
//	printf("OGinNum2=%d\n", inNum2);
	//Change inNum2 to it's key number.
	i=1;
	found=0;
	while(sectionkey[i]!=0){
		if(sectionkey[i]==inNum2){
			found=1;
			inNum2=i;
			break;
		}
		i++;
	}
	if(found==0){	//need to record inNum in sectionkey
		sectionkey[i]=inNum2;
		inNum2=i;
	}
//	printf("NEWinNum2=%d\n", inNum);

	unsigned long int	outNum2 = num_section_norder2(1);
//	printf("OGoutNum2=%d\n", outNum2);
	//reset outNum2 to it's key number.
	i=1;
	found=0;
	while(sectionkey[i]!=0){
		if(sectionkey[i]==outNum2){
			found=1;
			outNum2=i;
			break;
		}
		i++;
	}
	if(found==0){	//need to record inNum in sectionkey
		sectionkey[i]=outNum2;
		outNum2=i;
	}
//	printf("NEWoutNum2=%d\n", outNum2);


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
		if(test->outorder==outNum){	//if same 2 sections (for first poly)
			if(test->inorder2==inNum2 && test->outorder2==outNum2){	//if same 2 sections (second poly)
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
	}

//	printf("duplicate=%d\n", duplicate);
	if(duplicate==0){
		num_tspans++;

//		printf("actually recording\n");
		(*current_hinge_span[inNum]).nexthinge = newhinge(num_walks, num_walks2); /* add a newhinge to the linked list */
	
		current_hinge_span[inNum] = (*current_hinge_span[inNum]).nexthinge; /* let current_hinge_span[inNum] point to newhinge */
		(*current_hinge_span[inNum]).inorder2 = inNum2;	/* record the ordering on the out side */
		(*current_hinge_span[inNum]).outorder = outNum;	/* record the ordering on the out side */
		(*current_hinge_span[inNum]).outorder2 = outNum2;	/* record the ordering on the out side */
		(*current_hinge_span[inNum]).num_of_walks = num_walks;
		(*current_hinge_span[inNum]).num_of_walks2 = num_walks2;
		for(i=0; i<=vM*vL-1; i++){
			(*current_hinge_span[inNum]).hedges[i]=temp_hedges[i];
		}

		for(i=0; i<=num_walks-1; i++){
			for(j=0; j<=2; j++){
				current_hinge_span[inNum]->start[j][i] = curstart[j][i];
				current_hinge_span[inNum]->end[j][i] = curend[j][i];
			}
			for(j=0; j<=vM*vL; j++){
				current_hinge_span[inNum]->walks[j][i] = curwalks[j][i];
			}
		}

		for(i=0; i<=num_walks2-1; i++){
			for(j=0; j<=2; j++){
				current_hinge_span[inNum]->start2[j][i] = curstart2[j][i];
				current_hinge_span[inNum]->end2[j][i] = curend2[j][i];
			}
			for(j=0; j<=vM*vL; j++){
				current_hinge_span[inNum]->walks2[j][i] = curwalks2[j][i];
			}
		}

/*		printf("Recorded 2-span. Sections %d,%d to Sections %d,%d with hedges: ", inNum, inNum2, outNum, outNum2);
		for(i=0; i<=vM*vL-1; i++){
			printf("%d ", (*current_hinge_span[inNum]).hedges[i]);
		}
		printf("\n");
*/	}
	else{
		num_duplicate_tspans++;
	}
//	printf("finished recording\n");
	return;
}




///////////////ENDHINGE FUNCTIONS///////////////////////////////////////////////////////////////////


void
enterendhinge(int i, int j, int curlength)
{
	if( alreadyentered[i][j]!=1 ){
	//	printf("entered endhinge (i=%d, j=%d\n", i, j);

		curlength++;
		num_walks++;
		curstart[0][num_walks-1]=0;
		curstart[1][num_walks-1]=j;
		curstart[2][num_walks-1]=i;

		curwalks[curlength-1][num_walks-1] = 1;

		hingestatus[i][j]=1;
		endtemplate[i][j] = EndOrdNum[0];
		EndOrdNum[0]++;


		endrowedges(i, j, curlength);
			/* horizontal edges are explored */
		endcoledges(i, j, curlength);
			/* vertical edges are explored */

		//reset stuff
		EndOrdNum[0]--;
		endtemplate[i][j] = 0;
		hingestatus[i][j]=0;

		curstart[0][num_walks-1]=-1;
		curstart[1][num_walks-1]=-1;
		curstart[2][num_walks-1]=-1;
		curwalks[curlength-1][num_walks-1] = 0;

		curlength--;
		num_walks--;
	//	printf("finished enterhinge\n");
	}
	return;

}				/* end of function enterhinge */

/***************************************************************************/

void
leaveendhinge(int i, int j, int curlength)
/* called by: rowedges,coledges,enterhinge */
{
	if( alreadyentered[i][j]!=1 ){
	//	printf("left hinge (i=%d, j=%d\n", i, j);

		curlength++;
		curend[0][num_walks-1]=0;
		curend[1][num_walks-1]=j;
		curend[2][num_walks-1]=i;
		curwalks[curlength-1][num_walks-1] = 2;

		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		endtemplate[i][j] = EndOrdNum[0];
		EndOrdNum[0]++;

		if ( EndOrdNum[0]-1 <= vM*vL-2 ){	//make sure there's room for the 2nd poly to enter+exit
			for (ii = 0; ii <= M; ii++) { // these are the "vertical" locations of the vertices being considered
				for (jj = 0; jj <= L; jj++) { // these are the horizontal locations of the vertices being considered
					if (hingestatus[ii][jj] == 0) {
						// if the vertex is not occupied then explore entering the hinge
						// at this vertex with a second polygon
						enterendhinge2(ii, jj, 0);
					}
					alreadyentered2[ii][jj]=1;
				}
			}
			for (ii = 0; ii <= M; ii++) {
				for (jj = 0; jj <= L; jj++) {
					alreadyentered2[ii][jj]=0;
				}
			}
		}




		for (ii = 0; ii <= M; ii++) {
		/* these are the "vertical" locations of the vertices being considered */
			for (jj = 0; jj <= L; jj++) {
			/* these are the horizontal locations of the vertices being considered */
				if (hingestatus[ii][jj] == 0) {
				/* if the vertex is not occupied then explore entering the hinge at this vertex */
				/* this maintains the self avoiding nature of the walk */
					enterendhinge(ii, jj, 0);
				}
			}
		}

		EndOrdNum[0]--;
		endtemplate[i][j] = 0;

		curend[0][num_walks-1]=-1;
		curend[1][num_walks-1]=-1;
		curend[2][num_walks-1]=-1;
		curwalks[curlength-1][num_walks-1] = 0;
		curlength--;
	}
	return;
}
/***************************************************************************/

void endrowedges(int i, int j, int curlength)  //y-direction
{
//	printf("endrowedges called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(hingestatus[i][j-1]==0){
			hingestatus[i][j-1]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 4;

			EndOrdNum[1]++;
			rowhingeedges[i][j-1]=1;

			leaveendhinge(i, j-1, curlength);
			endrowedges(i, j-1, curlength);
			endcoledges(i, j-1, curlength);

			rowhingeedges[i][j-1]=0;
			EndOrdNum[1]--;
			hingestatus[i][j-1]=0;
			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
	if(j<L){
		if(hingestatus[i][j+1]==0){
			hingestatus[i][j+1]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 3;

			EndOrdNum[1]++;
			rowhingeedges[i][j]=1;

			leaveendhinge(i, j+1, curlength);
			endrowedges(i, j+1, curlength);
			endcoledges(i, j+1, curlength);

			rowhingeedges[i][j]=0;
			EndOrdNum[1]--;
			hingestatus[i][j+1]=0;
			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
//	printf("finished endrowedges\n");
	return;
}				/* end of function rowedges */

/***************************************************************************/

void endcoledges(int i, int j, int curlength)  //z-direction
{
//	printf("endcoledges called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(hingestatus[i-1][j]==0){
			hingestatus[i-1][j]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 6;

			EndOrdNum[1]++;
			colhingeedges[i-1][j]=1;

			leaveendhinge(i-1, j, curlength);
			endrowedges(i-1, j, curlength);
			endcoledges(i-1, j, curlength);

			colhingeedges[i-1][j]=0;
			EndOrdNum[1]--;
			hingestatus[i-1][j]=0;
			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
	if(i<M){
		if(hingestatus[i+1][j]==0){
			hingestatus[i+1][j]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 5;

			EndOrdNum[1]++;
			colhingeedges[i][j]=1;

			leaveendhinge(i+1, j, curlength);
			endrowedges(i+1, j, curlength);
			endcoledges(i+1, j, curlength);

			colhingeedges[i][j]=0;
			EndOrdNum[1]--;
			hingestatus[i+1][j]=0;
			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
//	printf("finished endcoledges\n");
	return;
}				/* end of function coledges */

/***************************************************************************/










void
enterendhinge2(int i, int j, int curlength2)
{
	if( alreadyentered2[i][j]!=1 ){
	//	printf("entered endhinge2 (i=%d, j=%d\n", i, j);

		curlength2++;
		num_walks2++;
		curstart2[0][num_walks2-1]=0;
		curstart2[1][num_walks2-1]=j;
		curstart2[2][num_walks2-1]=i;

		curwalks2[curlength2-1][num_walks2-1] = 1;

		hingestatus[i][j]=2;
		endtemplate3[i][j] = EndOrdNum2[0];
		EndOrdNum2[0]++;


		endrowedges2(i, j, curlength2);
			/* horizontal edges are explored */
		endcoledges2(i, j, curlength2);
			/* vertical edges are explored */

		//reset stuff
		EndOrdNum2[0]--;
		endtemplate3[i][j] = 0;
		hingestatus[i][j]=0;

		curstart2[0][num_walks2-1]=-1;
		curstart2[1][num_walks2-1]=-1;
		curstart2[2][num_walks2-1]=-1;
		curwalks2[curlength2-1][num_walks2-1] = 0;

		curlength2--;
		num_walks2--;
	//	printf("finished enterhinge2\n");
	}
	return;

}				/* end of function enterhinge */

/***************************************************************************/

void
leaveendhinge2(int i, int j, int curlength2)
/* called by: rowedges,coledges,enterhinge */
{
	if( alreadyentered2[i][j]!=1 ){
	//	printf("left hinge (i=%d, j=%d\n", i, j);

		curlength2++;
		curend2[0][num_walks2-1]=0;
		curend2[1][num_walks2-1]=j;
		curend2[2][num_walks2-1]=i;
		curwalks2[curlength2-1][num_walks2-1] = 2;

		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		endtemplate3[i][j] = EndOrdNum2[0];
		EndOrdNum2[0]++;

		num_potential_endhinges++;
		if (LFlag_endhinge2()){

			//Hamiltonian Check
			int isHam=1;
			for (ii = 0; ii <= M; ii++) {
				for (jj = 0; jj <= L; jj++) {
					if (hingestatus[ii][jj] == 0) {	//not ham
						isHam=0;
						break;
					}
				}
			}
			if(isHam==1){
				num_valid_endhinges++;
				fillrendtemplate();
				fillrendtemplate3();
				recordendtemplate(); /*the information contained in endtemplate is recorded*/
				//duplicate check is performed in recordtemplate.
			}


		}



		for (ii = 0; ii <= M; ii++) {
		/* these are the "vertical" locations of the vertices being considered */
			for (jj = 0; jj <= L; jj++) {
			/* these are the horizontal locations of the vertices being considered */
				if (hingestatus[ii][jj] == 0) {
				/* if the vertex is not occupied then explore entering the hinge at this vertex */
				/* this maintains the self avoiding nature of the walk */
					enterendhinge2(ii, jj, 0);
				}
			}
		}

		EndOrdNum2[0]--;
		endtemplate3[i][j] = 0;

		curend2[0][num_walks2-1]=-1;
		curend2[1][num_walks2-1]=-1;
		curend2[2][num_walks2-1]=-1;
		curwalks2[curlength2-1][num_walks2-1] = 0;
		curlength2--;
	}
	return;
}
/***************************************************************************/

void endrowedges2(int i, int j, int curlength2)  //y-direction
{
//	printf("endrowedges2 called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(hingestatus[i][j-1]==0){
			hingestatus[i][j-1]=2;

			curlength2++;
			curwalks2[curlength2-1][num_walks2-1] = 4;

			EndOrdNum2[1]++;
			rowhingeedges[i][j-1]=1;

			leaveendhinge2(i, j-1, curlength2);
			endrowedges2(i, j-1, curlength2);
			endcoledges2(i, j-1, curlength2);

			rowhingeedges[i][j-1]=0;
			EndOrdNum2[1]--;
			hingestatus[i][j-1]=0;
			curwalks2[curlength2-1][num_walks2-1] = 0;
			curlength2--;
		}
	}
	if(j<L){
		if(hingestatus[i][j+1]==0){
			hingestatus[i][j+1]=2;

			curlength2++;
			curwalks2[curlength2-1][num_walks2-1] = 3;

			EndOrdNum2[1]++;
			rowhingeedges[i][j]=1;

			leaveendhinge2(i, j+1, curlength2);
			endrowedges2(i, j+1, curlength2);
			endcoledges2(i, j+1, curlength2);

			rowhingeedges[i][j]=0;
			EndOrdNum2[1]--;
			hingestatus[i][j+1]=0;
			curwalks2[curlength2-1][num_walks2-1] = 0;
			curlength2--;
		}
	}
//	printf("finished endrowedges2\n");
	return;
}				/* end of function rowedges */

/***************************************************************************/

void endcoledges2(int i, int j, int curlength2)  //z-direction
{
//	printf("endcoledges2 called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(hingestatus[i-1][j]==0){
			hingestatus[i-1][j]=2;

			curlength2++;
			curwalks2[curlength2-1][num_walks2-1] = 6;

			EndOrdNum2[1]++;
			colhingeedges[i-1][j]=1;

			leaveendhinge2(i-1, j, curlength2);
			endrowedges2(i-1, j, curlength2);
			endcoledges2(i-1, j, curlength2);

			colhingeedges[i-1][j]=0;
			EndOrdNum2[1]--;
			hingestatus[i-1][j]=0;
			curwalks2[curlength2-1][num_walks2-1] = 0;
			curlength2--;
		}
	}
	if(i<M){
		if(hingestatus[i+1][j]==0){
			hingestatus[i+1][j]=2;

			curlength2++;
			curwalks2[curlength2-1][num_walks2-1] = 5;

			EndOrdNum2[1]++;
			colhingeedges[i][j]=1;

			leaveendhinge2(i+1, j, curlength2);
			endrowedges2(i+1, j, curlength2);
			endcoledges2(i+1, j, curlength2);

			colhingeedges[i][j]=0;
			EndOrdNum2[1]--;
			hingestatus[i+1][j]=0;
			curwalks2[curlength2-1][num_walks2-1] = 0;
			curlength2--;
		}
	}
//	printf("finished endcoledges2\n");
	return;
}				/* end of function coledges */

/***************************************************************************/












/***************************************************************************/
void
recordendtemplate()
{
//	printf("recording an endtemplate\n");
	int i,j;

	unsigned long int	secnum = num_section_endhinge_norder();
//	printf("OG secnum=%d\n", secnum);
	//Change secnum to it's key number.
	i=1;
	int found=0;
	while(sectionkey[i]!=0){
		if(sectionkey[i]==secnum){
			found=1;
			secnum=i;
			break;
		}
		i++;
	}
	if(found==0){	//This shouldn't happen, since all valid section should have been generated during 2-span geneartion.
		printf("another section was generated... problem. Exitting.\n");
		exit(1);
//		sectionkey[i]=secnum;
//		secnum=i;
	}
//	printf("NEW secnum=%lu\n", secnum);


	unsigned long int	secnum3 = num_section_endhinge_norder3();
//	printf("OG secnum3=%d\n", secnum3);
	//Change secnum to it's key number.
	i=1;
	found=0;
	while(sectionkey[i]!=0){
		if(sectionkey[i]==secnum3){
			found=1;
			secnum3=i;
			break;
		}
		i++;
	}
	if(found==0){	//This shouldn't happen, since all valid section should have been generated during 2-span geneartion.
		printf("another section was generated... problem. Exitting.\n");
		exit(1);
//		sectionkey[i]=secnum3;
//		secnum3=i;
	}
//	printf("NEW secnum3=%lu\n", secnum3);




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
/*	printf("there were %d edgenums\n", edgenum-1);
	printf("temp_hedges is: ");
	for(i=0; i<=vM*vL-1; i++){
		if(temp_hedges[i]==0){
			break;
		}
		printf("%d ", temp_hedges[i]);
	}
	printf("\n");
*/


	//check if this endhinge configuration has already been recorded
	int duplicate=0;
	int match;

	struct endhinge* test = firstendhinge[secnum];
	while(test->nextendhinge != NULL && duplicate==0){
//		printf("looking through linked list\n");
		test = test->nextendhinge;
		if(test->side==1 && test->sec2==secnum3){	//if right endhinge and same poly2 section.
//			printf("found one that needs to test\n");
			match=1;	//test if the hinge edges match
			for(i=0; i<=vM*vL-1; i++){
//				printf("temp_hedges[%d]=%d, test->hedges[%d]=%d\n", i, temp_hedges[i], i, test->hedges[i]);
				if(temp_hedges[i]!=test->hedges[i]){ //if an edge doesn't match
					match=0;
					break;
				}
			}
			if(match==1){
				duplicate=1;
				break;
			}
		}
	}
//	printf("duplicate=%d\n", duplicate);
	if(duplicate==0){
//		printf("actually recording right endhinge\n");
		tot_right_endhinges++;
		(*currentendhinge[secnum]).nextendhinge = newendhinge(num_walks, num_walks2); /* add a newendhinge to the linked list */
	
		currentendhinge[secnum] = (*currentendhinge[secnum]).nextendhinge;
		(*currentendhinge[secnum]).sec2 = secnum3;
		(*currentendhinge[secnum]).side = 1;
		(*currentendhinge[secnum]).num_of_walks = num_walks;
		(*currentendhinge[secnum]).num_of_walks2 = num_walks2;
		for(i=0; i<=vM*vL-1; i++){
			(*currentendhinge[secnum]).hedges[i]=temp_hedges[i];
		}

		for(i=0; i<=num_walks-1; i++){
			for(j=0; j<=2; j++){
				currentendhinge[secnum]->start[j][i] = curstart[j][i];
				currentendhinge[secnum]->end[j][i] = curend[j][i];
			}
			for(j=0; j<=vM*vL; j++){
				currentendhinge[secnum]->walks[j][i] = curwalks[j][i];
			}
		}

		for(i=0; i<=num_walks2-1; i++){
			for(j=0; j<=2; j++){
				currentendhinge[secnum]->start2[j][i] = curstart2[j][i];
				currentendhinge[secnum]->end2[j][i] = curend2[j][i];
			}
			for(j=0; j<=vM*vL; j++){
				currentendhinge[secnum]->walks2[j][i] = curwalks2[j][i];
			}
		}

	

		//now record the mirror image left endhinge version
		//first check if the left version will follow a lexicographical ordering on left endhinges: make poly1 be the first in the lexicographical ordering (of vertices)
		//lexi_left_endhinges counts number of left endhinges which have poly1 as the first (lexicographically) poly.

		int recordleft=0;	//0 if haven't found an occupied vertex yet. 1 if first vertex is from poly1. 2 if first vertex if from poly2

		for(i=0; i<=M; i++){
			for(j=0; j<=L; j++){
				if(recordleft==0){
					if(hingestatus[i][j]==1){
						recordleft=1;
					}
					else if(hingestatus[i][j]==2){
						recordleft=2;
					}
				}
				else{
					break;
				}
			}
		}
		

		if(recordleft==1){

			fillendtemplate2();
			fillrendtemplate2();
			fillendtemplate4();
			fillrendtemplate4();


			unsigned long int secnum2 = num_section_endhinge_norder2();
	//		printf("OG secnum2=%d\n", secnum2);
			//Change secnum2 to it's key number.
			i=1;
			found=0;
			while(sectionkey[i]!=0){
				if(sectionkey[i]==secnum2){
					found=1;
					secnum2=i;
					break;
				}
				i++;
			}
			if(found==0){	//This shouldn't happen, since all valid section should have been generated during 2-span geneartion.
				printf("another section was generated... problem. Exitting.\n");
				exit(1);
	//			sectionkey[i]=secnum2;
	//			secnum2=i;
			}
//			printf("NEW secnum2=%lu\n", secnum2);



			unsigned long int secnum4 = num_section_endhinge_norder4();
	//		printf("OG secnum4=%d\n", secnum4);
			//Change secnum4 to it's key number.
			i=1;
			found=0;
			while(sectionkey[i]!=0){
				if(sectionkey[i]==secnum4){
					found=1;
					secnum4=i;
					break;
				}
				i++;
			}
			if(found==0){	//This shouldn't happen, since all valid section should have been generated during 2-span geneartion.
				printf("another section was generated... problem. Exitting.\n");
				exit(1);
	//			sectionkey[i]=secnum4;
	//			secnum4=i;
			}
//			printf("NEW secnum4=%lu\n", secnum4);




			//Need to check if this left hinge will be a duplicate.

			//check if this endhinge configuration has already been recorded
			duplicate=0;
			test = firstendhinge[secnum2];
			while(test->nextendhinge != NULL && duplicate==0){
		//		printf("looking through linked list\n");
				test = test->nextendhinge;
				if(test->side==0 && test->sec2==secnum4){	//if left endhinge
		//			printf("found one that needs to test\n");
					match=1;	//test if the hinge edges match
					for(i=0; i<=vM*vL-1; i++){
		//				printf("temp_hedges[%d]=%d, test->hedges[%d]=%d\n", i, temp_hedges[i], i, test->hedges[i]);
						if(temp_hedges[i]!=test->hedges[i]){ //if an edge doesn't match
							match=0;
							break;
						}
					}
					if(match==1){
						duplicate=1;
						break;
					}
				}
			}
		//	printf("duplicate=%d\n", duplicate);

			if(duplicate==0){
				tot_left_endhinges++;
//				printf("recording left endhinge\n");
				(*currentendhinge[secnum2]).nextendhinge = newendhinge(num_walks, num_walks2); /* add a newendhinge to the linked list */

				currentendhinge[secnum2] = (*currentendhinge[secnum2]).nextendhinge;
				(*currentendhinge[secnum2]).sec2 = secnum4;
				(*currentendhinge[secnum2]).side = 0;
				(*currentendhinge[secnum2]).num_of_walks = num_walks;
				(*currentendhinge[secnum2]).num_of_walks2 = num_walks2;

				for(i=0; i<=vM*vL-1; i++){
					(*currentendhinge[secnum2]).hedges[i]=temp_hedges[i];
				}

				for(i=0; i<=num_walks-1; i++){
					currentendhinge[secnum2]->start[0][i] = 2;	//change
					currentendhinge[secnum2]->end[0][i] = 2;		//change
					for(j=1; j<=2; j++){
						currentendhinge[secnum2]->start[j][i] = curstart[j][i];
						currentendhinge[secnum2]->end[j][i] = curend[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						if(curwalks[j][i]==1){	//change
							currentendhinge[secnum2]->walks[j][i]=2;
						}
						else if(curwalks[j][i]==2){	//change
							currentendhinge[secnum2]->walks[j][i]=1;
						}
						else{
							currentendhinge[secnum2]->walks[j][i] = curwalks[j][i];
						}
					}
				}


				for(i=0; i<=num_walks2-1; i++){
					currentendhinge[secnum2]->start2[0][i] = 2;	//change
					currentendhinge[secnum2]->end2[0][i] = 2;		//change
					for(j=1; j<=2; j++){
						currentendhinge[secnum2]->start2[j][i] = curstart2[j][i];
						currentendhinge[secnum2]->end2[j][i] = curend2[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						if(curwalks2[j][i]==1){	//change
							currentendhinge[secnum2]->walks2[j][i]=2;
						}
						else if(curwalks2[j][i]==2){	//change
							currentendhinge[secnum2]->walks2[j][i]=1;
						}
						else{
							currentendhinge[secnum2]->walks2[j][i] = curwalks2[j][i];
						}
					}
				}
			}
			else{
				num_duplicate_left_endhinges++;
			}
		}
		else{	//poly2 was first poly accoridng to lexi ordering
			num_nonlexi_endhinges++;
		}		
	}
	else{
		num_duplicate_right_endhinges++;
	}
//	printf("finished recording\n");
	return;
}



///////////////////END OF ENDHINGE FUNCTIONS////////////////////////////////////////////////////////



/***************************************************************************/
struct hinge_span* newhinge(unsigned int numberofwalks, unsigned int numberofwalks2)		/* returns a pointer to a hinge_span structure */
/* this is a memory allocation function */
{
	struct hinge_span *nextnewhinge;

	nextnewhinge = ((struct hinge_span *) calloc(1, sizeof(struct hinge_span)));
	if (nextnewhinge == NULL) {
		fprintf(stderr, "unable to allocate memory");
		exit(1);
	}
	int i,j;
	for(i=0; i<=vM*vL-1; i++){
		nextnewhinge->hedges[i]=0;
	}

	//allocate space for start, end, and walks
	for(i=0; i<=2; i++){
		nextnewhinge->start[i] = int_vecalloc(0,numberofwalks-1);
		nextnewhinge->end[i] = int_vecalloc(0,numberofwalks-1);
		nextnewhinge->start2[i] = int_vecalloc(0,numberofwalks2-1);
		nextnewhinge->end2[i] = int_vecalloc(0,numberofwalks2-1);
	}
	for(i=0; i<=vM*vL; i++){
		nextnewhinge->walks[i] = int_vecalloc(0,numberofwalks-1);
		nextnewhinge->walks2[i] = int_vecalloc(0,numberofwalks2-1);
	}

	//initialize start, end, and walks to all zeros
	for(i=0; i<=numberofwalks-1; i++){
		for(j=0; j<=2; j++){
			nextnewhinge->start[j][i]=0;
			nextnewhinge->end[j][i]=0;
		}
		for(j=0; j<=vM*vL; j++){
			nextnewhinge->walks[j][i]=0;
		}
	}
	//initialize start, end, and walks to all zeros
	for(i=0; i<=numberofwalks2-1; i++){
		for(j=0; j<=2; j++){
			nextnewhinge->start2[j][i]=0;
			nextnewhinge->end2[j][i]=0;
		}
		for(j=0; j<=vM*vL; j++){
			nextnewhinge->walks2[j][i]=0;
		}
	}

	return nextnewhinge;
}


/***************************************************************************/
struct endhinge* newendhinge(unsigned int numberofwalks, unsigned int numberofwalks2)		/* returns a pointer to an endhinge structure */
/* this is a memory allocation function */
{
	struct endhinge *nextnewendhinge;

	nextnewendhinge = ((struct endhinge *) calloc(1, sizeof(struct endhinge)));
	if (nextnewendhinge == NULL) {
		fprintf(stderr, "unable to allocate memory");
		exit(1);
	}
	int i,j;

	for(i=0; i<=vM*vL-1; i++){
		nextnewendhinge->hedges[i]=0;
	}

	//allocate space for start, end, and walks
	for(i=0; i<=2; i++){
		nextnewendhinge->start[i] = int_vecalloc(0,numberofwalks-1);
		nextnewendhinge->end[i] = int_vecalloc(0,numberofwalks-1);
		nextnewendhinge->start2[i] = int_vecalloc(0,numberofwalks2-1);
		nextnewendhinge->end2[i] = int_vecalloc(0,numberofwalks2-1);
	}
	for(i=0; i<=vM*vL; i++){
		nextnewendhinge->walks[i] = int_vecalloc(0,numberofwalks-1);
		nextnewendhinge->walks2[i] = int_vecalloc(0,numberofwalks-1);
	}

	//initialize start, end, and walks to all zeros
	for(i=0; i<=numberofwalks-1; i++){
		for(j=0; j<=2; j++){
			nextnewendhinge->start[j][i]=0;
			nextnewendhinge->end[j][i]=0;
			nextnewendhinge->start2[j][i]=0;
			nextnewendhinge->end2[j][i]=0;
		}
		for(j=0; j<=vM*vL; j++){
			nextnewendhinge->walks[j][i]=0;
			nextnewendhinge->walks2[j][i]=0;
		}
	}

	return nextnewendhinge;
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

void fillreordertemplate2(int ledges, int redges){ //fills reordertemplate2 appropriately based on ordertemplate
	int i, j, ii, jj, s;
	int firstentry, secondentry;
	int connectingedge;

//	printf("fillreordertemplate2 called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			reordertemplate2[0][i][j]=0;
			reordertemplate2[1][i][j]=0;
		}
	}

//	printf("initialized\n");
//	printf("ledges=%d, redges=%d\n", ledges, redges);
	for(s=0; s<=1; s++){	//for each side
		firstentry=1;
		for(i=0; i<=M; i++){	
			for(j=0; j<=L; j++){
				if(ordertemplate2[s][i][j] % 2 != 0 && reordertemplate2[s][i][j]==0){	//if odd
//					printf("odd edge found in (%d,%d)\n", j, i);
					reordertemplate2[s][i][j] = firstentry;
//					printf("recorded edge %d in reorder\n", firstentry);
					if(ordertemplate2[s][i][j]==1){	//if 1, then have to find the location of ledges or redges.
						if(s==0){
							connectingedge=ledges;
						}
						else{
							connectingedge=redges;
						}
					}
					else{	//if not 1, then find location the edgenum-1
						connectingedge=ordertemplate2[s][i][j]-1;
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
							if(ordertemplate2[s][ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
								reordertemplate2[s][ii][jj] = secondentry;
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
				else if(ordertemplate2[s][i][j]>0 && reordertemplate2[s][i][j]==0){	// if even
//					printf("even edge found in (%d,%d)\n", j, i);
					reordertemplate2[s][i][j] = firstentry;
//					printf("recorded edge %d in reorder\n", firstentry);
					if(s==0 && ordertemplate2[s][i][j]==ledges){	//if max, then have to find the location of edge 1.
						connectingedge=1;
					}
					else if(s==1 && ordertemplate2[s][i][j]==redges){
						connectingedge=1;
					}
					else{	//if not max, then find location the edgenum+1
						connectingedge=ordertemplate2[s][i][j]+1;
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
							if(ordertemplate2[s][ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
								reordertemplate2[s][ii][jj] = secondentry;
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



void fillrendtemplate(){ //fills rendtemplate appropriately based on endtemplate
	int i, j, ii, jj;
	int firstentry, secondentry;
	int connectingedge;

//	printf("fillrendtemplate called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			rendtemplate[i][j]=0;
			rendtemplate[i][j]=0;
		}
	}

//	printf("initialized\n");
	firstentry=1;
	for(i=0; i<=M; i++){	
		for(j=0; j<=L; j++){
			if(endtemplate[i][j] % 2 != 0 && rendtemplate[i][j]==0){	//if odd
//				printf("odd edge found in (%d,%d)\n", j, i);
				rendtemplate[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(endtemplate[i][j]==1){	//if 1, then have to find the location of ledges or redges.
					connectingedge=EndOrdNum[0]-1;
				}
				else{	//if not 1, then find location the edgenum-1
					connectingedge=endtemplate[i][j]-1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=EndOrdNum[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(endtemplate[ii][jj]==connectingedge){
//							printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							rendtemplate[ii][jj] = secondentry;
//							printf("putting %d into that location\n", secondentry);
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
			else if(endtemplate[i][j]>0 && rendtemplate[i][j]==0){	// if even
//				printf("even edge found in (%d,%d)\n", j, i);
				rendtemplate[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(endtemplate[i][j]==EndOrdNum[0]-1){	//if max, then have to find the location of edge 1.
					connectingedge=1;
				}
				else{	//if not max, then find location the edgenum+1
					connectingedge=endtemplate[i][j]+1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=EndOrdNum[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(endtemplate[ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							rendtemplate[ii][jj] = secondentry;
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




void fillrendtemplate3(){ //fills rendtemplate3 appropriately based on endtemplate3
	int i, j, ii, jj;
	int firstentry, secondentry;
	int connectingedge;

//	printf("fillrendtemplate called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			rendtemplate3[i][j]=0;
			rendtemplate3[i][j]=0;
		}
	}

//	printf("initialized\n");
	firstentry=1;
	for(i=0; i<=M; i++){	
		for(j=0; j<=L; j++){
			if(endtemplate3[i][j] % 2 != 0 && rendtemplate3[i][j]==0){	//if odd
//				printf("odd edge found in (%d,%d)\n", j, i);
				rendtemplate3[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(endtemplate3[i][j]==1){	//if 1, then have to find the location of ledges or redges.
					connectingedge=EndOrdNum2[0]-1;
				}
				else{	//if not 1, then find location the edgenum-1
					connectingedge=endtemplate3[i][j]-1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=EndOrdNum2[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(endtemplate3[ii][jj]==connectingedge){
//							printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							rendtemplate3[ii][jj] = secondentry;
//							printf("putting %d into that location\n", secondentry);
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
			else if(endtemplate3[i][j]>0 && rendtemplate3[i][j]==0){	// if even
//				printf("even edge found in (%d,%d)\n", j, i);
				rendtemplate3[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(endtemplate3[i][j]==EndOrdNum2[0]-1){	//if max, then have to find the location of edge 1.
					connectingedge=1;
				}
				else{	//if not max, then find location the edgenum+1
					connectingedge=endtemplate3[i][j]+1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=EndOrdNum2[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(endtemplate3[ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							rendtemplate3[ii][jj] = secondentry;
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




void fillendtemplate2(){ //fills endtemplate2 (left endhinge) appropriately based on endtemplate (right endhinge)
	int i, j;
	int max=0;
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			endtemplate2[i][j] = 0;
			if(endtemplate[i][j]>max){
				max = endtemplate[i][j];
			}
		}
	}
	
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			if(endtemplate[i][j]>0 && endtemplate[i][j] < max){
				endtemplate2[i][j] = endtemplate[i][j] + 1;
			}
			else if(endtemplate[i][j]==max){
				endtemplate2[i][j]=1;
			}
		}
	}
}






void fillendtemplate4(){ //fills endtemplate4 (left endhinge) appropriately based on endtemplate3 (right endhinge)
	int i, j;
	int max=0;
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			endtemplate4[i][j] = 0;
			if(endtemplate3[i][j]>max){
				max = endtemplate3[i][j];
			}
		}
	}
	
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			if(endtemplate3[i][j]>0 && endtemplate3[i][j] < max){
				endtemplate4[i][j] = endtemplate3[i][j] + 1;
			}
			else if(endtemplate3[i][j]==max){
				endtemplate4[i][j]=1;
			}
		}
	}
}




void fillrendtemplate2(){ //fills rendtemplate2 appropriately based on endtemplate2
	int i, j, ii, jj;
	int firstentry, secondentry;
	int connectingedge;

//	printf("fillrendtemplate2 called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			rendtemplate2[i][j]=0;
			rendtemplate2[i][j]=0;
		}
	}

//	printf("initialized\n");
	firstentry=1;
	for(i=0; i<=M; i++){	
		for(j=0; j<=L; j++){
			if(endtemplate2[i][j] % 2 != 0 && rendtemplate2[i][j]==0){	//if odd
//				printf("odd edge found in (%d,%d)\n", j, i);
				rendtemplate2[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(endtemplate2[i][j]==1){	//if 1, then have to find the location of ledges or redges.
					connectingedge=EndOrdNum[0]-1;
				}
				else{	//if not 1, then find location the edgenum-1
					connectingedge=endtemplate2[i][j]-1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=EndOrdNum[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(endtemplate2[ii][jj]==connectingedge){
//							printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							rendtemplate2[ii][jj] = secondentry;
//							printf("putting %d into that location\n", secondentry);
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
			else if(endtemplate2[i][j]>0 && rendtemplate2[i][j]==0){	// if even
//				printf("even edge found in (%d,%d)\n", j, i);
				rendtemplate2[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(endtemplate2[i][j]==EndOrdNum[0]-1){	//if max, then have to find the location of edge 1.
					connectingedge=1;
				}
				else{	//if not max, then find location the edgenum+1
					connectingedge=endtemplate2[i][j]+1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=EndOrdNum[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(endtemplate2[ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							rendtemplate2[ii][jj] = secondentry;
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




void fillrendtemplate4(){ //fills rendtemplate4 appropriately based on endtemplate4
	int i, j, ii, jj;
	int firstentry, secondentry;
	int connectingedge;

//	printf("fillrendtemplate4 called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			rendtemplate4[i][j]=0;
			rendtemplate4[i][j]=0;
		}
	}

//	printf("initialized\n");
	firstentry=1;
	for(i=0; i<=M; i++){	
		for(j=0; j<=L; j++){
			if(endtemplate4[i][j] % 2 != 0 && rendtemplate4[i][j]==0){	//if odd
//				printf("odd edge found in (%d,%d)\n", j, i);
				rendtemplate4[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(endtemplate4[i][j]==1){	//if 1, then have to find the location of ledges or redges.
					connectingedge=EndOrdNum2[0]-1;
				}
				else{	//if not 1, then find location the edgenum-1
					connectingedge=endtemplate4[i][j]-1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=EndOrdNum2[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(endtemplate4[ii][jj]==connectingedge){
//							printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							rendtemplate4[ii][jj] = secondentry;
//							printf("putting %d into that location\n", secondentry);
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
			else if(endtemplate4[i][j]>0 && rendtemplate4[i][j]==0){	// if even
//				printf("even edge found in (%d,%d)\n", j, i);
				rendtemplate4[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(endtemplate4[i][j]==EndOrdNum2[0]-1){	//if max, then have to find the location of edge 1.
					connectingedge=1;
				}
				else{	//if not max, then find location the edgenum+1
					connectingedge=endtemplate4[i][j]+1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=EndOrdNum2[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(endtemplate4[ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							rendtemplate4[ii][jj] = secondentry;
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

void printbuiltwalks(){
	int i, j;
	printf("There are %d walks in built_walks:\n", num_built_walks);
	for(i=0; i<=num_built_walks-1; i++){
		printf("Start: (%d, %d, %d). ", built_walks_start[i][0], built_walks_start[i][1], built_walks_start[i][2]);
		printf("Direcs:");
		j=0;
		while(built_walks_direcs[i][j]!=0){
			printf(" %d", built_walks_direcs[i][j]);
			j++;
		}
		printf(". End: (%d, %d, %d)\n", built_walks_end[i][0], built_walks_end[i][1], built_walks_end[i][2]);
	}
}

void printbuiltwalks_all(){
	int i, j;
	printf("There are %d walks in built_walks (printingall):\n", num_built_walks);
	for(i=0; i<=num_built_walks-1; i++){
		printf("Start: (%d, %d, %d). ", built_walks_start[i][0], built_walks_start[i][1], built_walks_start[i][2]);
		printf("Direcs:");
		for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
			printf(" %d", built_walks_direcs[i][j]);
		}
		printf(". End: (%d, %d, %d)\n", built_walks_end[i][0], built_walks_end[i][1], built_walks_end[i][2]);
	}
}

void printbuiltwalks2_all(){
	int i, j;
	printf("There are %d walks in built_walks2 (printingall):\n", num_built_walks2);
	for(i=0; i<=num_built_walks2-1; i++){
		printf("Start: (%d, %d, %d). ", built_walks_start2[i][0], built_walks_start2[i][1], built_walks_start2[i][2]);
		printf("Direcs:");
		for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
			printf(" %d", built_walks_direcs2[i][j]);
		}
		printf(". End: (%d, %d, %d)\n", built_walks_end2[i][0], built_walks_end2[i][1], built_walks_end2[i][2]);
	}
}

void add_to_built_walks(unsigned long int secnum, int nth_tspan){
	//adds the walks in the tspan sec1->nth_tspan to built_walks

	//3 main cases for a walk in a 2span.
	//case 1a. starts left, ends right --> just append walk onto appropriate walk. num_built_walks doesn't change.
	//case 1b. starts right, ends left --> just append walk onto appropriate walk. num_built_walks doesn't change.
	//case 2. starts left, ends left --> will connect two pre-existing walks. num_built_walks decreases by one.
	//case 3. starts right, ends right --> will add a new walk to built_walks. num_built_walks increases by one.

	int i,j;

	int walknum;
	int walktoadd[vM*vL+1];
	int walktoadd_start[3];
	int walktoadd_end[3];

	int built_walknum;
	int length;

	//FIRST POLY
	for(walknum=0; walknum<=t_num_walks[secnum][nth_tspan]-1; walknum++){ //for each walk in the t-span
		//fill walktoadd
		for(i=0; i<=2; i++){
			walktoadd_start[i]=t_start[secnum][nth_tspan][i][walknum];
			walktoadd_end[i]=t_end[secnum][nth_tspan][i][walknum];
		}
		length=0;
		for(i=0; i<=vM*vL; i++){
			walktoadd[i] = t_walks[secnum][nth_tspan][i][walknum];
			if(walktoadd[i]!=0){
				length++;
			}
		}
/*
		printf("Adding this walk to built_walks: ");
		printf("Start: (%d, %d, %d). ", walktoadd_start[0], walktoadd_start[1], walktoadd_start[2]);
		printf("Direcs:");
		j=0;
		while(walktoadd[j]!=0){
			printf(" %d", walktoadd[j]);
			j++;
		}
		printf(". End: (%d, %d, %d,)", walktoadd_end[0], walktoadd_end[1], walktoadd_end[2]);
		printf("\n");
*/
		//3 cases.
		if(walktoadd_start[0]==0 && walktoadd_end[0]==2){	//case1a
//			printf("CASE 1a\n");
			//find which walk to append to
			for(built_walknum=0; built_walknum<=num_built_walks-1; built_walknum++){
				if(built_walks_start[built_walknum][1]==walktoadd_start[1] && built_walks_start[built_walknum][2]==walktoadd_start[2]){
					//built_walk starts at same place walktoadd starts (going in opposite direction)
					//add reversed version of walktoadd to the start of built_walks

					//shift built_walks_direcs up "length-1" entries.
					i=0;
					while(built_walks_direcs[built_walknum][i+1]!=0){
						i++;
					}
					while(i>=0){
						built_walks_direcs[built_walknum][i+length-1] = built_walks_direcs[built_walknum][i];
						i--;
					}
					//entries 0->length-2 are bogus now. Need to fill them with reverse version of walktoadd
					for(i=0; i<=length-2; i++){
						built_walks_direcs[built_walknum][i] = reverse_direc(walktoadd[length-1-i]);
					}

					//change start point
					built_walks_start[built_walknum][1] = walktoadd_end[1];
					built_walks_start[built_walknum][2] = walktoadd_end[2];
				}
				else if(built_walks_end[built_walknum][1]==walktoadd_start[1] && built_walks_end[built_walknum][2]==walktoadd_start[2]){
					//built_walk ends at same place walktoadd starts (going in same direction)

					i=0;
					while(built_walks_direcs[built_walknum][i]!=0){
						i++;
					}
					j=1;	//Note: skip overlapping first edge.
					while(walktoadd[j]!=0){
						built_walks_direcs[built_walknum][i] = walktoadd[j];
						i++;
						j++;
					}

					//change end point
					built_walks_end[built_walknum][1] = walktoadd_end[1];
					built_walks_end[built_walknum][2] = walktoadd_end[2];
				}
			}
		}
		else if(walktoadd_start[0]==2 && walktoadd_end[0]==0){	//case1b
//			printf("CASE 1b\n");
			//find which walk to append to
			for(built_walknum=0; built_walknum<=num_built_walks-1; built_walknum++){
				if(built_walks_start[built_walknum][1]==walktoadd_end[1] && built_walks_start[built_walknum][2]==walktoadd_end[2]){
					//built_walk starts at same place walktoadd ends (going in same direction)
					//add walktoadd to the start of built_walks

					//shift built_walks_direcs up "length-1" entries.
					i=0;
					while(built_walks_direcs[built_walknum][i+1]!=0){
						i++;
					}
					while(i>=0){
						built_walks_direcs[built_walknum][i+length-1] = built_walks_direcs[built_walknum][i];
						i--;
					}
					//entries 0->length-2 are bogus now. Need to fill them with walktoadd
					for(i=0; i<=length-2; i++){
						built_walks_direcs[built_walknum][i] = walktoadd[i];
					}

					//change start point
					built_walks_start[built_walknum][1] = walktoadd_start[1];
					built_walks_start[built_walknum][2] = walktoadd_start[2];
				}
				else if(built_walks_end[built_walknum][1]==walktoadd_end[1] && built_walks_end[built_walknum][2]==walktoadd_end[2]){
					//built_walk ends at same place walktoadd ends (going in opposite direction)

					i=0;
					while(built_walks_direcs[built_walknum][i]!=0){
						i++;
					}
					j=length-2; //Note: skip overlapping first edge.
					while(j>=0){
						built_walks_direcs[built_walknum][i] = reverse_direc(walktoadd[j]);
						i++;
						j--;
					}

					//change end point
					built_walks_end[built_walknum][1] = walktoadd_start[1];
					built_walks_end[built_walknum][2] = walktoadd_start[2];
				}
			}
		}
		else if(walktoadd_start[0]==0 && walktoadd_end[0]==0){	//case2
//			printf("CASE 2\n");	//most complicated case.
			//this will connect two existing walks. walk->walktoadd->walk. will reduce number of built_walks

			//ALGORITHM:
			//Find walk numbers that will connect (call them a,b) a connects to walktoadd_start, b connects to walktoadd_end.
			//Fill array "newwalk" which will hold the final product. must cover all cases for ordering.
				//fill with a, then walktoadd, then b.
			//replace a with newwalk.
			//move walk b to end of array, then delete
			//num_built_walks--

			int a=-1;
			int b=-1;
			int newwalk_start[3];
			int newwalk_end[3];
			int newwalk[vM*vL*(totalspan+1)];

			//find a and b.
			for(built_walknum=0; built_walknum<=num_built_walks-1; built_walknum++){
				if(built_walks_start[built_walknum][1]==walktoadd_start[1] && built_walks_start[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}
				else if(built_walks_end[built_walknum][1]==walktoadd_start[1] && built_walks_end[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}

				else if(built_walks_start[built_walknum][1]==walktoadd_end[1] && built_walks_start[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
				else if(built_walks_end[built_walknum][1]==walktoadd_end[1] && built_walks_end[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
			}


			if(a<0 || b<0){	//invalid--will close off a poly. return 0 (which means fail).
				printf("invalid connection. exitting\n");
				exit(1);
			}

			//fill newwalk with proper version of a. (start at opposite end of where walktoend starts)
			if(built_walks_end[a][1]==walktoadd_start[1] && built_walks_end[a][2]==walktoadd_start[2]){
				//correct direction. easy.
				for(i=0; i<=2; i++){
					newwalk_start[i]=built_walks_start[a][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					newwalk[i] = built_walks_direcs[a][i];
				}
			}
			else{
				//need to reverse a
				for(i=0; i<=2; i++){
					newwalk_start[i]=built_walks_end[a][i];
				}
				i=0;
				while(built_walks_direcs[a][i+1]!=0){
					i++;
				}
				j=0;
				while(i>=0){
					newwalk[j] = reverse_direc(built_walks_direcs[a][i]);
					j++;
					i--;
				}
				while(j<=vM*vL*(totalspan+1)-1){
					newwalk[j]=0;
					j++;
				}
			}
/*
printf("Newwalk Start: (%d, %d, %d). ", newwalk_start[0], newwalk_start[1], newwalk_start[2]);
printf("Newwalk Direcs:");
for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
	printf(" %d", newwalk[j]);
}
printf(". Newwalk End: (%d, %d, %d,)\n", newwalk_end[0], newwalk_end[1], newwalk_end[2]);
*/
			//append walktoadd to newwalk.
			i=0;
			while(newwalk[i]!=0){
				i++;
			}
			j=1;	//Note: skip overlapping first edge.
			while(walktoadd[j]!=0){
				newwalk[i] = walktoadd[j];
				i++;
				j++;
			}

			//append proper version of b to newwalk. (starting at newwalk[i])
			if(built_walks_start[b][1]==walktoadd_end[1] && built_walks_start[b][2]==walktoadd_end[2]){
				//correct direction. easy.
				//append b to newwalk. //i is already in the correct position.
				j=1;	//Note: skip overlapping first edge.
				while(built_walks_direcs[b][j]!=0){
					newwalk[i] = built_walks_direcs[b][j];
					i++;
					j++;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = built_walks_end[b][i];
				}
			}
			else{
				//need to append reversed version of b
				j=0;
				while(built_walks_direcs[b][j+1]!=0){
					j++;
				}
				j--; //Note: skip overlapping first edge.
				while(j>=0){
					newwalk[i] = reverse_direc(built_walks_direcs[b][j]);
					i++;
					j--;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = built_walks_start[b][i];
				}
			}

			//replace a with newwalk.
			for(i=0; i<=2; i++){
				built_walks_start[a][i] = newwalk_start[i];
				built_walks_end[a][i] = newwalk_end[i];
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				built_walks_direcs[a][i] = newwalk[i];
			}

			//move walk b to end of array
			while(b<num_built_walks-1){
				for(i=0; i<=2; i++){
					built_walks_start[b][i] = built_walks_start[b+1][i];
					built_walks_end[b][i] = built_walks_end[b+1][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					built_walks_direcs[b][i] = built_walks_direcs[b+1][i];
				}
				b++;
			}

			//delete b
			for(i=0; i<=2; i++){
				built_walks_start[b][i] = -1;
				built_walks_end[b][i] = -1;
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				built_walks_direcs[b][i] = 0;
			}

			num_built_walks--;
		}
		else if(walktoadd_start[0]==2 && walktoadd_end[0]==2){	//case3
//			printf("CASE 3\n");
			//add a new walk built_walks (walktoadd)
			num_built_walks++;
			for(i=0; i<=2; i++){
				built_walks_start[num_built_walks-1][i] = walktoadd_start[i];
				built_walks_end[num_built_walks-1][i] = walktoadd_end[i];
			}
			for(i=0; i<=vM*vL; i++){
				built_walks_direcs[num_built_walks-1][i] = walktoadd[i];
			}
		}
		else{
			printf("Problem when adding a tspan. Exitting\n");
			exit(1);
		}
//		printf("one walk has been added, here is the result:\n");
//		printbuiltwalks_all();
	}


////////////////////////////////////////////

	//SECOND POLY
	for(walknum=0; walknum<=t_num_walks2[secnum][nth_tspan]-1; walknum++){ //for each walk in the t-span
		//fill walktoadd
		for(i=0; i<=2; i++){
			walktoadd_start[i]=t_start2[secnum][nth_tspan][i][walknum];
			walktoadd_end[i]=t_end2[secnum][nth_tspan][i][walknum];
		}
		length=0;
		for(i=0; i<=vM*vL; i++){
			walktoadd[i] = t_walks2[secnum][nth_tspan][i][walknum];
			if(walktoadd[i]!=0){
				length++;
			}
		}
/*
		printf("Adding this walk to built_walks2: ");
		printf("Start: (%d, %d, %d). ", walktoadd_start[0], walktoadd_start[1], walktoadd_start[2]);
		printf("Direcs:");
		j=0;
		while(walktoadd[j]!=0){
			printf(" %d", walktoadd[j]);
			j++;
		}
		printf(". End: (%d, %d, %d,)", walktoadd_end[0], walktoadd_end[1], walktoadd_end[2]);
		printf("\n");
*/
		//3 cases.
		if(walktoadd_start[0]==0 && walktoadd_end[0]==2){	//case1a
//			printf("CASE 1a\n");
			//find which walk to append to
			for(built_walknum=0; built_walknum<=num_built_walks2-1; built_walknum++){
				if(built_walks_start2[built_walknum][1]==walktoadd_start[1] && built_walks_start2[built_walknum][2]==walktoadd_start[2]){
					//built_walk2 starts at same place walktoadd starts (going in opposite direction)
					//add reversed version of walktoadd to the start of built_walks2

					//shift built_walks_direcs2 up "length-1" entries.
					i=0;
					while(built_walks_direcs2[built_walknum][i+1]!=0){
						i++;
					}
					while(i>=0){
						built_walks_direcs2[built_walknum][i+length-1] = built_walks_direcs2[built_walknum][i];
						i--;
					}
					//entries 0->length-2 are bogus now. Need to fill them with reverse version of walktoadd
					for(i=0; i<=length-2; i++){
						built_walks_direcs2[built_walknum][i] = reverse_direc(walktoadd[length-1-i]);
					}

					//change start point
					built_walks_start2[built_walknum][1] = walktoadd_end[1];
					built_walks_start2[built_walknum][2] = walktoadd_end[2];
				}
				else if(built_walks_end2[built_walknum][1]==walktoadd_start[1] && built_walks_end2[built_walknum][2]==walktoadd_start[2]){
					//built_walk2 ends at same place walktoadd starts (going in same direction)

					i=0;
					while(built_walks_direcs2[built_walknum][i]!=0){
						i++;
					}
					j=1;	//Note: skip overlapping first edge.
					while(walktoadd[j]!=0){
						built_walks_direcs2[built_walknum][i] = walktoadd[j];
						i++;
						j++;
					}

					//change end point
					built_walks_end2[built_walknum][1] = walktoadd_end[1];
					built_walks_end2[built_walknum][2] = walktoadd_end[2];
				}
			}
		}
		else if(walktoadd_start[0]==2 && walktoadd_end[0]==0){	//case1b
//			printf("CASE 1b\n");
			//find which walk to append to
			for(built_walknum=0; built_walknum<=num_built_walks2-1; built_walknum++){
				if(built_walks_start2[built_walknum][1]==walktoadd_end[1] && built_walks_start2[built_walknum][2]==walktoadd_end[2]){
					//built_walk2 starts at same place walktoadd ends (going in same direction)
					//add walktoadd to the start of built_walks2

					//shift built_walks_direcs2 up "length-1" entries.
					i=0;
					while(built_walks_direcs2[built_walknum][i+1]!=0){
						i++;
					}
					while(i>=0){
						built_walks_direcs2[built_walknum][i+length-1] = built_walks_direcs2[built_walknum][i];
						i--;
					}
					//entries 0->length-2 are bogus now. Need to fill them with walktoadd
					for(i=0; i<=length-2; i++){
						built_walks_direcs2[built_walknum][i] = walktoadd[i];
					}

					//change start point
					built_walks_start2[built_walknum][1] = walktoadd_start[1];
					built_walks_start2[built_walknum][2] = walktoadd_start[2];
				}
				else if(built_walks_end2[built_walknum][1]==walktoadd_end[1] && built_walks_end2[built_walknum][2]==walktoadd_end[2]){
					//built_walk ends at same place walktoadd ends (going in opposite direction)

					i=0;
					while(built_walks_direcs2[built_walknum][i]!=0){
						i++;
					}
					j=length-2; //Note: skip overlapping first edge.
					while(j>=0){
						built_walks_direcs2[built_walknum][i] = reverse_direc(walktoadd[j]);
						i++;
						j--;
					}

					//change end point
					built_walks_end2[built_walknum][1] = walktoadd_start[1];
					built_walks_end2[built_walknum][2] = walktoadd_start[2];
				}
			}
		}
		else if(walktoadd_start[0]==0 && walktoadd_end[0]==0){	//case2
//			printf("CASE 2\n");	//most complicated case.
			//this will connect two existing walks. walk->walktoadd->walk. will reduce number of built_walks2

			//ALGORITHM:
			//Find walk numbers that will connect (call them a,b) a connects to walktoadd_start, b connects to walktoadd_end.
			//Fill array "newwalk" which will hold the final product. must cover all cases for ordering.
				//fill with a, then walktoadd, then b.
			//replace a with newwalk.
			//move walk b to end of array, then delete
			//num_built_walks2--

			int a=-1;
			int b=-1;
			int newwalk_start[3];
			int newwalk_end[3];
			int newwalk[vM*vL*(totalspan+1)];

			//find a and b.
			for(built_walknum=0; built_walknum<=num_built_walks2-1; built_walknum++){
				if(built_walks_start2[built_walknum][1]==walktoadd_start[1] && built_walks_start2[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}
				else if(built_walks_end2[built_walknum][1]==walktoadd_start[1] && built_walks_end2[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}

				else if(built_walks_start2[built_walknum][1]==walktoadd_end[1] && built_walks_start2[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
				else if(built_walks_end2[built_walknum][1]==walktoadd_end[1] && built_walks_end2[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
			}

//			printf("a=%d,b=%d\n", a, b);

			if(a<0 || b<0){	//invalid--will close off a seperate poly that doesn't span tube. return 0 (which means fail).
				printf("invalid connection (2nd poly). exitting\n");
				exit(1);
			}

			//fill newwalk with proper version of a. (start at opposite end of where walktoend starts)
			if(built_walks_end2[a][1]==walktoadd_start[1] && built_walks_end2[a][2]==walktoadd_start[2]){
				//correct direction. easy.
				for(i=0; i<=2; i++){
					newwalk_start[i]=built_walks_start2[a][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					newwalk[i] = built_walks_direcs2[a][i];
				}
			}
			else{
				//need to reverse a
				for(i=0; i<=2; i++){
					newwalk_start[i]=built_walks_end2[a][i];
				}
				i=0;
				while(built_walks_direcs2[a][i+1]!=0){
					i++;
				}
				j=0;
				while(i>=0){
					newwalk[j] = reverse_direc(built_walks_direcs2[a][i]);
					j++;
					i--;
				}
				while(j<=vM*vL*(totalspan+1)-1){
					newwalk[j]=0;
					j++;
				}
			}
/*
printf("Newwalk Start: (%d, %d, %d). ", newwalk_start[0], newwalk_start[1], newwalk_start[2]);
printf("Newwalk Direcs:");
for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
	printf(" %d", newwalk[j]);
}
printf(". Newwalk End: (%d, %d, %d,)\n", newwalk_end[0], newwalk_end[1], newwalk_end[2]);
*/
			//append walktoadd to newwalk.
			i=0;
			while(newwalk[i]!=0){
				i++;
			}
			j=1;	//Note: skip overlapping first edge.
			while(walktoadd[j]!=0){
				newwalk[i] = walktoadd[j];
				i++;
				j++;
			}

			//append proper version of b to newwalk. (starting at newwalk[i])
			if(built_walks_start2[b][1]==walktoadd_end[1] && built_walks_start2[b][2]==walktoadd_end[2]){
				//correct direction. easy.
				//append b to newwalk. //i is already in the correct position.
				j=1;	//Note: skip overlapping first edge.
				while(built_walks_direcs2[b][j]!=0){
					newwalk[i] = built_walks_direcs2[b][j];
					i++;
					j++;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = built_walks_end2[b][i];
				}
			}
			else{
				//need to append reversed version of b
				j=0;
				while(built_walks_direcs2[b][j+1]!=0){
					j++;
				}
				j--; //Note: skip overlapping first edge.
				while(j>=0){
					newwalk[i] = reverse_direc(built_walks_direcs2[b][j]);
					i++;
					j--;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = built_walks_start2[b][i];
				}
			}

			//replace a with newwalk.
			for(i=0; i<=2; i++){
				built_walks_start2[a][i] = newwalk_start[i];
				built_walks_end2[a][i] = newwalk_end[i];
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				built_walks_direcs2[a][i] = newwalk[i];
			}

			//move walk b to end of array
			while(b<num_built_walks2-1){
				for(i=0; i<=2; i++){
					built_walks_start2[b][i] = built_walks_start2[b+1][i];
					built_walks_end2[b][i] = built_walks_end2[b+1][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					built_walks_direcs2[b][i] = built_walks_direcs2[b+1][i];
				}
				b++;
			}

			//delete b
			for(i=0; i<=2; i++){
				built_walks_start2[b][i] = -1;
				built_walks_end2[b][i] = -1;
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				built_walks_direcs2[b][i] = 0;
			}

			num_built_walks2--;
		}
		else if(walktoadd_start[0]==2 && walktoadd_end[0]==2){	//case3
//			printf("CASE 3\n");
			//add a new walk built_walks2 (walktoadd)
			num_built_walks2++;
			for(i=0; i<=2; i++){
				built_walks_start2[num_built_walks2-1][i] = walktoadd_start[i];
				built_walks_end2[num_built_walks2-1][i] = walktoadd_end[i];
			}
			for(i=0; i<=vM*vL; i++){
				built_walks_direcs2[num_built_walks2-1][i] = walktoadd[i];
			}
		}
		else{
			printf("Problem when adding a tspan. Exitting\n");
			exit(1);
		}
//		printf("one walk has been added, here is the result:\n");
//		printbuiltwalks2_all();
	}


}

void add_right_endhinge(unsigned long int secnum, int nth_endhinge){
	//Each walk in the right endhinge will connect two existing walks. walk->walktoadd->walk.
	//Done until there is only one walk left.

	//ALGORITHM:
	//For num_built_walks>1:
		//Find walk numbers that will connect (call them a,b) a connects to walktoadd_start, b connects to walktoadd_end.
		//Fill array "newwalk" which will hold the final product. must cover all cases for ordering.
			//fill with a, then walktoadd, then b.
		//replace a with newwalk.
		//move walk b to end of array, then delete
		//num_built_walks--

	//When num_built_walks==1, finish up SAP.

	int i,j;

	int walknum;
	int walktoadd[vM*vL+1];
	int walktoadd_start[3];
	int walktoadd_end[3];

	int built_walknum;
	int length;

	int a=-1;
	int b=-1;
	int newwalk_start[3];
	int newwalk_end[3];
	int newwalk[vM*vL*(totalspan+1)];

	for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
		newwalk[i]=0;
	}
/*
	printf("Rend_num_walks=%d\n", Rend_num_walks[secnum][nth_endhinge]);
	printf("num_built_walks=%d\n", num_built_walks);
*/
	for(walknum=0; walknum<=Rend_num_walks[secnum][nth_endhinge]-1; walknum++){ //for each walk in the right endhinge
		//fill walktoadd
		for(i=0; i<=2; i++){
			walktoadd_start[i]=Rend_start[secnum][nth_endhinge][i][walknum];
			walktoadd_end[i]=Rend_end[secnum][nth_endhinge][i][walknum];
		}
		length=0;
		for(i=0; i<=vM*vL; i++){
			walktoadd[i] = Rend_walks[secnum][nth_endhinge][i][walknum];
			if(walktoadd[i]!=0){
				length++;
			}
		}

		if(num_built_walks>1){
			//find a and b.
			for(built_walknum=0; built_walknum<=num_built_walks-1; built_walknum++){
				if(built_walks_start[built_walknum][1]==walktoadd_start[1] && built_walks_start[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}
				else if(built_walks_end[built_walknum][1]==walktoadd_start[1] && built_walks_end[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}

				else if(built_walks_start[built_walknum][1]==walktoadd_end[1] && built_walks_start[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
				else if(built_walks_end[built_walknum][1]==walktoadd_end[1] && built_walks_end[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
			}

			if(a<0 || b<0){
				printf("didn't find a or b. Exitting\n");
				exit(1);
			}
			//fill newwalk with proper version of a. (start at opposite end of where walktoend starts)
			if(built_walks_end[a][1]==walktoadd_start[1] && built_walks_end[a][2]==walktoadd_start[2]){
				//correct direction. easy.
				for(i=0; i<=2; i++){
					newwalk_start[i]=built_walks_start[a][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					newwalk[i] = built_walks_direcs[a][i];
				}
			}
			else{
				//need to reverse a
				for(i=0; i<=2; i++){
					newwalk_start[i]=built_walks_end[a][i];
				}
				i=0;
				while(built_walks_direcs[a][i+1]!=0){
					i++;
				}
				j=0;
				while(i>=0){
					newwalk[j] = reverse_direc(built_walks_direcs[a][i]);
					j++;
					i--;
				}
				while(j<=vM*vL*(totalspan+1)-1){
					newwalk[j]=0;
					j++;
				}
			}

			//append walktoadd to newwalk.
			i=0;
			while(newwalk[i]!=0){
				i++;
			}
			j=1;	//Note: skip overlapping first edge.
			while(walktoadd[j]!=0){
				newwalk[i] = walktoadd[j];
				i++;
				j++;
			}

			//append proper version of b to newwalk. (starting at newwalk[i])
			if(built_walks_start[b][1]==walktoadd_end[1] && built_walks_start[b][2]==walktoadd_end[2]){
				//correct direction. easy.
				//append b to newwalk. //i is already in the correct position.
				j=1;	//Note: skip overlapping first edge.
				while(built_walks_direcs[b][j]!=0){
					newwalk[i] = built_walks_direcs[b][j];
					i++;
					j++;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = built_walks_end[b][i];
				}
			}
			else{
				//need to append reversed version of b
				j=0;
				while(built_walks_direcs[b][j+1]!=0){
					j++;
				}
				j--; //Note: skip overlapping first edge.
				while(j>=0){
					newwalk[i] = reverse_direc(built_walks_direcs[b][j]);
					i++;
					j--;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = built_walks_start[b][i];
				}
			}

			//replace a with newwalk.
			for(i=0; i<=2; i++){
				built_walks_start[a][i] = newwalk_start[i];
				built_walks_end[a][i] = newwalk_end[i];
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				built_walks_direcs[a][i] = newwalk[i];
			}

			//move walk b to end of array
			while(b<num_built_walks-1){
				for(i=0; i<=2; i++){
					built_walks_start[b][i] = built_walks_start[b+1][i];
					built_walks_end[b][i] = built_walks_end[b+1][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					built_walks_direcs[b][i] = built_walks_direcs[b+1][i];
				}
				b++;
			}

			//delete b
			for(i=0; i<=2; i++){
				built_walks_start[b][i] = -1;
				built_walks_end[b][i] = -1;
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				built_walks_direcs[b][i] = 0;
			}
			num_built_walks--;
		}
		else{	//close up SAP.
/*
printf("\ntime to close up the SAP:\n");
printbuiltwalks_all();
printf("walktoadd:");
i=0;
while(walktoadd[i]!=0){
	printf(" %d", walktoadd[i]);
	i++;
}
printf("\n");
*/			//either append walktoadd in same order or reverse
			i=0;
			while(built_walks_direcs[0][i]!=0){
				i++;
			}
			if(built_walks_end[0][1]==walktoadd_start[1] && built_walks_end[0][2]==walktoadd_start[2]){
				//just append
				j=1; //Note: skip overlapping first edge (which will be a 1).
				while(walktoadd[j]!=2){	//also skip overlapping last edge (which will be a 2)
					built_walks_direcs[0][i]=walktoadd[j];
					i++;
					j++;
				}
			}
			else{
				//append reverse of walktoadd
				j=0;
				while(walktoadd[j+1]!=2){	//Note: skip overlapping first edge (which will be a 2)
					j++;
				}
				while(j>=1){ //also skip overlapping last edge (which will be a 1)
					built_walks_direcs[0][i] = reverse_direc(walktoadd[j]);
					i++;
					j--;
				}
			}
			//set actual startpoint of SAP (and equivalent endpoint)
			built_walks_start[0][0]=totalspan;
			for(i=0; i<=2; i++){
				built_walks_end[0][i]=built_walks_start[0][i];
				built_walks_end[0][i]=built_walks_start[0][i];
			}
		}
/*
printf("%d-th part of endhinge added:\n", walknum+1);
printbuiltwalks_all();
*/
	}






	///////////2ND POLY

	for(walknum=0; walknum<=Rend_num_walks2[secnum][nth_endhinge]-1; walknum++){ //for each walk in the right endhinge
		//fill walktoadd
		for(i=0; i<=2; i++){
			walktoadd_start[i]=Rend_start2[secnum][nth_endhinge][i][walknum];

			walktoadd_end[i]=Rend_end2[secnum][nth_endhinge][i][walknum];
		}
		length=0;
		for(i=0; i<=vM*vL; i++){
			walktoadd[i] = Rend_walks2[secnum][nth_endhinge][i][walknum];
			if(walktoadd[i]!=0){
				length++;
			}
		}

		if(num_built_walks2>1){
			//find a and b.
			for(built_walknum=0; built_walknum<=num_built_walks2-1; built_walknum++){
				if(built_walks_start2[built_walknum][1]==walktoadd_start[1] && built_walks_start2[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}
				else if(built_walks_end2[built_walknum][1]==walktoadd_start[1] && built_walks_end2[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}

				else if(built_walks_start2[built_walknum][1]==walktoadd_end[1] && built_walks_start2[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
				else if(built_walks_end2[built_walknum][1]==walktoadd_end[1] && built_walks_end2[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
			}

			if(a<0 || b<0){
				printf("didn't find a or b (2ndpoly). Exitting\n");
				exit(1);
			}
			//fill newwalk with proper version of a. (start at opposite end of where walktoend starts)
			if(built_walks_end2[a][1]==walktoadd_start[1] && built_walks_end2[a][2]==walktoadd_start[2]){
				//correct direction. easy.
				for(i=0; i<=2; i++){
					newwalk_start[i]=built_walks_start2[a][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					newwalk[i] = built_walks_direcs2[a][i];
				}
			}
			else{
				//need to reverse a
				for(i=0; i<=2; i++){
					newwalk_start[i]=built_walks_end2[a][i];
				}
				i=0;
				while(built_walks_direcs2[a][i+1]!=0){
					i++;
				}
				j=0;
				while(i>=0){
					newwalk[j] = reverse_direc(built_walks_direcs2[a][i]);
					j++;
					i--;
				}
				while(j<=vM*vL*(totalspan+1)-1){
					newwalk[j]=0;
					j++;
				}
			}

			//append walktoadd to newwalk.
			i=0;
			while(newwalk[i]!=0){
				i++;
			}
			j=1;	//Note: skip overlapping first edge.
			while(walktoadd[j]!=0){
				newwalk[i] = walktoadd[j];
				i++;
				j++;
			}

			//append proper version of b to newwalk. (starting at newwalk[i])
			if(built_walks_start2[b][1]==walktoadd_end[1] && built_walks_start2[b][2]==walktoadd_end[2]){
				//correct direction. easy.
				//append b to newwalk. //i is already in the correct position.
				j=1;	//Note: skip overlapping first edge.
				while(built_walks_direcs2[b][j]!=0){
					newwalk[i] = built_walks_direcs2[b][j];
					i++;
					j++;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = built_walks_end2[b][i];
				}
			}
			else{
				//need to append reversed version of b
				j=0;
				while(built_walks_direcs2[b][j+1]!=0){
					j++;
				}
				j--; //Note: skip overlapping first edge.
				while(j>=0){
					newwalk[i] = reverse_direc(built_walks_direcs2[b][j]);
					i++;
					j--;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = built_walks_start2[b][i];
				}
			}

			//replace a with newwalk.
			for(i=0; i<=2; i++){
				built_walks_start2[a][i] = newwalk_start[i];
				built_walks_end2[a][i] = newwalk_end[i];
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				built_walks_direcs2[a][i] = newwalk[i];
			}

			//move walk b to end of array
			while(b<num_built_walks2-1){
				for(i=0; i<=2; i++){
					built_walks_start2[b][i] = built_walks_start2[b+1][i];
					built_walks_end2[b][i] = built_walks_end2[b+1][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					built_walks_direcs2[b][i] = built_walks_direcs2[b+1][i];
				}
				b++;
			}

			//delete b
			for(i=0; i<=2; i++){
				built_walks_start2[b][i] = -1;
				built_walks_end2[b][i] = -1;
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				built_walks_direcs2[b][i] = 0;
			}
			num_built_walks2--;
		}
		else{	//close up SAP.
/*
printf("\ntime to close up the SAP:\n");
printbuiltwalks_all();
printf("walktoadd:");
i=0;
while(walktoadd[i]!=0){
	printf(" %d", walktoadd[i]);
	i++;
}
printf("\n");
*/			//either append walktoadd in same order or reverse
			i=0;
			while(built_walks_direcs2[0][i]!=0){
				i++;
			}
			if(built_walks_end2[0][1]==walktoadd_start[1] && built_walks_end2[0][2]==walktoadd_start[2]){
				//just append
				j=1; //Note: skip overlapping first edge (which will be a 1).
				while(walktoadd[j]!=2){	//also skip overlapping last edge (which will be a 2)
					built_walks_direcs2[0][i]=walktoadd[j];
					i++;
					j++;
				}
			}
			else{
				//append reverse of walktoadd
				j=0;
				while(walktoadd[j+1]!=2){	//Note: skip overlapping first edge (which will be a 2)
					j++;
				}
				while(j>=1){ //also skip overlapping last edge (which will be a 1)
					built_walks_direcs2[0][i] = reverse_direc(walktoadd[j]);
					i++;
					j--;
				}
			}
			//set actual startpoint of SAP (and equivalent endpoint)
			built_walks_start2[0][0]=totalspan;
			for(i=0; i<=2; i++){
				built_walks_end2[0][i]=built_walks_start2[0][i];
				built_walks_end2[0][i]=built_walks_start2[0][i];
			}
		}
/*
printf("%d-th part of endhinge added:\n", walknum+1);
printbuiltwalks_all();
*/
	}
}


int reverse_direc(int direc){
	if(direc==1){
		return 2;
	}
	else if(direc==2){
		return 1;
	}
	else if(direc==3){
		return 4;
	}
	else if(direc==4){
		return 3;
	}
	else if(direc==5){
		return 6;
	}

	else if(direc==6){
		return 5;
	}
	else{
		printf("error in reverse_direc(). Exitting\n");
		exit(1);
	}
}


void printtofile(){
	if(filetotal>=maxpolys){
		//need to start new file
		printf("Finished filling file %lu. It contains %lu 2SAPs. Creating a new one.\n", filenum, filetotal);
		fprintf(fp, "-999\n");
		fclose(fp);
		filenum++;
		sprintf(filename, "MC2SAPsHamL%dM%dspan%drun%dnum%lu.txt", L, M, totalspan, runnum, filenum);
		fp = fopen(filename, "w");	//create or overwrite the file "filename

		if(fp != NULL){
			fprintf(fp, "UofS\n");	//first line in file is always "UofS"
			printf("printed UofS in file '%s'\n", filename);
		}
		else{
			printf("file pointer is pointing to NULL\n");
			exit(1);
		}
		filetotal=0;
		
	}
	//record
	fprintf(fp, "%d %d %d\n", built_walks_start[0][0], built_walks_start[0][1], built_walks_start[0][2]);	//starting point
	int i=0;
	while(built_walks_direcs[0][i] != 0){
		fprintf(fp, "%d\n", built_walks_direcs[0][i]);
		i++;
	}
	fprintf(fp, "-111\n");

	fprintf(fp, "%d %d %d\n", built_walks_start2[0][0], built_walks_start2[0][1], built_walks_start2[0][2]);	//starting point
	i=0;
	while(built_walks_direcs2[0][i] != 0){
		fprintf(fp, "%d\n", built_walks_direcs2[0][i]);
		i++;
	}
	fprintf(fp, "-111\n");
	filetotal++;
}



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
#else
#include "../topology/LFlag_norder2.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the left and returns 0 otherwise */
/* it also uses the global variable ordertemplate */
#include "../topology/LFlag_endhinge_norder2.c"

//#include "../topology/RFlag.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the right and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#if vM*vL<7
#include "../sections/Num_section_6V.c"  /*unsigned long int num_section(int side)*/
/* This function takes the arguement (int side) and assigns a unique number */
/* to the section on side = side */
#include "../sections/Num_section_6V_nonordered.c"
#include "../sections/Num_section_6V_nonordered2.c"
#include "../sections/Num_section_6V_endhinge_nonordered.c"
#include "../sections/Num_section_6V_endhinge_nonordered2.c"
#elif vM*vL<9
#include "../sections/Num_section_8V.c"
#include "../sections/Num_section_8V_nonordered.c"
#include "../sections/Num_section_8V_nonordered2.c"
#include "../sections/Num_section_8V_endhinge_nonordered.c"
#include "../sections/Num_section_8V_endhinge_nonordered2.c"
#elif vM*vL<11
#include "../sections/Num_section_10V.c"
#include "../sections/Num_section_10V_nonordered.c"
#include "../sections/Num_section_10V_nonordered2.c"
#include "../sections/Num_section_10V_endhinge_nonordered.c"
#include "../sections/Num_section_10V_endhinge_nonordered2.c"
#elif vM*vL<13
#include "../sections/Num_section_12V.c"
#include "../sections/Num_section_12V_nonordered.c"
#include "../sections/Num_section_12V_nonordered2.c"
#include "../sections/Num_section_12V_endhinge_nonordered.c"
#include "Num_section_12V_endhinge_nonordered2.c"
#endif

#endif



#include "../sections/printsection.c"

#include "../utils/int_vecalloc.c"
#include "../utils/unsgn_vecalloc.c"

#include "../utils/unsgnlong_vecalloc.c" /*unsigned long int unsgnlong_vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with unsigned long integer enteries */

#include "../utils/vecalloc.c"		/*vec_ent * vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with enteries of type vec_ent */

#include "../utils/matalloc.c"		/*mat_ent **matalloc(int rowlow, int rowhigh, int collow, int colhigh)*/
/* used for dynamically allocating memory for matrices with enteries of type mat_ent */





	









