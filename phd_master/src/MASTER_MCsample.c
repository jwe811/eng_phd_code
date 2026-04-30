//This program will sample SAPs of a certain span uniformly at random (by span).
//uses TMcalc_PrintEvectors.c, which prints out eigenvectors to:
//L_Evector_TS_L%dM%d.txt
//R_Evector_TS_L%dM%d.txt
//These files must be in same folder.

//compiled with: gcc -lm -Wall -O3 -o mc_master MASTER_MCsample.c
//ran using ./MCsample.o

//sample SAPs will go to files: MCpolysL%dM%dspan%dnum%d.txt

#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <unistd.h>
#include <getopt.h>
#include "../include/marsaglia.h"


/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/
// These are now runtime parameters
int L = 2;
int M = 1;
int totalspan = 100;
int samplesize = 100;
int runnum = 1;
unsigned int seednum = 227001;

int max_sections;
int max_tspans;
double dom_evalue;

int maxpolys = 10000;	// number of SAPs per file.

/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/

// System variables set at runtime
int vec_length;
int vM;
int vL;
int CS_mode = 0; // Replacement for CS macro
int mode = 0;    // 0: SAP, 1: HamSAP, 2: 2SAP, 3: Ham2SAP
int ham_check = 0; // 1 if Hamiltonian mode is active
char* output_dir = "data/SAPs";
double fval = 0.0; // Force value for spectral solver

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
	unsigned long int	outorder;		// the section number of the righthand section
	struct hinge_span	*nexthinge;	// pointer to another hinge_span data structure
	unsigned int		*hedges;	//contains info about which hinge-edges are filled

	int		*start[3];	//start[0][i] contains the starting x-coordinate of the i'th walk.
						//start[1][i] contains the starting y-coordinate of the i'th walk.
						//start[2][i] contains the starting z-coordinate of the i'th walk.
	int		*end[3];		//end[0][i] contains the ending x-coordinate of the i'th walk.
						//end[1][i] contains the ending y-coordinate of the i'th walk.
						//end[2][i] contains the ending z-coordinate of the i'th walk.
	int		**walks;	//walks[j][i] contains the j'th step in the i'th walk.	rest zero's.
	unsigned int	num_of_walks;	//number of walks in this 2-span.
	unsigned int	edgecount;	//number of edges in the hinge and first section

}	**first_hinge_span, **current_hinge_span; /* these are global pointers now */

	/* first_hinge_span is an array of pointers to hinge_span data structures where each pointer	*/
	/* in the array points to the first element of a linked list of hinge_span data structures	*/
	/* The array index implicely defines the "inorder" (or the lefthand section of the two_span)  	*/

	/* current_hinge_span is a dummy pointer used for traversing the linked lists mentioned above	*/
	/* an array of dummy pointers is needed because there is more than one linked list to consider	*/


struct endhinge {
	unsigned int		*hedges;	//contains info about which hinge-edges are filled
	int		*start[3];	//start[0][i] contains the starting x-coordinate of the i'th walk.
						//start[1][i] contains the starting y-coordinate of the i'th walk.
						//start[2][i] contains the starting z-coordinate of the i'th walk.
	int		*end[3];		//end[0][i] contains the ending x-coordinate of the i'th walk.
						//end[1][i] contains the ending y-coordinate of the i'th walk.
						//end[2][i] contains the ending z-coordinate of the i'th walk.
	int		**walks;	//walks[j][i] contains the j'th step in the i'th walk.	rest zero's.
	unsigned int	num_of_walks;	//number of walks in this endhinge
	unsigned int	edgecount;	//number of edges in the endhinge section
	unsigned short int side;		//0 if a left endhinge, 1 if a right endhinge
	struct endhinge* nextendhinge;// pointer to the next endhinge structure in the linked list
};

struct endhinge **firstendhinge;
struct endhinge **currentendhinge;



/***************************************************************************/
/****************** end of structures **************************************/
/***************************************************************************/


/***************************************************************************/
/*************** Global variables used in this program *********************/
/***************************************************************************/

unsigned int		***ordertemplate;	/* The first index is for side (left = 0 ; right = 1) */
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

////////////////////// arrays which hold the tspan info ////////////////////////////////

unsigned long int	*num_outsections;	/* For a given section (which has a number section_num assigned to it)*/
										/* num_outsections[section_num] is the number of two spans with this  */
										/* section as a first section                                         */
 
unsigned long int	**t_outsection;	/* An array of pointers to integers 					*/
										/* records the "outsection" of the two-span				*/
										/* t_outsection[insection][nth_tspan]				*/
										/* nth_tspan = nth two-span with firstsection = insection	*/

unsigned long int	**t_nrr;		/*An array of pointers to integers*/
										/*records the tspan "number"      */
										/* defined by the array index */
										/*example:		*/
										/*t_nrr[1][1] gives the two-span number */
										/* of the two-span with first section = 1 */
										/* and the second section is the first such */
										/* section which connects with first section = 1 */

unsigned int		**t_num_walks;	//t_num_walks[insection][nth_tspan] holds the number of walks in this tspan

int				****t_start;	//t_start[insection][nth_tspan][i][nth_walk] holds the i-coordinate (0=x, 1=y, 2=z) of the start point of the nth_walk

int				****t_end;	//t_end[insection][nth_tspan][i][nth_walk] holds the i-coordinate (0=x, 1=y, 2=z) of the end point of the nth_walk

int				****t_walks;	//t_end[insection][nth_tspan][direc][nth_walk] holds the "direc'th" direction the nth_walk
unsigned int		**tspans_edges;	//tspans_edges[insection][nth_tspan] holds the edgecount of the tspan

double			*L_Evector[2];
double			*R_Evector_solve[2]; // Renamed to avoid conflict with main's local R_Evector



////////////////////// arrays which hold the endhinge info ////////////////////////////////////

unsigned long int	*num_left_endhinges;

unsigned int		**Lend_num_walks;

int				****Lend_start;

int				****Lend_end;

int				****Lend_walks;


unsigned long int	*num_right_endhinges;

unsigned int		**Rend_num_walks;

int				****Rend_start;

int				****Rend_end;

int				****Rend_walks;


////////////////

unsigned int ***reordertemplate;	//holds the re-ordered ordertemplate to reduce number of unique sections.

unsigned short int **hingestatus;		//keeps track of which vertices in the hinge are occupied. 0=empty, 1=filled.
unsigned short int **alreadyentered;	//keeps track of where we've already initially entered. 0=no, 1=yes.

unsigned short int **colhingeedges; //keeps track of which column (i/M) edges in the hinge are occupied. 0=empty, 1=filled.
unsigned short int **rowhingeedges; //keeps track of which row (j/L) edges in the hinge are occupied. 0=empty, 1=filled.

unsigned long int *sectionkey;	//will hold valid (re-ordered) section nums. index is the key num
									//note: sectionkey[0] is garbage. Key numbers start at 1.

unsigned long int dupcounter=0;	//will count the number of duplicate 2-spans created and not stored.

int **curstart;		//will hold the start-points of the walks being built
int **curend;		//will hold the end-points of the walks being built
int **curwalks;	//will hold the walks as they are being built.
						//curwalks[j][i] contains the j'th step in the i'th walk.	rest zero's.

unsigned int num_walks=0;	//number of walks in the 2-span

//endhinge stuff
int EndOrdNum[2];		//contains info on # of edges in left section and endhinge

int **endtemplate;		//contains the section associated with the end hinge
int **rendtemplate;	//holds the re-ordered endtemplate to reduce number of unique sections.
int **endtemplate2;
int **rendtemplate2;

unsigned long int tot_right_endhinges=0;
unsigned long int tot_left_endhinges=0;
unsigned long int num_duplicate_right_endhinges=0;
unsigned long int num_duplicate_left_endhinges=0;
unsigned long int num_tspans=0;
unsigned long int num_duplicate_tspans=0;

//integer division rounds down in case of vM*vL is odd
int **built_walks_start;
int **built_walks_end;
int **built_walks_direcs;

int num_built_walks;

char 			filename[100];	//this will hold the filename of the text file that we want to write to
FILE*			fp;			//file pointer. will point to the file that is to be written to.
unsigned long int filetotal=0;	//number of polys in the file 
unsigned long int filenum=1;		//number of polygon files created.





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

void            	enterhinge(int, int, int, int (*)[3], int);
void            	leavehinge(int, int, int, int (*)[3], int);
void            	rowedges(int, int, int (*)[3], int);
void            	coledges(int, int, int (*)[3], int);

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

/****  Memory allocation functions  ********************************************************/

struct hinge_span	*newhinge(unsigned int);
struct endhinge	*newendhinge(unsigned int);
int				*int_vecalloc(int low, int high);
unsigned int		*unsgn_vecalloc(int low, int high);
unsigned long int	*unsgnlong_vecalloc(int low, int high);
vec_ent			*vecalloc(int low, int high);
mat_ent			**matalloc(int rowlow, int highrow, int lowcol, int highcol);


/**** Added Functions  ********************************************************************/
void				findsection(double, double, double);
void				printsection(int);

void 			fillreordertemplate(int, int);	//fills reordertemplate appropriately based on ordertemplate
void				printordtemp();
void				printreordtemp();

void enterendhinge(int i, int j, int curlength);
void leaveendhinge(int i, int j, int curlength);
void endrowedges(int i, int j, int curlength);
void endcoledges(int i, int j, int curlength);

unsigned short int LFlag_endhinge(int (*EndOrdNum)[2]);
void recordendtemplate();
void fillrendtemplate();

void fillendtemplate2();
void fillrendtemplate2();

unsigned long int num_section_endhinge_norder(void);
unsigned long int num_section_endhinge_norder2(void);

void printbuiltwalks();
void printbuiltwalks_all();

void add_to_built_walks(unsigned long int secnum, int nth_tspan);
void add_right_endhinge(unsigned long int secnum, int nth_endhinge);

int reverse_direc(int direc);

void printtofile();

void generate_evectors() {
	printf("Starting integrated spectral solver...\n");
	double fugacity = 1.0; // Uniform sampling uses eigenvectors at z=1.0
	
	// The integrated power method uses L_Evector and R_Evector_solve
	// We need to provide them as expected by pw_meth_ts_LRvec_fcheck.c
	extern double max_eval_LRvec(double fugacity);
	double calculated_evalue = max_eval_LRvec(fugacity) + 1.0;
	printf("Calculated dominant eigenvalue: %.15f (Expected: %.15f)\n", calculated_evalue, dom_evalue);
	
	// VERIFICATION: Check against archival file if it exists
	char r_filename[100];
	if (ham_check) {
		sprintf(r_filename, "data/R_EvectorHam_TS_L%dM%d.txt", L, M);
	} else {
		sprintf(r_filename, "data/R_Evector_TS_L%dM%d.txt", L, M);
	}

	FILE *r_fp = fopen(r_filename, "r");
	if (r_fp != NULL) {
		printf("Verifying calculated eigenvectors against archival data in %s...\n", r_filename);
		double max_diff = 0.0;
		int diff_count = 0;
		for (int i = 1; i <= (int)num_tspans; i++) {
			double file_val;
			if (fscanf(r_fp, "%lf", &file_val) == 1) {
				double diff = fabs(R_Evector_solve[0][i] - file_val);
				if (diff > max_diff) max_diff = diff;
				if (diff > 1e-9) diff_count++;
			}
		}
		printf("Verification complete. Max difference: %e, Significant differences (>1e-9): %d\n", max_diff, diff_count);
		fclose(r_fp);
	} else {
		printf("No archival file found for verification (%s). Proceeding with calculated values.\n", r_filename);
	}
}

/***************************************************************************/
void set_system_params() {
    vM = M + 1;
    vL = L + 1;
    if (M == 0) CS_mode = 1;
    
    if (ham_check) {
        if (M == 1 && L == 1) {
            max_sections = 8;
            max_tspans = 22;
            dom_evalue = 3.732050810014727;
        } else if (M == 1 && L == 2) {
            max_sections = 73;
            max_tspans = 649;
            dom_evalue = 14.076438172824952;
        } else if (M == 1 && L == 3) {
            max_sections = 742;
            max_tspans = 19554;
            dom_evalue = 49.643407510907970;
        } else if (M == 1 && L == 4) {
            max_sections = 9309;
            max_tspans = 728925;
            dom_evalue = 172.714480164060546;
        } else if (M == 1 && L == 5) {
            max_sections = 138038;
            max_tspans = 32294131;
            dom_evalue = 596.477267188525502;
        } else if (M == 2 && L == 2) {
            max_sections = 2619;
            max_tspans = 184574;
            dom_evalue = 104.488979151088131;
        } else if (M == 2 && L == 3) {
            max_sections = 138322;
            max_tspans = 39671908;
            dom_evalue = 880.559222610926781;
        } else {
            fprintf(stderr, "Unsupported Hamiltonian L and M values (%d, %d).\n", L, M);
            exit(1);
        }
    } else {
        if (M == 1 && L == 1) {
            max_sections = 8;
            max_tspans = 48;
            dom_evalue = 7.0;
        } else if (M == 1 && L == 2) {
            max_sections = 73;
            max_tspans = 1829;
            dom_evalue = 34.360180657516501;
        } else if (M == 1 && L == 3) {
            max_sections = 742;
            max_tspans = 70306;
            dom_evalue = 165.169003076277392;
        } else if (M == 1 && L == 4) {
            max_sections = 9309;
            max_tspans = 3165653;
            dom_evalue = 770.927136914359608;
        } else if (M == 1 && L == 5) {
            max_sections = 138038;
            max_tspans = 165637127;
            dom_evalue = 3562.911364811740896;
        } else if (M == 2 && L == 2) {
            max_sections = 2619;
            max_tspans = 513585;
            dom_evalue = 416.870158713340345;
        } else if (M == 2 && L == 3) {
            max_sections = 138322;
            max_tspans = 201423784;
            dom_evalue = 4945.864451078689854;
        } else {
            fprintf(stderr, "Unsupported Standard L and M values (%d, %d).\n", L, M);
            exit(1);
        }
    }
    vec_length = max_sections + 1;
}

unsigned int*** allocate_3d_uint(int d1, int d2, int d3) {
    unsigned int ***arr = (unsigned int***)malloc(d1 * sizeof(unsigned int**));
    for (int i = 0; i < d1; i++) {
        arr[i] = (unsigned int**)malloc(d2 * sizeof(unsigned int*));
        for (int j = 0; j < d2; j++) {
            arr[i][j] = (unsigned int*)calloc(d3, sizeof(unsigned int));
        }
    }
    return arr;
}

unsigned short int** allocate_2d_ushort(int d1, int d2) {
    unsigned short int **arr = (unsigned short int**)malloc(d1 * sizeof(unsigned short int*));
    for (int i = 0; i < d1; i++) {
        arr[i] = (unsigned short int*)calloc(d2, sizeof(unsigned short int));
    }
    return arr;
}

int** allocate_2d_int(int d1, int d2) {
    int **arr = (int**)malloc(d1 * sizeof(int*));
    for (int i = 0; i < d1; i++) {
        arr[i] = (int*)calloc(d2, sizeof(int));
    }
    return arr;
}

void allocate_globals() {
    ordertemplate = allocate_3d_uint(2, vM, vL);
    reordertemplate = allocate_3d_uint(2, vM, vL);
    hingestatus = allocate_2d_ushort(vM, vL);
    alreadyentered = allocate_2d_ushort(vM, vL);
    colhingeedges = allocate_2d_ushort(vM, vL); // Used as colhingeedges[M][vL], vM is enough
    rowhingeedges = allocate_2d_ushort(vM, vL);
    sectionkey = (unsigned long int*)calloc(vec_length, sizeof(unsigned long int));
    
    num_outsections = (unsigned long int*)calloc(vec_length, sizeof(unsigned long int));
    t_outsection = (unsigned long int**)calloc(vec_length, sizeof(unsigned long int*));
    t_nrr = (unsigned long int**)calloc(vec_length, sizeof(unsigned long int*));
    t_num_walks = (unsigned int**)calloc(vec_length, sizeof(unsigned int*));
    t_start = (int****)calloc(vec_length, sizeof(int***));
    t_end = (int****)calloc(vec_length, sizeof(int***));
    t_walks = (int****)calloc(vec_length, sizeof(int***));
    tspans_edges = (unsigned int**)calloc(vec_length, sizeof(unsigned int*));
    
    L_Evector[0] = (double*)calloc(max_tspans + 1, sizeof(double));
    L_Evector[1] = (double*)calloc(max_tspans + 1, sizeof(double));
    R_Evector_solve[0] = (double*)calloc(max_tspans + 1, sizeof(double));
    R_Evector_solve[1] = (double*)calloc(max_tspans + 1, sizeof(double));
    
    num_left_endhinges = (unsigned long int*)calloc(vec_length, sizeof(unsigned long int));
    Lend_num_walks = (unsigned int**)calloc(vec_length, sizeof(unsigned int*));
    Lend_start = (int****)calloc(vec_length, sizeof(int***));
    Lend_end = (int****)calloc(vec_length, sizeof(int***));
    Lend_walks = (int****)calloc(vec_length, sizeof(int***));
    
    num_right_endhinges = (unsigned long int*)calloc(vec_length, sizeof(unsigned long int));
    Rend_num_walks = (unsigned int**)calloc(vec_length, sizeof(unsigned int*));
    Rend_start = (int****)calloc(vec_length, sizeof(int***));
    Rend_end = (int****)calloc(vec_length, sizeof(int***));
    Rend_walks = (int****)calloc(vec_length, sizeof(int***));
    
    first_hinge_span = (struct hinge_span**)calloc(vec_length, sizeof(struct hinge_span*));
    current_hinge_span = (struct hinge_span**)calloc(vec_length, sizeof(struct hinge_span*));
    firstendhinge = (struct endhinge**)calloc(vec_length, sizeof(struct endhinge*));
    currentendhinge = (struct endhinge**)calloc(vec_length, sizeof(struct endhinge*));
    
    curstart = allocate_2d_int(3, vM * vL);
    curend = allocate_2d_int(3, vM * vL);
    curwalks = allocate_2d_int(vM * vL + 1, vM * vL);
    
    endtemplate = allocate_2d_int(vM, vL);
    rendtemplate = allocate_2d_int(vM, vL);
    endtemplate2 = allocate_2d_int(vM, vL);
    rendtemplate2 = allocate_2d_int(vM, vL);
    
    built_walks_start = allocate_2d_int(vM * vL / 2 + 1, 3);
    built_walks_end = allocate_2d_int(vM * vL / 2 + 1, 3);
    built_walks_direcs = allocate_2d_int(vM * vL / 2 + 1, vM * vL * (totalspan + 1) + 1);
}

void check_reachability(void);

/******************* end of functions used in this program *****************/
/***************************************************************************/


int main(int argc, char *argv[]) {

	int opt;
	while ((opt = getopt(argc, argv, "L:M:s:n:r:S:m:")) != -1) {
		switch (opt) {
			case 'L': L = atoi(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 's': totalspan = atoi(optarg); break;
			case 'n': samplesize = atoi(optarg); break;
			case 'r': runnum = atoi(optarg); break;
			case 'S': seednum = (unsigned int)atoi(optarg); break;
			case 'm': 
				mode = atoi(optarg);
				if (mode == 1 || mode == 3) ham_check = 1;
				break;
			default:
				fprintf(stderr, "Usage: %s [-L L] [-M M] [-s totalspan] [-n samplesize] [-r runnum] [-S seednum] [-m mode]\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	set_system_params();
	allocate_globals();

	if (mode == 1) output_dir = "data/HamSAPs";
	else if (mode == 2) output_dir = "data/2SAPs";
	else if (mode == 3) output_dir = "data/Ham2SAPs";
	else output_dir = "data/SAPs";

	unsigned int seed=seednum;
	initran_(&seed);

/**** variable used in main ************************************************/

	/*int	           	(*pointordNum)[3]; *//*pointer to an array of int*/
	int				ordNum[3],side = 0;
	int             	i,j,k;	   /*used in for loops */

/**** functions called by main *********************************************/

	clock_t         clock(void);
	void            conv_to_array(void);
	void            conv_endhinges_to_array(void);

/**************************************************************************/
/**** start of main program ***********************************************/
/**************************************************************************/

	if(system("clear")){} /* clears the screen */
	if(system("date")){}  /* prints the date and time */

	printf("section-section version\n");

	printf("\nThe variables are as follows M: %i L: %i\n", M, L);


	clock();	/* Start clock to time program */

	ordNum[0] = 1;	/* first edge in section 0 gets numbered 1 (left side of 2-span)	*/
	ordNum[1] = 1;	/* first edge in section 1 gets numbered 1 (right side of 2-span)	*/
	ordNum[2] = 0;	/* to start there are no edges in the 2-span (edges in hinge)		*/

	for (i=1;i<=max_sections;i++){
		current_hinge_span[i] = newhinge(1);		/* free up space of hinge structure */
		first_hinge_span[i] = current_hinge_span[i];	/* make first_hinge_span[i] point to this first one */
		currentendhinge[i] = newendhinge(1);		/* free up space of endhinge structure */
		firstendhinge[i] = currentendhinge[i];		/* make firstendhinge[i] point to this first one */
		sectionkey[i]=0;
	}
	printf("space created for hinge structures, endhinge structures, and sectionkey initialized\n");

	//initializations//////////////////////////////////////////////////
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

	num_walks=0;
	for(i=0; i<=vM*vL-1; i++){
		for(j=0; j<=2; j++){
			curstart[j][i]=-1;
			curend[j][i]=-1;
		}
		for(j=0; j<=vM*vL; j++){
			curwalks[j][i]=0;
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

	printf("%lu 2-spans recorded\n", num_tspans);
	printf("%lu duplicate 2spans\n", num_duplicate_tspans);


	//generate the endhinges
	//we are going to create "right endhinges" (entering from the left, can't leave to the right)
	//When we record these right endhinges, then we will make appropriate changes to also record it's "mirror" left endhinge

	//initialize stuff to zero.
	num_walks=0;
	for(i=0; i<=vM*vL-1; i++){
		for(j=0; j<=2; j++){
			curstart[j][i]=-1;
			curend[j][i]=-1;
		}
		for(j=0; j<=vM*vL; j++){
			curwalks[j][i]=0;
		}
	}

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

	EndOrdNum[0]=1;	//first edge in section will get the label 1
	EndOrdNum[1]=0;	//no edges currently in the endhinge

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

	printf("%lu right endhinges recorded\n", tot_right_endhinges);
	printf("%lu left endhinges recorded\n", tot_left_endhinges);
	printf("%lu duplicate right endhinges\n", num_duplicate_right_endhinges);
	printf("%lu duplicate left endhinges\n", num_duplicate_left_endhinges);

//tests
/*
	int k;
	for(i=1; i<=max_sections; i++){
		currentendhinge[i]=firstendhinge[i];
		while(currentendhinge[i]->nextendhinge!=NULL){
			currentendhinge[i]=currentendhinge[i]->nextendhinge;
			if(currentendhinge[i]->side==0){
				printf("Endhinge with sectionkey =%d, side=%d, num_walks=%d\n", i, currentendhinge[i]->side, currentendhinge[i]->num_of_walks);
				for(j=0; j<=currentendhinge[i]->num_of_walks-1; j++){
					printf("walk %d. Start=(%d,%d,%d), Directions: ", j, currentendhinge[i]->start[0][j], currentendhinge[i]->start[1][j], currentendhinge[i]->start[2][j]);
					for(k=0; k<=vM*vL; k++){
						printf("%d ", currentendhinge[i]->walks[k][j]);
					}
					printf(". End=(%d,%d,%d)\n", currentendhinge[i]->end[0][j], currentendhinge[i]->end[1][j], currentendhinge[i]->end[2][j]);
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

	unsigned long int totalsecs = 0;
	for (i = 1; i <= (unsigned long int)max_sections; i++) {
		if (num_outsections[i] > 0) {
			totalsecs++;
		}
	}
	printf("Final: Total sections=%lu, total 2spans=%lu\n", totalsecs, num_tspans);
	printf("max_sections was set = %d\n", max_sections);
	printf("max_tspans was set = %d\n", max_tspans);
	printf("There are a total of %lu potential two-span patterns.\n", num_2_spans - 1);
	printf("There are a total of %lu valid two-span patterns.\n", num_tspans);
	printf("Number of duplicate 2-spans generated that were not recorded=%lu\n", dupcounter);
    printf("Dominant Eigenvalue: %.15f\n", dom_evalue);

	if (ham_check) {
		check_reachability();
	}


/*
	int k,l;
	for(i=1; i<=max_sections; i++){
		for(j=1; j<=num_outsections[i]; j++){
			printf("Section %d connects to section %lu. It consists of %d walks:\n", i, t_outsection[i][j], t_num_walks[i][j]);
			for(k=0; k<=t_num_walks[i][j]-1; k++){
				printf("Walk %d. Start (%d,%d,%d). Directions: ", k, t_start[i][j][0][k], t_start[i][j][1][k], t_start[i][j][2][k]);
				for(l=0; l<=vM*vL; l++){
					printf("%d ", t_walks[i][j][l][k]);
				}
				printf(". End (%d,%d,%d)\n", t_end[i][j][0][k], t_end[i][j][1][k], t_end[i][j][2][k]);
			}
		}
	}
*/



	printf("Finished generating all 2-spans and endhinges.\n");

	if(totalspan<2){
		printf("totalspan=%d. This program only works for span>=2. Exitting\n", totalspan);
		exit(1);
	}

	printf("\nNOW SAMPLING: %d samples from L=%d, M=%d, span=%d%s\n", samplesize, L, M, totalspan, (ham_check ? " (Hamiltonian)" : ""));

	const char *file_prefix = (ham_check) ? "MCpolysHam" : "MCpolys";
	sprintf(filename, "%s/%sL%dM%dspan%drun%dnum%lu.txt", output_dir, file_prefix, L, M, totalspan, runnum, filenum);
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

	// Always calculate eigenvectors at runtime
	generate_evectors();

	double* R_Evector;
	R_Evector = R_Evector_solve[0]; // Point directly to the calculated global vector
/*
	for(i=1; i<= max_tspans; i++){
		printf("L_Evector[%d] = %.15f\n", i, L_Evector[i]);
		printf("R_Evector[%d] = %.15f\n", i, R_Evector[i]);
	}

	exit(1);
*/

	double maxt_one=0.0;	//used in rejection
	double t_one[tot_left_endhinges];	//each left end hinge has a t_one value.
	curLEH=0;
	for(i=1; i<=max_sections; i++){
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


	for(i=1; i<=max_sections; i++){
		for(j=1; j<=num_outsections[i]; j++){
			t_two[t_nrr[i][j]] = num_right_endhinges[t_outsection[i][j]] / R_Evector[t_nrr[i][j]];
			if(t_two[t_nrr[i][j]] > maxt_two){
				maxt_two = t_two[t_nrr[i][j]];
			}
//			printf("filled t_two[%d]=%f. maxt_two=%f\n", t_nrr[i][j], t_two[t_nrr[i][j]], maxt_two);
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
			}
			for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
				built_walks_direcs[i][j]=0;
			}
		}

		chosenLEH = floor(ran1real_()*tot_left_endhinges); //random number [0,tot_left_endhinges-1]

//	printf("\ntot_left_endhinges=%d\n", tot_left_endhinges);
//	chance = (double)1/tot_left_endhinges;
//	printf("chose left endhinge.\n");
//	printf("numerator=%d. denominator=%d\n", 1, tot_left_endhinges);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);

//		printf("Sample=%lu. chosenLEH=%d, from range [%d,%lu].\n", curSample, chosenLEH, 0, tot_left_endhinges-1);
		curLEH=0;
		for(secnum=1; secnum<=max_sections; secnum++){
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
						//built_walks is filled with a left endhinge.
						curspan=1;
	/*					printf("\nLeft endhinge has been added to built_walks:\n");
						printbuiltwalks_all();
						printf("\n");
	*/
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
//	printf("first block chosen with prob=%f. It was the %d-th choice.\n", R_Evector[t_nrr[sec1][nth_tspan]] / t_one[curLEH], nth_tspan);
//	printf("numerator=%f. denominator=%f\n", R_Evector[t_nrr[sec1][nth_tspan]], t_one[curLEH]);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);
//						printf("first tspan chosen uniformly: section=%lu, nth_tspan=%d (possible [%d,%lu]), outsection=%lu\n", sec1, nth_tspan, 1, num_outsections[sec1], t_outsection[sec1][nth_tspan]);
						add_to_built_walks(sec1, nth_tspan);	//add nth_tspan from sec1 to built_walks
//						printf("tspan has added, here is the result:\n");
//						printbuiltwalks_all();
//						printf("\n");
						sec2 = t_outsection[sec1][nth_tspan];	//set sec1
						cur_tspan_num = t_nrr[sec1][nth_tspan];
						curspan++;

						while(curspan<totalspan){
//							printf("adding another tspan to section %lu. Probabilities:\n", sec2);
							sumofprobs=0.0;
							for(i=1; i<=num_outsections[sec2]; i++){
								sumofprobs = sumofprobs + R_Evector[t_nrr[sec2][i]] / R_Evector[cur_tspan_num] / dom_evalue;
//								printf("sec %lu -> sec %lu. tspan_num=%lu. probability=%f. cumulative prob=%f\n", sec2, t_outsection[sec2][i], t_nrr[sec2][i], R_Evector[t_nrr[sec2][i]] / R_Evector[cur_tspan_num] / dom_evalue, sumofprobs);
							}

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
//	printf("next block chosen with prob=%f.\n", R_Evector[t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue);
//	printf("numerator=%f. denominator=%f. dom_evalue=%f.\n", R_Evector[t_nrr[sec2][nth_tspan]], R_Evector[cur_tspan_num], dom_evalue);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);
//							printf("sec %lu -> sec %lu. tspan_num=%lu. probability=%f. cumulative prob=%f\n", sec2, t_outsection[sec2][nth_tspan], t_nrr[sec2][nth_tspan], R_Evector[t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue, sumofprobs);
							sec1=sec2;
							add_to_built_walks(sec1, nth_tspan);	//add nth_tspan from sec2 to built_walks
	//						printf("tspan has added, here is the result:\n");
	//						printbuiltwalks_all();
	//						printf("\n");
							sec2 = t_outsection[sec1][nth_tspan];
							cur_tspan_num = t_nrr[sec1][nth_tspan];
							curspan++;
						}

						//CHECK FOR REJECTION
//						printf("accept prob_two=%f\n", t_two[cur_tspan_num] / maxt_two);
						if( ran1real_() < t_two[cur_tspan_num] / maxt_two ){	//do not reject
//	chance = chance * t_two[cur_tspan_num] / maxt_two;
//	printf("passed 2nd rejection with prob=%f.\n", t_two[cur_tspan_num] / maxt_two);
//	printf("numerator=%f. denominator=%f\n", t_two[cur_tspan_num], maxt_two);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);
							chosenREH = 1 + floor(ran1real_()*num_right_endhinges[sec2]); //rando number [1, num_right_endhinges[sec2]]
//	chance = chance * (double)1/num_right_endhinges[sec2];
//	printf("chose final endhinge with prob=%f.\n", (double)1/num_right_endhinges[sec2]);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);
	//						printf("closing off SAP with endhinge: section=%lu, nth_endhinge=%d (possible [%d,%lu]).\n", sec2, chosenREH, 1, num_right_endhinges[sec2]);
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



	if(system("date")){}  /* prints the date and time */
	printf("\nProgram Complete.\n\n");

	return 0;

}				/* end of main */
/***************************************************************************/


/***** start of function conv_to_array *************************************/
/* called only by main */


void conv_to_array(void)
{
	unsigned long int section_num;
	unsigned long int arr_ent;	//array entry
	struct hinge_span *hinge_to_free;

	for (section_num = 1; section_num <= max_sections; section_num++) {
		num_outsections[section_num] = 0;	// initialize to zero
		current_hinge_span[section_num] = first_hinge_span[section_num];

		while ((*current_hinge_span[section_num]).nexthinge != NULL) {
			current_hinge_span[section_num] = (*current_hinge_span[section_num]).nexthinge;
			num_outsections[section_num]++;
		}

		t_outsection[section_num] = unsgnlong_vecalloc(1, num_outsections[section_num]);
		t_num_walks[section_num] = unsgn_vecalloc(1, num_outsections[section_num]);
		tspans_edges[section_num] = (unsigned int*)calloc(num_outsections[section_num] + 1, sizeof(unsigned int));
		t_nrr[section_num] = unsgnlong_vecalloc(1, num_outsections[section_num]);/*allocate memory n=number*/
	}

	//fill t_outsection and t_num_walks
	for (section_num = 1; section_num <= max_sections; section_num++) {
		arr_ent = 1;
		current_hinge_span[section_num] = first_hinge_span[section_num];
		while ((*current_hinge_span[section_num]).nexthinge != NULL) {
			current_hinge_span[section_num] = (*current_hinge_span[section_num]).nexthinge;
			t_outsection[section_num][arr_ent] = current_hinge_span[section_num]->outorder;
			t_num_walks[section_num][arr_ent] = current_hinge_span[section_num]->num_of_walks;
			tspans_edges[section_num][arr_ent] = current_hinge_span[section_num]->edgecount;
			arr_ent++;
		}
	}



	//allocate space for t_start, t_end, and t_walks
	printf("allocating space for t_start, t_end, and t_walks\n");
	int i, j;
	for(section_num=1; section_num<=max_sections; section_num++){
		t_start[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_start[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		t_end[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_end[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		t_walks[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_walks[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		for (i = 1; i <= num_outsections[section_num]; i++){
			t_start[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (t_start[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			t_end[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (t_end[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			t_walks[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (t_walks[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			for (j = 0; j < 3; j++){
				t_start[section_num][i][j] = (int*)malloc(t_num_walks[section_num][i] * sizeof(int));
		   		if (t_start[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				t_end[section_num][i][j] = (int*)malloc(t_num_walks[section_num][i] * sizeof(int));
		   		if (t_end[section_num][i][j] == NULL) {
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
			}
		}
	}
	printf("Done allocating space for t_start, t_end, and t_walks\n");


	///Fill t_start, t_end, t_walks
	printf("Filling t_start, t_end, and t_walks\n");
	for (section_num = 1; section_num <= max_sections; section_num++) {
//		printf("section_num=%lu\n", section_num);
		arr_ent = 1;
		current_hinge_span[section_num] = first_hinge_span[section_num];
		while ((*current_hinge_span[section_num]).nexthinge != NULL) {
//			printf("arr_ent=%lu\n\n", arr_ent);
			hinge_to_free = current_hinge_span[section_num];
			current_hinge_span[section_num] = (*current_hinge_span[section_num]).nexthinge;
			free(hinge_to_free);
			for(i=0; i<=(current_hinge_span[section_num]->num_of_walks)-1; i++){
				for(j=0; j<=2; j++){
//					printf("recording %d-th coordinate\n", j);
//					printf("section_num=%lu, arr_ent=%lu, j=%d, i=%d\n", section_num, arr_ent, j, i);
					t_start[section_num][arr_ent][j][i] = current_hinge_span[section_num]->start[j][i];
//					printf("t_start[%lu][%lu][%d][%d] = %d\n", section_num, arr_ent, j, i, t_start[section_num][arr_ent][j][i]);
					t_end[section_num][arr_ent][j][i] = current_hinge_span[section_num]->end[j][i];
//					printf("t_end[%lu][%lu][%d][%d] = %d\n", section_num, arr_ent, j, i, t_end[section_num][arr_ent][j][i]);
				}
				for(j=0; j<=vM*vL; j++){
					t_walks[section_num][arr_ent][j][i] = current_hinge_span[section_num]->walks[j][i];
				}
			}
			arr_ent++;
		}
	}
	printf("Done filling t_start, t_end, and t_walks\n");

	//fill t_nrr
	unsigned long int tspan_num=0;
	unsigned long int nth_outsection;
	for (section_num=1; section_num<=max_sections;section_num++){
		for (nth_outsection=1;nth_outsection<=num_outsections[section_num];nth_outsection++){
			t_nrr[section_num][nth_outsection] = ++tspan_num;
		}
	}
			
	return;
}




/***** start of function conv_to_array *************************************/
/* called only by main */


void conv_endhinges_to_array(void)
{
	unsigned long int section_num;
	unsigned long int L_arr_ent;	//array entry
	unsigned long int R_arr_ent;	//array entry
	struct endhinge *hinge_to_free;

	for (section_num = 1; section_num <= max_sections; section_num++) {
		num_left_endhinges[section_num] = 0;	// initialize to zero
		num_right_endhinges[section_num] = 0;	// initialize to zero
		currentendhinge[section_num] = firstendhinge[section_num];

		while ((*currentendhinge[section_num]).nextendhinge != NULL) {
			currentendhinge[section_num] = (*currentendhinge[section_num]).nextendhinge;
			if(currentendhinge[section_num]->side==0){
				num_left_endhinges[section_num]++;
			}
			else{
				num_right_endhinges[section_num]++;
			}
		}

		Lend_num_walks[section_num] = unsgn_vecalloc(1, num_left_endhinges[section_num]);
		Rend_num_walks[section_num] = unsgn_vecalloc(1, num_right_endhinges[section_num]);
	}

	//allocate space for Lend_start, Lend_end, and Lend_walks
	int i, j;
	for(section_num=1; section_num<=max_sections; section_num++){
/////////
		Lend_start[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_start[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		Lend_end[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_end[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		Lend_walks[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_walks[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
/////////
		Rend_start[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_start[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		Rend_end[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_end[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}

		Rend_walks[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_walks[section_num] == NULL) {
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
			Lend_end[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Lend_end[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			Lend_walks[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (Lend_walks[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			for (j = 0; j < 3; j++){
				Lend_start[section_num][i][j] = (int*)malloc(Lend_num_walks[section_num][i] * sizeof(int));
		   		if (Lend_start[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				Lend_end[section_num][i][j] = (int*)malloc(Lend_num_walks[section_num][i] * sizeof(int));
		   		if (Lend_end[section_num][i][j] == NULL) {
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
		}
/////////
		for (i = 1; i <= num_right_endhinges[section_num]; i++){
			Rend_start[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Rend_start[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			Rend_end[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Rend_end[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
			Rend_walks[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (Rend_walks[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}

			for (j = 0; j < 3; j++){
				Rend_start[section_num][i][j] = (int*)malloc(Rend_num_walks[section_num][i] * sizeof(int));
		   		if (Rend_start[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(0);
				}
				Rend_end[section_num][i][j] = (int*)malloc(Rend_num_walks[section_num][i] * sizeof(int));
		   		if (Rend_end[section_num][i][j] == NULL) {
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
			}
		}
	}


	///Fill L/Rend_start, L/Rend_end, L/Rend_walks
	for (section_num = 1; section_num <= max_sections; section_num++) {
		L_arr_ent = 1;
		R_arr_ent = 1;
		currentendhinge[section_num] = firstendhinge[section_num];
		while ((*currentendhinge[section_num]).nextendhinge != NULL) {
			hinge_to_free = currentendhinge[section_num];
			currentendhinge[section_num] = (*currentendhinge[section_num]).nextendhinge;
			free(hinge_to_free);
			if(currentendhinge[section_num]->side==0){	//leftendhinge
				Lend_num_walks[section_num][L_arr_ent] = currentendhinge[section_num]->num_of_walks;
				for(i=0; i<=(currentendhinge[section_num]->num_of_walks)-1; i++){
					for(j=0; j<=2; j++){
						Lend_start[section_num][L_arr_ent][j][i] = currentendhinge[section_num]->start[j][i];
						Lend_end[section_num][L_arr_ent][j][i] = currentendhinge[section_num]->end[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						Lend_walks[section_num][L_arr_ent][j][i] = currentendhinge[section_num]->walks[j][i];
					}
				}
				L_arr_ent++;
			}
			else{ //rightendhinge
				Rend_num_walks[section_num][R_arr_ent] = currentendhinge[section_num]->num_of_walks;
				for(i=0; i<=(currentendhinge[section_num]->num_of_walks)-1; i++){
					for(j=0; j<=2; j++){
						Rend_start[section_num][R_arr_ent][j][i] = currentendhinge[section_num]->start[j][i];
						Rend_end[section_num][R_arr_ent][j][i] = currentendhinge[section_num]->end[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						Rend_walks[section_num][R_arr_ent][j][i] = currentendhinge[section_num]->walks[j][i];
					}
				}
				R_arr_ent++;
			}
		}
	}
			
	return;
}



/***************************************************************************/

void
enterhinge(int i, int j, int side, int (*pointordNum)[3], int curlength)
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
leavehinge(int i, int j, int side, int (*pointordNum)[3], int curlength)
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
		/* this may be a valid two span */
			if (LFlag(pointordNum) /*&& RFlag(pointordNum)*/){
				int isHam = 1;
				if (ham_check) {
					for (int ii = 0; ii <= M; ii++) {
						for (int jj = 0; jj <= L; jj++) {
							if (hingestatus[ii][jj] == 0) { isHam = 0; break; }
						}
						if (!isHam) break;
					}
				}

				if (isHam) {
					fillreordertemplate((*pointordNum)[0]-1, (*pointordNum)[1]-1);
					valid_2_spans++; /*This is a valid 2 span so it is counted. Includes duplicates. */
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

void rowedges(int i, int j, int (*pointordNum)[3], int curlength)  //y-direction
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

void coledges(int i, int j, int (*pointordNum)[3], int curlength)  //z-direction
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
recordtemplate( int (*pointordNum)[3])
/* This function records the pertinent information about a two-span */
/* which has been verified as being valid into a linked list */
{
	unsigned short int	Num_edges= (*pointordNum)[0] + (*pointordNum)[2] - 1;
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
		num_tspans++;

//		printf("actually recording\n");
		(*current_hinge_span[inNum]).nexthinge = newhinge(num_walks); /* add a newhinge to the linked list */
	
		current_hinge_span[inNum] = (*current_hinge_span[inNum]).nexthinge; /* let current_hinge_span[inNum] point to newhinge */
		(*current_hinge_span[inNum]).outorder = outNum;	/* record the ordering on the out side */
		(*current_hinge_span[inNum]).num_of_walks = num_walks;
		(*current_hinge_span[inNum]).edgecount = Num_edges;
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

/*		printf("Recorded 2-span. Section %d to %d, with hedges: ", inNum, outNum);
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

		if (LFlag_endhinge(&EndOrdNum)){
			int isHam = 1;
			if (ham_check) {
				for (int ii = 0; ii <= M; ii++) {
					for (int jj = 0; jj <= L; jj++) {
						if (hingestatus[ii][jj] == 0) { isHam = 0; break; }
					}
					if (!isHam) break;
				}
			}

			if (isHam) {
				fillrendtemplate();
				recordendtemplate(); /*the information contained in endtemplate is recorded*/
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



/***************************************************************************/
void
recordendtemplate()
{
	unsigned short int	Num_edges= EndOrdNum[0] + EndOrdNum[1] - 1;
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
//	printf("NEW secnum=%d\n", secnum);


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
		if(test->side==1){	//if right endhinge
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
//		printf("actually recording\n");
		tot_right_endhinges++;
		(*currentendhinge[secnum]).nextendhinge = newendhinge(num_walks); /* add a newendhinge to the linked list */
	
		currentendhinge[secnum] = (*currentendhinge[secnum]).nextendhinge;
		(*currentendhinge[secnum]).side = 1;
		(*currentendhinge[secnum]).num_of_walks = num_walks;
		(*currentendhinge[secnum]).edgecount = Num_edges;
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

/*		printf("Recorded 2-span. Section %d to %d, with hedges: ", inNum, outNum);
		for(i=0; i<=vM*vL-1; i++){
			printf("%d ", (*current_hinge_span[inNum]).hedges[i]);
		}
		printf("\n");
*/	

		//now record the mirror image left endhinge version.
		fillendtemplate2();
		fillrendtemplate2();

		secnum = num_section_endhinge_norder2();
//		printf("OG secnum2=%d\n", secnum);
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
//			sectionkey[i]=secnum;
//			secnum=i;
		}
//		printf("NEW secnum=%d\n", secnum);

		//Need to check if this left hinge will be a duplicate.

		//check if this endhinge configuration has already been recorded
		duplicate=0;
		test = firstendhinge[secnum];
		while(test->nextendhinge != NULL && duplicate==0){
	//		printf("looking through linked list\n");
			test = test->nextendhinge;
			if(test->side==0){	//if left endhinge
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
	//		printf("recording left endhinge\n");
			(*currentendhinge[secnum]).nextendhinge = newendhinge(num_walks); /* add a newendhinge to the linked list */
	
			currentendhinge[secnum] = (*currentendhinge[secnum]).nextendhinge;
			(*currentendhinge[secnum]).side = 0;
			(*currentendhinge[secnum]).num_of_walks = num_walks;
			(*currentendhinge[secnum]).edgecount = Num_edges;
			for(i=0; i<=vM*vL-1; i++){
				(*currentendhinge[secnum]).hedges[i]=temp_hedges[i];
			}

			for(i=0; i<=num_walks-1; i++){
				currentendhinge[secnum]->start[0][i] = 2;	//change
				currentendhinge[secnum]->end[0][i] = 2;		//change
				for(j=1; j<=2; j++){
					currentendhinge[secnum]->start[j][i] = curstart[j][i];
					currentendhinge[secnum]->end[j][i] = curend[j][i];
				}
				for(j=0; j<=vM*vL; j++){
					if(curwalks[j][i]==1){	//change
						currentendhinge[secnum]->walks[j][i]=2;
					}
					else if(curwalks[j][i]==2){	//change
						currentendhinge[secnum]->walks[j][i]=1;
					}
					else{
						currentendhinge[secnum]->walks[j][i] = curwalks[j][i];
					}
				}
			}
		}
		else{
			num_duplicate_left_endhinges++;
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
struct hinge_span* newhinge(unsigned int numberofwalks)		/* returns a pointer to a hinge_span structure */
/* this is a memory allocation function */
{
	struct hinge_span *nextnewhinge;

	nextnewhinge = ((struct hinge_span *) calloc(1, sizeof(struct hinge_span)));
	if (nextnewhinge == NULL) {
		fprintf(stderr, "unable to allocate memory");
		exit(1);
	}
	int i,j;
	nextnewhinge->hedges = (unsigned int*)calloc(vM * vL, sizeof(unsigned int));
	nextnewhinge->walks = (int**)calloc(vM * vL + 1, sizeof(int*));
	
	for(i=0; i<=vM*vL-1; i++){
		nextnewhinge->hedges[i]=0;
	}

	//allocate space for start, end, and walks
	for(i=0; i<=2; i++){
		nextnewhinge->start[i] = int_vecalloc(0,numberofwalks-1);
		nextnewhinge->end[i] = int_vecalloc(0,numberofwalks-1);
	}
	for(i=0; i<=vM*vL; i++){
		nextnewhinge->walks[i] = int_vecalloc(0,numberofwalks-1);
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

	return nextnewhinge;
}


/***************************************************************************/
struct endhinge* newendhinge(unsigned int numberofwalks)		/* returns a pointer to an endhinge structure */
/* this is a memory allocation function */
{
	struct endhinge *nextnewendhinge;

	nextnewendhinge = ((struct endhinge *) calloc(1, sizeof(struct endhinge)));
	if (nextnewendhinge == NULL) {
		fprintf(stderr, "unable to allocate memory");
		exit(1);
	}
	int i,j;
	nextnewendhinge->hedges = (unsigned int*)calloc(vM * vL, sizeof(unsigned int));
	nextnewendhinge->walks = (int**)calloc(vM * vL + 1, sizeof(int*));

	for(i=0; i<=vM*vL-1; i++){
		nextnewendhinge->hedges[i]=0;
	}

	//allocate space for start, end, and walks
	for(i=0; i<=2; i++){
		nextnewendhinge->start[i] = int_vecalloc(0,numberofwalks-1);
		nextnewendhinge->end[i] = int_vecalloc(0,numberofwalks-1);
	}
	for(i=0; i<=vM*vL; i++){
		nextnewendhinge->walks[i] = int_vecalloc(0,numberofwalks-1);
	}

	//initialize start, end, and walks to all zeros
	for(i=0; i<=numberofwalks-1; i++){
		for(j=0; j<=2; j++){
			nextnewendhinge->start[j][i]=0;
			nextnewendhinge->end[j][i]=0;
		}
		for(j=0; j<=vM*vL; j++){
			nextnewendhinge->walks[j][i]=0;
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

//			printf("a=%d,b=%d\n", a, b);

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

void check_reachability() {
    if (!ham_check) return;
    
    printf("\n--- HAMILTONIAN REACHABILITY CHECK ---\n");
    unsigned long int testsection = 1; 
    unsigned char *testvector = (unsigned char*)calloc(max_sections + 1, sizeof(unsigned char));
    unsigned char *temptestvector = (unsigned char*)calloc(max_sections + 1, sizeof(unsigned char));
    unsigned char *testvector2 = NULL;
    unsigned char *temptestvector2 = NULL;
    
    if (M == 2 && L == 2) {
        testvector2 = (unsigned char*)calloc(max_sections + 1, sizeof(unsigned char));
        temptestvector2 = (unsigned char*)calloc(max_sections + 1, sizeof(unsigned char));
    }

    testvector[testsection] = 1;
    unsigned long int testpower = 200; // Sufficient for most reachable spaces
    
    printf("Performing reachability analysis (power=%lu)...\n", testpower);
    for (unsigned long int k = 1; k <= testpower; k++) {
        for (int i = 1; i <= max_sections; i++) temptestvector[i] = 0;
        for (int i = 1; i <= max_sections; i++) {
            if (testvector[i]) {
                for (int j = 1; j <= num_outsections[i]; j++) {
                    unsigned long int out = t_outsection[i][j];
                    if (out > 0) temptestvector[out] = 1;
                }
            }
        }
        int changed = 0;
        for (int i = 1; i <= max_sections; i++) {
            if (testvector[i] != temptestvector[i]) changed = 1;
            testvector[i] = temptestvector[i];
        }
        if (!changed && k > 10) break; // Converged
    }
    
    if (M == 2 && L == 2) {
        testvector2[testsection] = 1;
        for (unsigned long int k = 1; k <= testpower + 1; k++) {
            for (int i = 1; i <= max_sections; i++) temptestvector2[i] = 0;
            for (int i = 1; i <= max_sections; i++) {
                if (testvector2[i]) {
                    for (int j = 1; j <= num_outsections[i]; j++) {
                        unsigned long int out = t_outsection[i][j];
                        if (out > 0) temptestvector2[out] = 1;
                    }
                }
            }
            int changed = 0;
            for (int i = 1; i <= max_sections; i++) {
                if (testvector2[i] != temptestvector2[i]) changed = 1;
                testvector2[i] = temptestvector2[i];
            }
            if (!changed && k > 10) break;
        }
    }

    unsigned long int numconnect = 0;
    unsigned long int numnoconnect = 0;
    for (int i = 1; i <= max_sections; i++) {
        if (num_outsections[i] > 0) {
            int reachable = 0;
            if (M == 2 && L == 2) reachable = (testvector[i] || testvector2[i]);
            else reachable = testvector[i];
            
            if (!reachable) {
                numnoconnect++;
                // Invalidate transitions
                for (int j = 1; j <= num_outsections[i]; j++) {
                    t_outsection[i][j] = 0;
                    t_num_walks[i][j] = 0;
                }
                num_outsections[i] = 0;
                num_left_endhinges[i] = 0;
                num_right_endhinges[i] = 0;
            } else {
                numconnect++;
            }
        }
    }
    
    // Re-index transitions (t_nrr)
    unsigned long int new_tspan_num = 0;
    for (int i = 1; i <= max_sections; i++) {
        for (int j = 1; j <= num_outsections[i]; j++) {
            t_nrr[i][j] = ++new_tspan_num;
        }
    }
    
    printf("Reachable Hamiltonian sections: %lu (Removed %lu isolated states)\n", numconnect, numnoconnect);
    printf("Total valid Hamiltonian 2-spans: %lu\n", new_tspan_num);
    num_tspans = new_tspan_num;

    free(testvector); free(temptestvector);
    if (testvector2) { free(testvector2); free(temptestvector2); }
}

void printtofile(){
	if(filetotal>=maxpolys){
		//need to start new file
		printf("Finished filling file %lu. It contains %lu polys. Creating a new one.\n", filenum, filetotal);
		fprintf(fp, "-999\n");
		fclose(fp);
		filenum++;
		const char *file_prefix = (ham_check) ? "MCpolysHam" : "MCpolys";
		sprintf(filename, "%s/%sL%dM%dspan%drun%dnum%lu.txt", output_dir, file_prefix, L, M, totalspan, runnum, filenum);
		fp = fopen(filename, "w");	//create or overwrite the file "filename"

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
	filetotal++;
}



#if defined(CS)
#include "archive_deps/topology/LFlag_0.c"		/* include if M = 0  */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the left and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

//#include "../topology/RFlag_0.c"		/* include if M = 0  */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the right and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#include "archive_deps/utils/noncrossing.c"	/* include if M = 0 */
/* This function takes arguements (int side, int a, int b, int c, int d) */
/* and returns 1 if the walk connecting a and b does not cross */
/* the walk connecting c and d it returns 0 otherwise*/
#include "archive_deps/sections/cstatenum.c"
#else
#include "archive_deps/topology/LFlag_norder.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the left and returns 0 otherwise */
/* it also uses the global variable ordertemplate */
#include "archive_deps/topology/LFlag_endhinge_norder.c"

#define tspans_outsection t_outsection
#define tspans_nrr t_nrr
#define R_Evector R_Evector_solve
#include "archive_deps/transfer_matrix/pw_meth_ts_LRvec_fcheck.c"
#undef tspans_outsection
#undef tspans_nrr
#undef R_Evector

//#include "../topology/RFlag.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the right and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#if vM*vL<7
#include "archive_deps/sections/Num_section_6V.c"  /*unsigned long int num_section(int side)*/
/* This function takes the arguement (int side) and assigns a unique number */
/* to the section on side = side */
#include "archive_deps/sections/Num_section_6V_nonordered.c"
#include "archive_deps/sections/Num_section_6V_endhinge_nonordered.c"
#elif vM*vL<9
#include "archive_deps/sections/Num_section_8V.c"
#include "archive_deps/sections/Num_section_8V_nonordered.c"
#include "archive_deps/sections/Num_section_8V_endhinge_nonordered.c"
#elif vM*vL<11
#include "archive_deps/sections/Num_section_10V.c"
#include "archive_deps/sections/Num_section_10V_nonordered.c"
#include "archive_deps/sections/Num_section_10V_endhinge_nonordered.c"
#elif vM*vL<13
#include "archive_deps/sections/Num_section_12V.c"
#include "archive_deps/sections/Num_section_12V_nonordered.c"
#include "archive_deps/sections/Num_section_12V_endhinge_nonordered.c"
#endif

#endif



#include "archive_deps/sections/printsection.c"

#include "archive_deps/utils/int_vecalloc.c"
#include "archive_deps/utils/unsgn_vecalloc.c"

#include "archive_deps/utils/unsgnlong_vecalloc.c" /*unsigned long int unsgnlong_vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with unsigned long integer enteries */

#include "archive_deps/utils/vecalloc.c"		/*vec_ent * vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with enteries of type vec_ent */

#include "archive_deps/utils/matalloc.c"		/*mat_ent **matalloc(int rowlow, int rowhigh, int collow, int colhigh)*/
/* used for dynamically allocating memory for matrices with enteries of type mat_ent */





	









