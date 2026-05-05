#ifndef MC_GLOBALS_H
#define MC_GLOBALS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>

// Macros
#define vec_ent double
#define mat_ent double
#define newline printf("\n")
#define VALID 1
#define NOTVALID !VALID

// Global variable declarations
extern int L;
extern int M;
extern int totalspan;
extern int samplesize;
extern int runnum;
extern unsigned int seednum;
extern int max_sections;
extern unsigned long int max_keynum;
extern int max_tspans;
extern double dom_evalue;
extern int maxpolys;
extern int vec_length;
extern int vM;
extern int vL;
extern int CS_mode;
extern int mode;
extern int ham_check;
extern char* output_dir;
extern double fval;

typedef struct PolygonState {
	unsigned int ***ordertemplate;
	unsigned int ***reordertemplate;
	unsigned short int **hingestatus;
	unsigned short int **alreadyentered;
	unsigned short int **colhingeedges;
	unsigned short int **rowhingeedges;
	int **curstart;
	int **curend;
	int **curwalks;
	unsigned int *num_walks;
} PolygonState;

extern PolygonState primary_polygon_state;

// Structs
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

}	; /* these are global pointers now */

extern struct hinge_span **first_hinge_span, **current_hinge_span;

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

extern struct endhinge **firstendhinge;
extern struct endhinge **currentendhinge;

// Arrays
extern unsigned int ***ordertemplate;
extern unsigned long int num_2_spans;
extern unsigned long int valid_2_spans;

extern unsigned long int *num_outsections;
extern unsigned long int **t_outsection;
extern unsigned long int **t_nrr;
extern unsigned int **t_num_walks;
extern int ****t_start;
extern int ****t_end;
extern int ****t_walks;
extern unsigned int **tspans_edges;
extern double *L_Evector[2];
extern double *R_Evector_solve[2];

extern unsigned long int *num_left_endhinges;
extern unsigned int **Lend_num_walks;
extern int ****Lend_start;
extern int ****Lend_end;
extern int ****Lend_walks;

extern unsigned long int *num_right_endhinges;
extern unsigned int **Rend_num_walks;
extern int ****Rend_start;
extern int ****Rend_end;
extern int ****Rend_walks;

extern unsigned int ***reordertemplate;
extern unsigned short int **hingestatus;
extern unsigned short int **alreadyentered;
extern unsigned short int **colhingeedges;
extern unsigned short int **rowhingeedges;
extern unsigned long int *sectionkey;
extern unsigned long int dupcounter;

extern int **curstart;
extern int **curend;
extern int **curwalks;
extern unsigned int num_walks;

extern int EndOrdNum[2];
extern int **endtemplate;
extern int **rendtemplate;
extern int **endtemplate2;
extern int **rendtemplate2;

extern unsigned long int tot_right_endhinges;
extern unsigned long int tot_left_endhinges;
extern unsigned long int num_duplicate_right_endhinges;
extern unsigned long int num_duplicate_left_endhinges;
extern unsigned long int num_tspans;
extern unsigned long int num_duplicate_tspans;

extern int **built_walks_start;
extern int **built_walks_end;
extern int **built_walks_direcs;
extern int num_built_walks;

extern char filename[100];
extern FILE* fp;
extern unsigned long int filetotal;
extern unsigned long int filenum;

// Function prototypes
void set_system_params();
void allocate_globals();
int run_integrated_2sap_sampler(int argc, char *argv[]);
int run_integrated_2sap_ham_sampler(int argc, char *argv[]);
void generate_evectors();
void conv_to_array(void);
void conv_endhinges_to_array(void);
void enterhinge(PolygonState *poly, int i, int j, int side, int (*pointordNum)[3], int curlength);
void leavehinge(PolygonState *poly, int i, int j, int side, int (*pointordNum)[3], int curlength);
void rowedges(PolygonState *poly, int i, int j, int (*pointordNum)[3], int curlength);
void coledges(PolygonState *poly, int i, int j, int (*pointordNum)[3], int curlength);
void recordtemplate( int (*pointordNum)[3]);
void enterendhinge(int i, int j, int curlength);
void leaveendhinge(int i, int j, int curlength);
void endrowedges(int i, int j, int curlength);
void endcoledges(int i, int j, int curlength);
void recordendtemplate();
struct hinge_span *newhinge(unsigned int numberofwalks);
struct endhinge *newendhinge(unsigned int numberofwalks);
void fillreordertemplate(int ledges, int redges);
void printordtemp();
void printreordtemp();
void fillrendtemplate();
void fillendtemplate2();
void fillrendtemplate2();
void printbuiltwalks();
void printbuiltwalks_all();
void add_to_built_walks(unsigned long int secnum, int nth_tspan);
void add_right_endhinge(unsigned long int secnum, int nth_endhinge);
int reverse_direc(int direc);
void check_reachability();
void printtofile();

unsigned long int num_section(int side);
unsigned long int num_section_norder(int side);
unsigned long int num_section_endhinge_norder(void);
unsigned long int num_section_endhinge_norder2(void);
unsigned short int LFlag( int (*)[3]);
unsigned short int LFlag_endhinge(int (*EndOrdNum)[2]);
int *int_vecalloc(int low, int high);
unsigned int *unsgn_vecalloc(int low, int high);
unsigned long int *unsgnlong_vecalloc(int low, int high);

#endif
