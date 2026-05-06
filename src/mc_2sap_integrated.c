//This program will sample SAPs of a certain span uniformly at random (by span).
//uses 2SAP_TMcalc_PrintEvectors.c, which prints out eigenvectors to:
//2SAP_L_Evector_TS_L%dM%d.txt
//2SAP_R_Evector_TS_L%dM%d.txt
//These files must be in same folder.

//compliled with: gcc -lm -Wall -O3 -o 2SAP_MCsample.o 2SAP_MCsample.c
//ran using ./2SAP_MCsample.o

//sample 2SAPs will go to files: MC2SAPsL%dM%dspan%dnum%d.txt

#include <stdio.h>

#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include "mc_runtime.h"
#include "mc_2sap_common.h"
#include "mc_spectral.h"
#include "run_metadata.h"

extern int L;
extern int M;
extern int totalspan;
extern int samplesize;
extern int runnum;
extern unsigned int seednum;
extern int maxpolys;
extern int mode;
extern int ham_check;

extern int max_sections;
extern int max_tspans;
extern unsigned long int max_keynum;
extern double dom_evalue;
extern double fval;

static double mc2_max_eval_LRvec(double fugacity);

double *mc2_MC_L_Evector[2];
double *mc2_MC_R_Evector[2];
unsigned long int **mc2_MC_tspans_edges;


#define MAX_vM 5
#define MAX_vL 5
#define MAX_vMvL 25
#define MAX_SPAN 20
#define MAX_KEYNUM_ARR 140000

#define vM (M+1)
#define vL (L+1)
		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
/* Source-level namespace for compatibility 2SAP helper includes. */
#define EndOrdNum mc2_EndOrdNum
#define EndOrdNum2 mc2_EndOrdNum2
#define LFlag2 mc2_LFlag2
#define LFlag_endhinge2 mc2_LFlag_endhinge2
#define Lend_end mc2_Lend_end
#define Lend_end2 mc2_Lend_end2
#define Lend_num_walks mc2_Lend_num_walks
#define Lend_num_walks2 mc2_Lend_num_walks2
#define Lend_start mc2_Lend_start
#define Lend_start2 mc2_Lend_start2
#define Lend_walks mc2_Lend_walks
#define Lend_walks2 mc2_Lend_walks2
#define MC_L_Evector mc2_MC_L_Evector
#define MC_R_Evector mc2_MC_R_Evector
#define MC_tspans_edges mc2_MC_tspans_edges
#define Rend_end mc2_Rend_end
#define Rend_end2 mc2_Rend_end2
#define Rend_num_walks mc2_Rend_num_walks
#define Rend_num_walks2 mc2_Rend_num_walks2
#define Rend_start mc2_Rend_start
#define Rend_start2 mc2_Rend_start2
#define Rend_walks mc2_Rend_walks
#define Rend_walks2 mc2_Rend_walks2
#define add_right_endhinge mc2_add_right_endhinge
#define add_to_built_walks mc2_add_to_built_walks
#define alreadyentered mc2_alreadyentered
#define alreadyentered2 mc2_alreadyentered2
#define built_walks_direcs mc2_built_walks_direcs
#define built_walks_direcs2 mc2_built_walks_direcs2
#define built_walks_end mc2_built_walks_end
#define built_walks_end2 mc2_built_walks_end2
#define built_walks_start mc2_built_walks_start
#define built_walks_start2 mc2_built_walks_start2
#define coledges mc2_coledges
#define coledges2 mc2_coledges2
#define colhingeedges mc2_colhingeedges
#define conv_endhinges_to_array mc2_conv_endhinges_to_array
#define conv_to_array mc2_conv_to_array
#define curend mc2_curend
#define curend2 mc2_curend2
#define current_hinge_span mc2_current_hinge_span
#define currentendhinge mc2_currentendhinge
#define curstart mc2_curstart
#define curstart2 mc2_curstart2
#define curwalks mc2_curwalks
#define curwalks2 mc2_curwalks2
#define dupcounter mc2_dupcounter
#define endcoledges mc2_endcoledges
#define endcoledges2 mc2_endcoledges2
#define endrowedges mc2_endrowedges
#define endrowedges2 mc2_endrowedges2
#define endtemplate mc2_endtemplate
#define endtemplate2 mc2_endtemplate2
#define endtemplate3 mc2_endtemplate3
#define endtemplate4 mc2_endtemplate4
#define enterendhinge mc2_enterendhinge
#define enterendhinge2 mc2_enterendhinge2
#define enterhinge mc2_enterhinge
#define enterhinge2 mc2_enterhinge2
#define filename mc2_filename
#define filenum mc2_filenum
#define filetotal mc2_filetotal
#define fillendtemplate2 mc2_fillendtemplate2
#define fillendtemplate4 mc2_fillendtemplate4
#define fillrendtemplate mc2_fillrendtemplate
#define fillrendtemplate2 mc2_fillrendtemplate2
#define fillrendtemplate3 mc2_fillrendtemplate3
#define fillrendtemplate4 mc2_fillrendtemplate4
#define fillreordertemplate mc2_fillreordertemplate
#define fillreordertemplate2 mc2_fillreordertemplate2
#define first_hinge_span mc2_first_hinge_span
#define firstendhinge mc2_firstendhinge
#define fp mc2_fp
#define hingestatus mc2_hingestatus
#define initran_ mc2_initran_
#define int_vecalloc mc2_int_vecalloc
#define leaveendhinge mc2_leaveendhinge
#define leaveendhinge2 mc2_leaveendhinge2
#define leavehinge mc2_leavehinge
#define leavehinge2 mc2_leavehinge2
#define matalloc mc2_matalloc
#define max_eval_LRvec mc2_max_eval_LRvec
#define newendhinge mc2_newendhinge
#define newhinge mc2_newhinge
#define noncrossing mc2_noncrossing
#define num_2_spans mc2_num_2_spans
#define num_built_walks mc2_num_built_walks
#define num_built_walks2 mc2_num_built_walks2
#define num_duplicate_left_endhinges mc2_num_duplicate_left_endhinges
#define num_duplicate_right_endhinges mc2_num_duplicate_right_endhinges
#define num_duplicate_tspans mc2_num_duplicate_tspans
#define num_left_endhinges mc2_num_left_endhinges
#define num_nonlexi_endhinges mc2_num_nonlexi_endhinges
#define num_outsections mc2_num_outsections
#define num_potential_endhinges mc2_num_potential_endhinges
#define num_right_endhinges mc2_num_right_endhinges
#define num_section mc2_num_section
#define num_section_0V mc2_num_section_0V
#define num_section_10V mc2_num_section_10V
#define num_section_12V mc2_num_section_12V
#define num_section_6V mc2_num_section_6V
#define num_section_8V mc2_num_section_8V
#define num_section_endhinge_norder mc2_num_section_endhinge_norder
#define num_section_endhinge_norder2 mc2_num_section_endhinge_norder2
#define num_section_endhinge_norder2_10V mc2_num_section_endhinge_norder2_10V
#define num_section_endhinge_norder2_12V mc2_num_section_endhinge_norder2_12V
#define num_section_endhinge_norder2_6V mc2_num_section_endhinge_norder2_6V
#define num_section_endhinge_norder2_8V mc2_num_section_endhinge_norder2_8V
#define num_section_endhinge_norder3 mc2_num_section_endhinge_norder3
#define num_section_endhinge_norder3_10V mc2_num_section_endhinge_norder3_10V
#define num_section_endhinge_norder3_12V mc2_num_section_endhinge_norder3_12V
#define num_section_endhinge_norder3_6V mc2_num_section_endhinge_norder3_6V
#define num_section_endhinge_norder3_8V mc2_num_section_endhinge_norder3_8V
#define num_section_endhinge_norder4 mc2_num_section_endhinge_norder4
#define num_section_endhinge_norder4_10V mc2_num_section_endhinge_norder4_10V
#define num_section_endhinge_norder4_12V mc2_num_section_endhinge_norder4_12V
#define num_section_endhinge_norder4_6V mc2_num_section_endhinge_norder4_6V
#define num_section_endhinge_norder4_8V mc2_num_section_endhinge_norder4_8V
#define num_section_endhinge_norder_10V mc2_num_section_endhinge_norder_10V
#define num_section_endhinge_norder_12V mc2_num_section_endhinge_norder_12V
#define num_section_endhinge_norder_6V mc2_num_section_endhinge_norder_6V
#define num_section_endhinge_norder_8V mc2_num_section_endhinge_norder_8V
#define num_section_norder mc2_num_section_norder
#define num_section_norder2 mc2_num_section_norder2
#define num_section_norder2_10V mc2_num_section_norder2_10V
#define num_section_norder2_12V mc2_num_section_norder2_12V
#define num_section_norder2_6V mc2_num_section_norder2_6V
#define num_section_norder2_8V mc2_num_section_norder2_8V
#define num_section_norder_10V mc2_num_section_norder_10V
#define num_section_norder_12V mc2_num_section_norder_12V
#define num_section_norder_6V mc2_num_section_norder_6V
#define num_section_norder_8V mc2_num_section_norder_8V
#define num_tspans mc2_num_tspans
#define num_valid_endhinges mc2_num_valid_endhinges
#define num_walks mc2_num_walks
#define num_walks2 mc2_num_walks2
#define ordertemplate mc2_ordertemplate
#define ordertemplate2 mc2_ordertemplate2
#define paircomp mc2_paircomp
#define printbuiltwalks mc2_printbuiltwalks
#define printbuiltwalks2_all mc2_printbuiltwalks2_all
#define printbuiltwalks_all mc2_printbuiltwalks_all
#define printordtemp mc2_printordtemp
#define printreordtemp mc2_printreordtemp
#define printsection mc2_printsection
#define printtofile mc2_printtofile
#define ran1real_ mc2_ran1real_
#define recordendtemplate mc2_recordendtemplate
#define recordtemplate mc2_recordtemplate
#define rendtemplate mc2_rendtemplate
#define rendtemplate2 mc2_rendtemplate2
#define rendtemplate3 mc2_rendtemplate3
#define rendtemplate4 mc2_rendtemplate4
#define reordertemplate mc2_reordertemplate
#define reordertemplate2 mc2_reordertemplate2
#define reverse_direc mc2_reverse_direc
#define rowedges mc2_rowedges
#define rowedges2 mc2_rowedges2
#define rowhingeedges mc2_rowhingeedges
#define sectionkey mc2_sectionkey
#define sectionkey2SAP mc2_sectionkey2SAP
#define t_end mc2_t_end
#define t_end2 mc2_t_end2
#define t_nrr mc2_t_nrr
#define t_num_walks mc2_t_num_walks
#define t_num_walks2 mc2_t_num_walks2
#define t_outsection mc2_t_outsection
#define t_start mc2_t_start
#define t_start2 mc2_t_start2
#define t_walks mc2_t_walks
#define t_walks2 mc2_t_walks2
#define tot_left_endhinges mc2_tot_left_endhinges
#define tot_right_endhinges mc2_tot_right_endhinges
#define unsgn_vecalloc mc2_unsgn_vecalloc
#define unsgnlong_vecalloc mc2_unsgnlong_vecalloc
#define valid_2_spans mc2_valid_2_spans
#define vecalloc mc2_vecalloc
#include "../include/marsaglia.h"


/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/
//for now, need L,M>0, totalspan>0,

				// number of horizontal edges
				// number of vertical edges
		// span of sampled SAPs
	// number of samples desired



	// number of 2SAPs per file.

/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/


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

	unsigned int		hedges[MAX_vMvL];	//contains info about which hinge-edges are filled

	int		*start[3];	//start[0][i] contains the starting x-coordinate of the i'th walk.
						//start[1][i] contains the starting y-coordinate of the i'th walk.
						//start[2][i] contains the starting z-coordinate of the i'th walk.
	int		*end[3];		//end[0][i] contains the ending x-coordinate of the i'th walk.
						//end[1][i] contains the ending y-coordinate of the i'th walk.
						//end[2][i] contains the ending z-coordinate of the i'th walk.
	int		*walks[MAX_vMvL+1];	//walks[j][i] contains the j'th step in the i'th walk.	rest zero's.
	unsigned int	num_of_walks;	//number of walks in this 2-span.

	int		*start2[3];	//same but for 2nd poly
	int		*end2[3];		//same but for 2nd poly
	int		*walks2[MAX_vMvL+1];	//same but for 2nd poly
	unsigned int	num_of_walks2;	//same but for 2nd poly

	struct hinge_span	*nexthinge;	// pointer to another hinge_span data structure

} **mc2_first_hinge_span, **mc2_current_hinge_span; /* these are global */

	/* mc2_first_hinge_span is an array of pointers to hinge_span data structures where each pointer	*/
	/* in the array points to the first element of a linked list of hinge_span data structures	*/
	/* The array index implicely defines the "inorder" (or the lefthand section of the two_span)  	*/

	/* mc2_current_hinge_span is a dummy pointer used for traversing the linked lists mentioned above	*/
	/* an array of dummy pointers is needed because there is more than one linked list to consider	*/


struct endhinge {
	unsigned long int	sec2;			//section number of poly2
	unsigned int		hedges[MAX_vMvL];	//contains info about which hinge-edges are filled
	unsigned short int	side;		//0 if a left endhinge, 1 if a right endhinge

	int		*start[3];	//start[0][i] contains the starting x-coordinate of the i'th walk.
						//start[1][i] contains the starting y-coordinate of the i'th walk.
						//start[2][i] contains the starting z-coordinate of the i'th walk.
	int		*end[3];		//end[0][i] contains the ending x-coordinate of the i'th walk.
						//end[1][i] contains the ending y-coordinate of the i'th walk.
						//end[2][i] contains the ending z-coordinate of the i'th walk.
	int		*walks[MAX_vMvL+1];	//walks[j][i] contains the j'th step in the i'th walk.	rest zero's.
	unsigned int	num_of_walks;	//number of walks in this endhinge

	int		*start2[3];		//same but for 2nd poly
	int		*end2[3];			//same but for 2nd poly
	int		*walks2[MAX_vMvL+1];	//same but for 2nd poly
	unsigned int	num_of_walks2;	//same but for 2nd poly

	struct endhinge* nextendhinge;// pointer to the next endhinge structure in the linked list
};

struct endhinge **mc2_firstendhinge;		/* Array of pointers to endhinge structures, where each	*/
										/* pointer in the array points to the first element in		*/
										/* the linked list of endhinge structures				*/
struct endhinge **mc2_currentendhinge;	/* Array of dummy pointers that are used to traverse		*/
										/* the endhinge linked lists.							*/



/***************************************************************************/
/****************** end of structures **************************************/
/***************************************************************************/


/***************************************************************************/
/*************** Global variables used in this program *********************/
/***************************************************************************/

unsigned int		mc2_ordertemplate[2][MAX_vM][MAX_vL];	/* The first index is for side (left = 0 ; right = 1) */
										/* the second index is for the location of the vertex from the top */
										/* the third index is for the location of the vertex from the left */
										/* This global array is used as a "template" to store information  */
										/* on the edges in the two sections that make up a two span        */
										/* The array is initiallized to zero and numbers are stored at     */
										/* the appropriate locations (corresponding to given vertices) as  */
										/* the potential two spans are built up algorithmically by the     */
										/* program. The stored numbers represent the order in which the    */
										/* edge is traversed when the two span is part of a polygonal walk */
unsigned int		mc2_ordertemplate2[2][MAX_vM][MAX_vL]; //same but for 2nd poly

unsigned long int 	mc2_num_2_spans = 1;	/* This variable keeps track of the total number of potential two spans   */

unsigned long int	mc2_valid_2_spans = 0;	/* This is used to count the number of two spans which are actually valid */

////////////////////// arrays which hold the tspan info ////////////////////////////////

unsigned long int	mc2_num_outsections[MAX_KEYNUM_ARR+1];	/* For a given section (which has a number section_num assigned to it)*/
										/* mc2_num_outsections[section_num] is the number of two spans with this  */
										/* section as a first section                                         */
 
unsigned long int	*mc2_t_outsection[MAX_KEYNUM_ARR+1];	/* An array of pointers to integers 					*/
										/* records the "outsection" of the two-span				*/
										/* mc2_t_outsection[insection][nth_tspan]				*/
										/* nth_tspan = nth two-span with firstsection = insection	*/

unsigned long int	*mc2_t_nrr[MAX_KEYNUM_ARR+1];		/*An array of pointers to integers*/
										/*records the tspan "number"      */
										/* defined by the array index */
										/*example:		*/
										/*mc2_t_nrr[1][1] gives the two-span number */
										/* of the two-span with first section = 1 */
										/* and the second section is the first such */
										/* section which connects with first section = 1 */

unsigned int		*mc2_t_num_walks[MAX_KEYNUM_ARR+1];	//mc2_t_num_walks[insection][nth_tspan] holds the number of walks in this tspan

int				***mc2_t_start[MAX_KEYNUM_ARR+1];	//mc2_t_start[insection][nth_tspan][i][nth_walk] holds the i-coordinate (0=x, 1=y, 2=z) of the start point of the nth_walk

int				***mc2_t_end[MAX_KEYNUM_ARR+1];		//mc2_t_end[insection][nth_tspan][i][nth_walk] holds the i-coordinate (0=x, 1=y, 2=z) of the end point of the nth_walk

int				***mc2_t_walks[MAX_KEYNUM_ARR+1];	//mc2_t_end[insection][nth_tspan][direc][nth_walk] holds the "direc'th" direction the nth_walk

unsigned int		*mc2_t_num_walks2[MAX_KEYNUM_ARR+1];	//same but for poly2
int				***mc2_t_start2[MAX_KEYNUM_ARR+1];	//same but for poly2
int				***mc2_t_end2[MAX_KEYNUM_ARR+1];		//same but for poly2
int				***mc2_t_walks2[MAX_KEYNUM_ARR+1];	//same but for poly2



////////////////////// arrays which hold the endhinge info ////////////////////////////////////

unsigned long int	mc2_num_left_endhinges[MAX_KEYNUM_ARR+1];
unsigned int		*mc2_Lend_num_walks[MAX_KEYNUM_ARR+1];
int				***mc2_Lend_start[MAX_KEYNUM_ARR+1];
int				***mc2_Lend_end[MAX_KEYNUM_ARR+1];
int				***mc2_Lend_walks[MAX_KEYNUM_ARR+1];

unsigned int		*mc2_Lend_num_walks2[MAX_KEYNUM_ARR+1];
int				***mc2_Lend_start2[MAX_KEYNUM_ARR+1];
int				***mc2_Lend_end2[MAX_KEYNUM_ARR+1];
int				***mc2_Lend_walks2[MAX_KEYNUM_ARR+1];



unsigned long int	mc2_num_right_endhinges[MAX_KEYNUM_ARR+1];
unsigned int		*mc2_Rend_num_walks[MAX_KEYNUM_ARR+1];
int				***mc2_Rend_start[MAX_KEYNUM_ARR+1];
int				***mc2_Rend_end[MAX_KEYNUM_ARR+1];
int				***mc2_Rend_walks[MAX_KEYNUM_ARR+1];

unsigned int		*mc2_Rend_num_walks2[MAX_KEYNUM_ARR+1];
int				***mc2_Rend_start2[MAX_KEYNUM_ARR+1];
int				***mc2_Rend_end2[MAX_KEYNUM_ARR+1];
int				***mc2_Rend_walks2[MAX_KEYNUM_ARR+1];



////////////////

unsigned int mc2_reordertemplate[2][MAX_vM][MAX_vL];	//holds the re-ordered mc2_ordertemplate to reduce number of unique sections.
unsigned int mc2_reordertemplate2[2][MAX_vM][MAX_vL];	//same but for poly2

unsigned short int mc2_hingestatus[MAX_vM][MAX_vL];		//keeps track of which vertices in the hinge are occupied. 0=empty, 1=filled.
unsigned short int mc2_alreadyentered[MAX_vM][MAX_vL];	//keeps track of where we've already initially entered. 0=no, 1=yes.
unsigned short int mc2_alreadyentered2[MAX_vM][MAX_vL];	//same but for poly2

unsigned short int mc2_colhingeedges[MAX_vM][MAX_vL]; //keeps track of which column (i/M) edges in the hinge are occupied. 0=empty, 1=filled.
unsigned short int mc2_rowhingeedges[MAX_vM][MAX_vL]; //keeps track of which row (j/L) edges in the hinge are occupied. 0=empty, 1=filled.

unsigned long int mc2_sectionkey[MAX_KEYNUM_ARR+1];	//will hold valid (re-ordered) section nums. index is the key num
									//note: mc2_sectionkey[0] is garbage. Key numbers start at 1.

unsigned long int mc2_sectionkey2SAP[MAX_KEYNUM_ARR+1][2];

unsigned long int mc2_dupcounter=0;	//will count the number of duplicate 2-spans created and not stored.

int mc2_curstart[3][MAX_vMvL];		//will hold the start-points of the walks being built
int mc2_curend[3][MAX_vMvL];		//will hold the end-points of the walks being built
int mc2_curwalks[MAX_vMvL+1][MAX_vMvL];	//will hold the walks as they are being built.
						//mc2_curwalks[j][i] contains the j'th step in the i'th walk.	rest zero's.

int mc2_curstart2[3][MAX_vMvL];			//same but for poly2
int mc2_curend2[3][MAX_vMvL];			//same but for poly2
int mc2_curwalks2[MAX_vMvL+1][MAX_vMvL];	//same but for poly2

unsigned int mc2_num_walks=0;	//number of walks in the 2-span
unsigned int mc2_num_walks2=0;	//same but for poly2

//endhinge stuff
int mc2_EndOrdNum[2];		//contains info on # of edges in left section and endhinge
int mc2_EndOrdNum2[2];		//contains info on # of edges in left section and endhinge

int mc2_endtemplate[MAX_vM][MAX_vL];		//contains the section associated with the end hinge
int mc2_rendtemplate[MAX_vM][MAX_vL];	//holds the re-ordered mc2_endtemplate to reduce number of unique sections.
int mc2_endtemplate2[MAX_vM][MAX_vL];
int mc2_rendtemplate2[MAX_vM][MAX_vL];

int mc2_endtemplate3[MAX_vM][MAX_vL];	//same but for 2nd poly
int mc2_rendtemplate3[MAX_vM][MAX_vL];
int mc2_endtemplate4[MAX_vM][MAX_vL];
int mc2_rendtemplate4[MAX_vM][MAX_vL];

unsigned long int mc2_tot_right_endhinges=0;
unsigned long int mc2_tot_left_endhinges=0;
unsigned long int mc2_num_duplicate_right_endhinges=0;
unsigned long int mc2_num_duplicate_left_endhinges=0;
unsigned long int mc2_num_tspans=0;
unsigned long int mc2_num_duplicate_tspans=0;

unsigned long int mc2_num_potential_endhinges=0;
unsigned long int mc2_num_valid_endhinges=0;

unsigned long int mc2_num_nonlexi_endhinges=0;


//integer division rounds down in case of vM*vL is odd
int mc2_built_walks_start[MAX_vMvL/2][3];
int mc2_built_walks_end[MAX_vMvL/2][3];
int mc2_built_walks_direcs[MAX_vMvL/2][MAX_vMvL*(MAX_SPAN+1)];

int mc2_num_built_walks;

int mc2_built_walks_start2[MAX_vMvL/2][3];	//2nd poly
int mc2_built_walks_end2[MAX_vMvL/2][3];	//2nd poly
int mc2_built_walks_direcs2[MAX_vMvL/2][MAX_vMvL*(MAX_SPAN+1)];	//2nd poly

int mc2_num_built_walks2;	//2nd poly

char 			mc2_filename[100];	//this will hold the mc2_filename of the text file that we want to write to
FILE*			mc2_fp;			//file pointer. will point to the file that is to be written to.
unsigned long int mc2_filetotal=0;	//number of polys in the file 
unsigned long int mc2_filenum=1;		//number of polygon files created.





/* note: The following are also global: mc2_first_hinge_span, mc2_current_hinge_span */
/* They are discussed in the structures section above */

/***************************************************************************/
/******************* end of global variables  ******************************/
/***************************************************************************/

/***************************************************************************/
/******************* functions used in this program  ***********************/
/***************************************************************************/

int			system(const char *string);	/* called by main */
			/* used solely for printing the date and timing the program */

/* mc2_conv_to_array	is only called by main and is declared in main */
/* mc2_conv_endhinges_to_array	is only called by main and is declared in main */

/****  two_span builder functions  **********************************************************/

void            	mc2_enterhinge(int, int, int, int (*)[6], int);
void            	mc2_leavehinge(int, int, int, int (*)[6], int);
void            	mc2_rowedges(int, int, int (*)[6], int);
void            	mc2_coledges(int, int, int (*)[6], int);

void            	mc2_enterhinge2(int, int, int, int (*)[6], int);
void            	mc2_leavehinge2(int, int, int, int (*)[6], int);
void            	mc2_rowedges2(int, int, int (*)[6], int);
void            	mc2_coledges2(int, int, int (*)[6], int);

/****  two_span validity verifier functions  ******************* *****************************/

unsigned short int	mc2_LFlag2( int (*)[6]);
//unsigned short int	RFlag( int (*)[3]);
unsigned short int	mc2_noncrossing(int, int, int, int, int);	/* used in 0_L cases only   */

/****  one_span/section recorder functions  *************************************************/

void            	mc2_recordtemplate( int (*)[6]);
unsigned long int	mc2_num_section(int);
unsigned long int	mc2_num_section_norder(int);
unsigned long int	mc2_num_section_norder2(int);
unsigned short int	mc2_paircomp(unsigned short int, unsigned short int);
unsigned long int	choose(unsigned short int, unsigned short int);

/****  Memory allocation functions  ********************************************************/

struct hinge_span	*mc2_newhinge(unsigned int, unsigned int);
struct endhinge	*mc2_newendhinge(unsigned int, unsigned int);
int				*mc2_int_vecalloc(int low, int high);
unsigned int		*mc2_unsgn_vecalloc(int low, int high);
unsigned long int	*mc2_unsgnlong_vecalloc(int low, int high);
vec_ent			*mc2_vecalloc(int low, int high);
mat_ent			**mc2_matalloc(int rowlow, int highrow, int lowcol, int highcol);


/**** Added Functions  ********************************************************************/
void				findsection(double, double, double);
void				mc2_printsection(int);

void 			mc2_fillreordertemplate(int, int);	//fills mc2_reordertemplate appropriately based on mc2_ordertemplate
void 			mc2_fillreordertemplate2(int, int);	//fills mc2_reordertemplate appropriately based on mc2_ordertemplate
void				mc2_printordtemp();
void				mc2_printreordtemp();

void mc2_enterendhinge(int i, int j, int curlength);
void mc2_leaveendhinge(int i, int j, int curlength);
void mc2_endrowedges(int i, int j, int curlength);
void mc2_endcoledges(int i, int j, int curlength);

void mc2_enterendhinge2(int i, int j, int curlength2);
void mc2_leaveendhinge2(int i, int j, int curlength2);
void mc2_endrowedges2(int i, int j, int curlength2);
void mc2_endcoledges2(int i, int j, int curlength2);

unsigned short int mc2_LFlag_endhinge2();
void mc2_recordendtemplate();

void mc2_fillrendtemplate();
void mc2_fillendtemplate2();
void mc2_fillrendtemplate2();

void mc2_fillrendtemplate3();
void mc2_fillendtemplate4();
void mc2_fillrendtemplate4();


unsigned long int mc2_num_section_endhinge_norder(void);
unsigned long int mc2_num_section_endhinge_norder2(void);
unsigned long int mc2_num_section_endhinge_norder3(void);
unsigned long int mc2_num_section_endhinge_norder4(void);

void mc2_printbuiltwalks();
void mc2_printbuiltwalks_all();
void mc2_printbuiltwalks2_all();

void mc2_add_to_built_walks(unsigned long int secnum, int nth_tspan);
void mc2_add_right_endhinge(unsigned long int secnum, int nth_endhinge);

int mc2_reverse_direc(int direc);

void mc2_printtofile();
static int mc2_run_creator_all_after_build(void);

static int mc2_creator_all_mode = 0;
static unsigned long int mc2_creator_limit = 100000UL;
static int mc2_creator_force = 0;

/***************************************************************************/
/******************* end of functions used in this program *****************/
/***************************************************************************/


int run_integrated_2sap_sampler(int argc, char *argv[])
{

    int opt;
    while ((opt = getopt(argc, argv, "L:M:s:n:r:S:m:")) != -1) {
        switch (opt) {
            case 'L': L = atoi(optarg); break;
            case 'M': M = atoi(optarg); break;
            case 's': totalspan = atoi(optarg); break;
            case 'n': samplesize = atoi(optarg); break;
            case 'r': runnum = atoi(optarg); break;
            case 'S': seednum = (unsigned int)atoi(optarg); break;
            case 'm': break;
        }
    }
    mode = 2;
    ham_check = 0;
    mc_2sap_set_system_params(0);
    
    mc2_first_hinge_span = (struct hinge_span**)mc_xcalloc(max_sections + 1, sizeof(struct hinge_span*), "first hinge span table");
    mc2_current_hinge_span = (struct hinge_span**)mc_xcalloc(max_sections + 1, sizeof(struct hinge_span*), "current hinge span table");
    mc2_firstendhinge = (struct endhinge**)mc_xcalloc(max_sections + 1, sizeof(struct endhinge*), "first endhinge table");
    mc2_currentendhinge = (struct endhinge**)mc_xcalloc(max_sections + 1, sizeof(struct endhinge*), "current endhinge table");
    mc2_MC_L_Evector[0] = (double*)mc_xcalloc(max_tspans+1, sizeof(double), "left eigenvector 0");
    mc2_MC_L_Evector[1] = (double*)mc_xcalloc(max_tspans+1, sizeof(double), "left eigenvector 1");
    mc2_MC_R_Evector[0] = (double*)mc_xcalloc(max_tspans+1, sizeof(double), "right eigenvector 0");
    mc2_MC_R_Evector[1] = (double*)mc_xcalloc(max_tspans+1, sizeof(double), "right eigenvector 1");
    mc2_MC_tspans_edges = (unsigned long int**)mc_xcalloc(max_keynum+1, sizeof(unsigned long int*), "2SAP edge table");


	unsigned int seed=seednum;
	mc2_initran_(&seed);

//	time_t t;
//	srand((unsigned) time(&t));

/**** variable used in main ************************************************/

	/*int	           	(*pointordNum)[3]; *//*pointer to an array of int*/
	int				ordNum[6],side = 0;
	int             	i,j,k;	   /*used in for loops */

/**** functions called by main *********************************************/

	clock_t         clock(void);
	void            mc2_conv_to_array(void);
	void            mc2_conv_endhinges_to_array(void);

/**************************************************************************/
/**** start of main program ***********************************************/
/**************************************************************************/

	mc_2sap_ignore_system_result("clear"); /* clears the screen */
	mc_2sap_ignore_system_result("date");  /* prints the date and time */

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
		mc2_current_hinge_span[i] = mc2_newhinge(1,1);		/* free up space of hinge structure */
		mc2_first_hinge_span[i] = mc2_current_hinge_span[i];	/* make mc2_first_hinge_span[i] point to this first one */
		mc2_currentendhinge[i] = mc2_newendhinge(1, 1);		/* free up space of endhinge structure */
		mc2_firstendhinge[i] = mc2_currentendhinge[i];		/* make mc2_firstendhinge[i] point to this first one */
		mc2_sectionkey[i]=0;
	}
	printf("space created for hinge structures, endhinge structures, and mc2_sectionkey initialized\n");

	//initializations//////////////////////////////////////////////////
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			mc2_hingestatus[i][j]=0;
			mc2_alreadyentered[i][j]=0;
			mc2_alreadyentered2[i][j]=0;
		}
	}

	for(i=0; i<=M-1; i++){
		for(j=0; j<=L; j++){
			mc2_colhingeedges[i][j]=0;
		}
	}
	for(i=0; i<=M; i++){
		for(j=0; j<=L-1; j++){
			mc2_rowhingeedges[i][j]=0;
		}
	}

	mc2_num_walks=0;
	for(i=0; i<=vM*vL-1; i++){
		for(j=0; j<=2; j++){
			mc2_curstart[j][i]=-1;
			mc2_curstart2[j][i]=-1;
			mc2_curend[j][i]=-1;
			mc2_curend2[j][i]=-1;
		}
		for(j=0; j<=vM*vL; j++){
			mc2_curwalks[j][i]=0;
			mc2_curwalks2[j][i]=0;
		}
	}
	//initializations end///////////////////////////////////////////////	


	//create 2-spans
	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if( !(i==M && j==L) ){	//don't need to enter at last spot
				printf("2-span: entering at i=%d,j=%d\n", i, j);
				mc2_enterhinge(i, j, side, &ordNum, 0);
			}
			mc2_alreadyentered[i][j]=1;	//don't need to explore any 2-spans that enter/exit here.
		}
	}



/*************************************************************************************************/
	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	newline; /* results are printed to the terminal */
	printf("There are a total of %9lu potential two-span patterns.\n", --mc2_num_2_spans);
	printf("There are a total of %9lu   valid   two-span patterns.\n", mc2_valid_2_spans);
/************************************************************************************************/
	printf("%lu recorded tspans\n", mc2_num_tspans);
	printf("%lu duplicate tspans that were 'valid', but not recorded\n", mc2_dupcounter);

	printf("\n");


////////////////////////////////////////////////////////////////////////


	//generate the endhinges
	//we are going to create "right endhinges" (entering from the left, can't leave to the right)
	//When we record these right endhinges, then we will make appropriate changes to also record it's "mirror" left endhinge

	//initialize stuff to zero.
	mc2_num_walks=0;
	mc2_num_walks2=0;
	for(i=0; i<=vM*vL-1; i++){
		for(j=0; j<=2; j++){
			mc2_curstart[j][i]=-1;
			mc2_curend[j][i]=-1;
			mc2_curstart2[j][i]=-1;
			mc2_curend2[j][i]=-1;
		}
		for(j=0; j<=vM*vL; j++){
			mc2_curwalks[j][i]=0;
			mc2_curwalks2[j][i]=0;
		}
	}

	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			mc2_hingestatus[i][j]=0;
			mc2_alreadyentered[i][j]=0;
			mc2_alreadyentered2[i][j]=0;
		}
	}

	for(i=0; i<=M-1; i++){
		for(j=0; j<=L; j++){
			mc2_colhingeedges[i][j]=0;
		}
	}
	for(i=0; i<=M; i++){
		for(j=0; j<=L-1; j++){
			mc2_rowhingeedges[i][j]=0;
		}
	}

	mc2_EndOrdNum[0]=1;	//first edge in section will get the label 1
	mc2_EndOrdNum[1]=0;	//no edges currently in the endhinge
	mc2_EndOrdNum2[0]=1;	//first edge in section will get the label 1
	mc2_EndOrdNum2[1]=0;	//no edges currently in the endhinge

	//create endhinges
	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if( !(i==M && j==L) ){	//don't need to enter at last spot
				printf("Endhinge: entering at i=%d,j=%d\n", i, j);
				mc2_enterendhinge(i, j, 0);
			}
			mc2_alreadyentered[i][j]=1;	//don't need to explore any 2-spans that enter/exit here.
		}
	}

	printf("number of potential endhinges: %lu\n", mc2_num_potential_endhinges);
	printf("number of valid endhinges: %lu\n", mc2_num_valid_endhinges);

	printf("%lu right endhinges recorded\n", mc2_tot_right_endhinges);
	printf("%lu left endhinges recorded\n", mc2_tot_left_endhinges);
	printf("%lu non-lexi left endhinges\n", mc2_tot_left_endhinges);
	printf("%lu duplicate right endhinges\n", mc2_num_duplicate_right_endhinges);
	printf("%lu duplicate left endhinges\n", mc2_num_duplicate_left_endhinges);




//tests
/*
	for(i=1; i<=max_sections; i++){
		mc2_currentendhinge[i]=mc2_firstendhinge[i];
		while(mc2_currentendhinge[i]->nextendhinge!=NULL){
			mc2_currentendhinge[i]=mc2_currentendhinge[i]->nextendhinge;
			if(mc2_currentendhinge[i]->side==1){
				printf("Endhinge with mc2_sectionkey =%d,%lu. side=%d, mc2_num_walks=%d, mc2_num_walks2=%d\n", i, mc2_currentendhinge[i]->sec2, mc2_currentendhinge[i]->side, mc2_currentendhinge[i]->num_of_walks, mc2_currentendhinge[i]->num_of_walks2);
				for(j=0; j<=mc2_currentendhinge[i]->num_of_walks-1; j++){
					printf("walk1: %d. Start=(%d,%d,%d), Directions: ", j, mc2_currentendhinge[i]->start[0][j], mc2_currentendhinge[i]->start[1][j], mc2_currentendhinge[i]->start[2][j]);
					for(k=0; k<=vM*vL; k++){
						printf("%d ", mc2_currentendhinge[i]->walks[k][j]);
					}
					printf(". End=(%d,%d,%d)\n", mc2_currentendhinge[i]->end[0][j], mc2_currentendhinge[i]->end[1][j], mc2_currentendhinge[i]->end[2][j]);
				}
				for(j=0; j<=mc2_currentendhinge[i]->num_of_walks2-1; j++){
					printf("walk2: %d. Start=(%d,%d,%d), Directions: ", j, mc2_currentendhinge[i]->start2[0][j], mc2_currentendhinge[i]->start2[1][j], mc2_currentendhinge[i]->start2[2][j]);
					for(k=0; k<=vM*vL; k++){
						printf("%d ", mc2_currentendhinge[i]->walks2[k][j]);
					}
					printf(". End=(%d,%d,%d)\n", mc2_currentendhinge[i]->end2[0][j], mc2_currentendhinge[i]->end2[1][j], mc2_currentendhinge[i]->end2[2][j]);
				}
				printf("hedges = ");
				for(j=0; j<=vM*vL-1; j++){
					printf("%d ", mc2_currentendhinge[i]->hedges[j]);
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
	mc2_conv_to_array();
	printf("done converting hinge_spans to arrays.\n\n");

	printf("converting endhinge linked lists to arrays...\n");
	mc2_conv_endhinges_to_array();
	printf("done converting endhinges to arrays.\n\n");



//tests

/*
	int l;
	for(i=1; i<=max_keynum; i++){
		for(j=1; j<=mc2_num_outsections[i]; j++){
			printf("\nsectionkey %d (%lu,%lu) connects to mc2_sectionkey %lu (%lu, %lu). It consists of %d walks in poly1, and %d walks in poly2.\n", i, mc2_sectionkey2SAP[i][0], mc2_sectionkey2SAP[i][1], mc2_t_outsection[i][j], mc2_sectionkey2SAP[mc2_t_outsection[i][j]][0], mc2_sectionkey2SAP[mc2_t_outsection[i][j]][1], mc2_t_num_walks[i][j], mc2_t_num_walks2[i][j]);
			printf("poly1:\n");
			for(k=0; k<=mc2_t_num_walks[i][j]-1; k++){
				printf("Walk %d. Start (%d,%d,%d). Directions: ", k, mc2_t_start[i][j][0][k], mc2_t_start[i][j][1][k], mc2_t_start[i][j][2][k]);
				for(l=0; l<=vM*vL; l++){
					printf("%d ", mc2_t_walks[i][j][l][k]);
				}
				printf(". End (%d,%d,%d)\n", mc2_t_end[i][j][0][k], mc2_t_end[i][j][1][k], mc2_t_end[i][j][2][k]);
			}
			printf("poly2:\n");
			for(k=0; k<=mc2_t_num_walks2[i][j]-1; k++){
				printf("Walk %d. Start (%d,%d,%d). Directions: ", k, mc2_t_start2[i][j][0][k], mc2_t_start2[i][j][1][k], mc2_t_start2[i][j][2][k]);
				for(l=0; l<=vM*vL; l++){
					printf("%d ", mc2_t_walks2[i][j][l][k]);
				}
				printf(". End (%d,%d,%d)\n", mc2_t_end2[i][j][0][k], mc2_t_end2[i][j][1][k], mc2_t_end2[i][j][2][k]);
			}

		}
	}


	for(i=1; i<=max_keynum; i++){
		printf("mc2_sectionkey %d has %lu left endhinges attached\n", i, mc2_num_left_endhinges[i]);
		for(j=1; j<=mc2_num_left_endhinges[i]; j++){
			printf("left endhinge %d has %d walk(s) in poly1 and %d walk(s) in poly2\n", j, mc2_Lend_num_walks[i][j], mc2_Lend_num_walks2[i][j]);
			printf("poly1:\n");
			for(k=0; k<=mc2_Lend_num_walks[i][j]-1; k++){
				printf("Walk %d. Start (%d,%d,%d). Directions: ", k, mc2_Lend_start[i][j][0][k], mc2_Lend_start[i][j][1][k], mc2_Lend_start[i][j][2][k]);
				for(l=0; l<=vM*vL; l++){
					printf("%d ", mc2_Lend_walks[i][j][l][k]);
				}
				printf(". End (%d,%d,%d)\n", mc2_Lend_end[i][j][0][k], mc2_Lend_end[i][j][1][k], mc2_Lend_end[i][j][2][k]);
			}
			printf("poly2:\n");
			for(k=0; k<=mc2_Lend_num_walks2[i][j]-1; k++){
				printf("Walk %d. Start (%d,%d,%d). Directions: ", k, mc2_Lend_start2[i][j][0][k], mc2_Lend_start2[i][j][1][k], mc2_Lend_start2[i][j][2][k]);
				for(l=0; l<=vM*vL; l++){
					printf("%d ", mc2_Lend_walks2[i][j][l][k]);
				}
				printf(". End (%d,%d,%d)\n", mc2_Lend_end2[i][j][0][k], mc2_Lend_end2[i][j][1][k], mc2_Lend_end2[i][j][2][k]);
			}
		}
	}

*/






	printf("Finished generating all 2-spans and endhinges.\n");

	if (mc2_creator_all_mode) {
		return mc2_run_creator_all_after_build();
	}

	if(totalspan<2){
		printf("totalspan=%d. This program only works for span>=2. Exitting\n", totalspan);
		exit(1);
	}

	printf("\nNOW SAMPLING 2SAPs: %d samples from L=%d, M=%d, span=%d\n", samplesize, L, M, totalspan);


	mkdir("data", 0775);
	mkdir("data/MonteCarlo", 0775);
	mkdir("data/MonteCarlo/2SAPs", 0775);
	mc_checked_snprintf(mc2_filename, sizeof(mc2_filename), "data/MonteCarlo/2SAPs/MC2SAPsL%dM%dspan%drun%dnum%lu.txt", L, M, totalspan, runnum, mc2_filenum);
	mc2_fp = fopen(mc2_filename, "w");	//create or overwrite the file "mc2_filename"

	if(mc2_fp != NULL){
		fprintf(mc2_fp, "UofS\n");	//first line in file is always "UofS"
		run_metadata_write(mc2_filename, "mc_master", "samples_uofs", mode, L, M, totalspan, seednum, runnum, dom_evalue);
		printf("printed UofS in file '%s'\n", mc2_filename);
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
//		exit(EXIT_FAILURE);
//	}

//	static double R_Evector[max_tspans+1];

	double* R_Evector;

    for(int k=1; k<=max_keynum; k++){
        if(mc2_num_outsections[k] == 0) continue;
        mc2_MC_tspans_edges[k] = mc2_unsgnlong_vecalloc(1, mc2_num_outsections[k]);
        for(int m=1; m<=mc2_num_outsections[k]; m++) mc2_MC_tspans_edges[k][m] = 0;
    }
    double calculated_dom_evalue = mc2_max_eval_LRvec(1.0) + 1.0;
    printf("Calculated 2SAP dominant eigenvalue: %.15f (Expected: %.15f)\n", calculated_dom_evalue, dom_evalue);
    dom_evalue = calculated_dom_evalue;

    mkdir("data", 0775);
    mkdir("data/MonteCarlo", 0775);
    mkdir("data/MonteCarlo/MC_Evectors", 0775);
    char export_fn[128];
    mc_checked_snprintf(export_fn, sizeof(export_fn), "data/MonteCarlo/MC_Evectors/2SAP_R_Evector%s_TS_L%dM%d.txt", (ham_check ? "Ham" : ""), L, M);
    FILE *export_fp = fopen(export_fn, "w");
    if (export_fp != NULL) {
        for (int i = 1; i <= max_tspans; i++) {
            fprintf(export_fp, "%.15f\n", mc2_MC_R_Evector[0][i]);
        }
        fclose(export_fp);
        run_metadata_write(export_fn, "mc_master", "right_eigenvector", mode, L, M, totalspan, seednum, runnum, calculated_dom_evalue);
        printf("Calculated eigenvectors exported to %s\n", export_fn);
    }

	R_Evector = mc2_MC_R_Evector[0];

/*	for(i=1; i<= max_tspans; i++){
//		printf("L_Evector[%d] = %.15f\n", i, L_Evector[i]);
		printf("R_Evector[%d] = %.15f\n", i, R_Evector[i]);
	}
*/


	double maxt_one=0.0;	//used in rejection
	double t_one[mc2_tot_left_endhinges];	//each left end hinge has a t_one value.
	curLEH=0;
	for(i=1; i<=max_keynum; i++){
		for(j=1; j<=mc2_num_left_endhinges[i]; j++){
			t_one[curLEH]=0;
			for(k=1; k<=mc2_num_outsections[i]; k++){
				t_one[curLEH] = t_one[curLEH] + R_Evector[mc2_t_nrr[i][k]];
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
	t_two = (double*)mc_xcalloc((size_t)max_tspans + 1, sizeof(*t_two), "2SAP second-step weights");


	for(i=1; i<=max_keynum; i++){
		for(j=1; j<=mc2_num_outsections[i]; j++){
//			printf("keynum=%d, outsection_ent=%d, mc2_num_right_endhinges[%lu]=%lu, R_Evector[%lu]=%f\n", i, j, mc2_t_outsection[i][j], mc2_num_right_endhinges[mc2_t_outsection[i][j]], mc2_t_nrr[i][j], R_Evector[mc2_t_nrr[i][j]]);
			if(R_Evector[mc2_t_nrr[i][j]]>0.000000000000001){
				t_two[mc2_t_nrr[i][j]] = mc2_num_right_endhinges[mc2_t_outsection[i][j]] / R_Evector[mc2_t_nrr[i][j]];
			}
			else{
				t_two[mc2_t_nrr[i][j]]=0;
			}
			if(t_two[mc2_t_nrr[i][j]] > maxt_two){
				maxt_two = t_two[mc2_t_nrr[i][j]];
			}
//			printf("filled t_two[%lu]=%f. maxt_two=%f\n", mc2_t_nrr[i][j], t_two[mc2_t_nrr[i][j]], maxt_two);
		}
	}
	unsigned long int reject_two=0;




//	double chance=0.0;
	while(curSample<=samplesize-1){

		//initialize built stuff to zeros or -1
		for(i=0; i<=vM*vL/2-1; i++){
			for(j=0; j<=2; j++){
				mc2_built_walks_start[i][j]=-1;
				mc2_built_walks_end[i][j]=-1;
				mc2_built_walks_start2[i][j]=-1;
				mc2_built_walks_end2[i][j]=-1;
			}
			for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
				mc2_built_walks_direcs[i][j]=0;
				mc2_built_walks_direcs2[i][j]=0;
			}
		}

		chosenLEH = floor(mc2_ran1real_()*mc2_tot_left_endhinges); //random number [0,mc2_tot_left_endhinges-1]

//	printf("\ntot_left_endhinges=%lu\n", mc2_tot_left_endhinges);
//	chance = (double)1/mc2_tot_left_endhinges;
//	printf("chose left endhinge.\n");
//	printf("numerator=%d. denominator=%d\n", 1, mc2_tot_left_endhinges);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);

//		printf("\nSample=%lu. chosenLEH=%d, from range [%d,%lu].\n", curSample, chosenLEH, 0, mc2_tot_left_endhinges-1);
		curLEH=0;
		for(secnum=1; secnum<=max_keynum; secnum++){
			for(nth_endhinge=1; nth_endhinge<=mc2_num_left_endhinges[secnum]; nth_endhinge++){
//				printf("curLEH=%d. secnum=%lu, nth_endhinge=%d\n", curLEH, secnum, nth_endhinge);
				if(curLEH==chosenLEH){
//					printf("chosenLEH=%d was found. Section=%lu, nth_endhinge=%d\n", chosenLEH, secnum, nth_endhinge);

					//CHECK FOR REJECTION
//					printf("accept prob_one=%f\n", t_one[curLEH] / maxt_one);
					if( mc2_ran1real_() < t_one[curLEH] / maxt_one ){	//do not reject
//	chance = chance * t_one[curLEH] / maxt_one;
//	printf("passed first rejection with prob=%f.\n", t_one[curLEH] / maxt_one);
//	printf("numerator=%f. denominator=%f\n", t_one[curLEH], maxt_one);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);

						//fill built_walks with left endhinge
						mc2_num_built_walks=mc2_Lend_num_walks[secnum][nth_endhinge];
						mc2_num_built_walks2=mc2_Lend_num_walks2[secnum][nth_endhinge];
						for(nth_walk=0; nth_walk<=mc2_num_built_walks-1; nth_walk++){
							for(i=0; i<=2; i++){
								mc2_built_walks_start[nth_walk][i] = mc2_Lend_start[secnum][nth_endhinge][i][nth_walk];
								mc2_built_walks_end[nth_walk][i] = mc2_Lend_end[secnum][nth_endhinge][i][nth_walk];
							}
							for(i=0; i<=vM*vL; i++){
								mc2_built_walks_direcs[nth_walk][i] = mc2_Lend_walks[secnum][nth_endhinge][i][nth_walk];
							}
							for(i=vM*vL+1; i<=vM*vL*(totalspan+1)-1; i++){
								mc2_built_walks_direcs[nth_walk][i]=0;
							}
						}
						for(nth_walk=0; nth_walk<=mc2_num_built_walks2-1; nth_walk++){
							for(i=0; i<=2; i++){
								mc2_built_walks_start2[nth_walk][i] = mc2_Lend_start2[secnum][nth_endhinge][i][nth_walk];
								mc2_built_walks_end2[nth_walk][i] = mc2_Lend_end2[secnum][nth_endhinge][i][nth_walk];
							}
							for(i=0; i<=vM*vL; i++){
								mc2_built_walks_direcs2[nth_walk][i] = mc2_Lend_walks2[secnum][nth_endhinge][i][nth_walk];
							}
							for(i=vM*vL+1; i<=vM*vL*(totalspan+1)-1; i++){
								mc2_built_walks_direcs2[nth_walk][i]=0;
							}
						}
						//built_walks is filled with a left endhinge.
						curspan=1;
//						printf("\nLeft endhinge has been added to built_walks:\n");
//						mc2_printbuiltwalks_all();
//						mc2_printbuiltwalks2_all();
//						printf("\n");

						sec1=secnum;

						//choose first tspan with probability proportional to it's right eigenvector value (from those t-spans who have first section = sec1)
						prob = mc2_ran1real_();
//						printf("prob=%f\n", prob);
						nth_tspan=1;
						sumofprobs = R_Evector[mc2_t_nrr[sec1][nth_tspan]] / t_one[curLEH];
						while(sumofprobs<prob){
							nth_tspan++;
							sumofprobs = sumofprobs + R_Evector[mc2_t_nrr[sec1][nth_tspan]] / t_one[curLEH];
						}
//	chance = chance * R_Evector[mc2_t_nrr[sec1][nth_tspan]] / t_one[curLEH];
//	printf("first block chosen with prob=%f. It was the %d-th choice (out of %lu choices). It consists of sections %lu->%lu (tspan number %lu).\n", R_Evector[mc2_t_nrr[sec1][nth_tspan]] / t_one[curLEH], nth_tspan, mc2_num_outsections[sec1], sec1, mc2_t_outsection[sec1][nth_tspan], mc2_t_nrr[sec1][nth_tspan]);
//	printf("numerator=%f. denominator=%f\n", R_Evector[mc2_t_nrr[sec1][nth_tspan]], t_one[curLEH]);
//	printf("chance = %f\n", chance);
//	printf("1/chance = %f\n\n", 1/chance);

						mc2_add_to_built_walks(sec1, nth_tspan);	//add nth_tspan from sec1 to built_walks

//							printf("first block has been added, here is the result:\n");
//							mc2_printbuiltwalks_all();
//							mc2_printbuiltwalks2_all();
//							printf("\n");


						sec2 = mc2_t_outsection[sec1][nth_tspan];	//set sec2
						cur_tspan_num = mc2_t_nrr[sec1][nth_tspan];
						if(R_Evector[cur_tspan_num] > 0.000000000000001){
							curspan++;

							while(curspan<totalspan){
//									printf("adding another tspan to section %lu. Probabilities:\n", sec2);
								sumofprobs=0.0;
								for(i=1; i<=mc2_num_outsections[sec2]; i++){
									sumofprobs = sumofprobs + R_Evector[mc2_t_nrr[sec2][i]] / R_Evector[cur_tspan_num] / dom_evalue;
//										printf("sec %lu -> sec %lu. tspan_num=%lu. R_Evector[%lu]=%f. probability=%f. cumulative prob=%f\n", sec2, mc2_t_outsection[sec2][i], mc2_t_nrr[sec2][i], mc2_t_nrr[sec2][i], R_Evector[mc2_t_nrr[sec2][i]], R_Evector[mc2_t_nrr[sec2][i]] / R_Evector[cur_tspan_num] / dom_evalue, sumofprobs);
								}
//								printf("sumofprobs=%f\n", sumofprobs);
								if(sumofprobs<0.9 || sumofprobs>1.1){
									printf("second sumofprobs=%f\n", sumofprobs);
									printf("problem with second sumofprobs when leaving tspan number %lu, (made up of %lu and %lu)\n", cur_tspan_num, sec1, sec2);
									printf("This was on sample %lu. Exitting\n", curSample);
									exit(1);
								}
	

								//choose next tspan according to correct probability
								prob = mc2_ran1real_();
								prob = prob * sumofprobs; //in case sum of probs doesn't add up to one.
	//							printf("prob=%f\n", prob);
								nth_tspan=1;
								sumofprobs = R_Evector[mc2_t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue;
								while(sumofprobs<prob){
									nth_tspan++;
									sumofprobs = sumofprobs + R_Evector[mc2_t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue;
								}
	//	chance = chance * R_Evector[mc2_t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue;
//		printf("next block chosen with prob=%f. It is sections %lu->%lu (tspan number %lu)\n", R_Evector[mc2_t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue, sec2, mc2_t_outsection[sec2][nth_tspan], mc2_t_nrr[sec2][nth_tspan]);
	//	printf("numerator=%f. denominator=%f. dom_evalue=%f.\n", R_Evector[mc2_t_nrr[sec2][nth_tspan]], R_Evector[cur_tspan_num], dom_evalue);
	//	printf("chance = %f\n", chance);
	//	printf("1/chance = %f\n\n", 1/chance);
	//							printf("sec %lu -> sec %lu. tspan_num=%lu. probability=%f. cumulative prob=%f\n", sec2, mc2_t_outsection[sec2][nth_tspan], mc2_t_nrr[sec2][nth_tspan], R_Evector[mc2_t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue, sumofprobs);
								sec1=sec2;
								mc2_add_to_built_walks(sec1, nth_tspan);	//add nth_tspan from sec2 to built_walks

//								printf("tspan has added, here is the result:\n");
//								mc2_printbuiltwalks_all();
//								mc2_printbuiltwalks2_all();
//								printf("\n");
								sec2 = mc2_t_outsection[sec1][nth_tspan];
								cur_tspan_num = mc2_t_nrr[sec1][nth_tspan];
								if(R_Evector[cur_tspan_num] < 0.000000000000001){
									printf("tspan chosen (during intermediate tspans) that will not result in a poly. Shouldn't ever happen. Exitting\n");
									printf("cur_tspan_num=%lu, sections %lu->%lu (%d-th outsec).\n", cur_tspan_num, sec1, mc2_t_outsection[sec1][nth_tspan], nth_tspan);
									exit(1);
								}
								curspan++;
							}
							//CHECK FOR REJECTION
//									printf("accept prob_two=%f\n", t_two[cur_tspan_num] / maxt_two);
							if( mc2_ran1real_() < t_two[cur_tspan_num] / maxt_two ){	//do not reject
	//	chance = chance * t_two[cur_tspan_num] / maxt_two;
//				printf("passed 2nd rejection with prob=%f.\n", t_two[cur_tspan_num] / maxt_two);
	//	printf("numerator=%f. denominator=%f\n", t_two[cur_tspan_num], maxt_two);
	//	printf("chance = %f\n", chance);
	//	printf("1/chance = %f\n\n", 1/chance);
							chosenREH = 1 + floor(mc2_ran1real_()*mc2_num_right_endhinges[sec2]); //rando number [1, mc2_num_right_endhinges[sec2]]
	//	chance = chance * (double)1/mc2_num_right_endhinges[sec2];
//				printf("chose final endhinge with prob=%f.\n", (double)1/mc2_num_right_endhinges[sec2]);
	//	printf("chance = %f\n", chance);
	//	printf("1/chance = %f\n\n", 1/chance);
//										printf("closing off SAP with endhinge: section=%lu, nth_endhinge=%d (possible [%d,%lu]).\n", sec2, chosenREH, 1, mc2_num_right_endhinges[sec2]);
								mc2_add_right_endhinge(sec2, chosenREH);

		//						printf("FINAL POLY:\n");
		//						mc2_printbuiltwalks();

								//sample created
	//	printf("FINAL CHANCE OF CHOOSING:\n");
	//	printf("chance = %f\n", chance);
	//	printf("1/chance = %f\n", 1/chance);
								curSample++;

								//PRINT TO FILE.
								mc2_printtofile();
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


	fprintf(mc2_fp, "-999\n");
	fclose(mc2_fp);
	printf("File %lu contains %lu polys.\n", mc2_filenum, mc2_filetotal);

	printf("\nSampling Complete.\n");

	printf("\nRESULTS: %lu samples created from L=%d, M=%d, span=%d\n", curSample, L, M, totalspan);
	printf("Rejected at first step %lu times\n", reject_one);
	printf("Rejected at last step %lu times\n", reject_two);


	printf("\n");	



	mc_2sap_ignore_system_result("date");  /* prints the date and time */
	printf("\nProgram Complete.\n\n");

	return 0;

}				/* end of main */
/***************************************************************************/


/***** start of function mc2_conv_to_array *************************************/
/* called only by main */


void mc2_conv_to_array(void)
{
	//Fill mc2_sectionkey2SAP[key][2]
		//mc2_sectionkey2SAP[key][0] = secnum of poly1
		//mc2_sectionkey2SAP[key][1] = secnum of poly2
		//key=key number assigned to this ordered pairing.
	//Also fill mc2_num_outsections at the same time.
	int i, j;
	int found;

	for(i=1; i<=max_keynum; i++){
		mc2_sectionkey2SAP[i][0]=0;
		mc2_sectionkey2SAP[i][1]=0;
		mc2_num_outsections[i]=0;
	}

	int keycounter=0;

	for(i=1; i<=max_sections; i++){
		mc2_current_hinge_span[i] = mc2_first_hinge_span[i];
		while(mc2_current_hinge_span[i]->nexthinge != NULL) {
			mc2_current_hinge_span[i] = mc2_current_hinge_span[i]->nexthinge;
//			printf("looking at section composed of sections %d and %lu (on the left)\n", i, mc2_current_hinge_span[i]->inorder2);
			j=1;
			found=0;
			while(mc2_sectionkey2SAP[j][0]!=0){
				if(mc2_sectionkey2SAP[j][0]==i && mc2_sectionkey2SAP[j][1]==mc2_current_hinge_span[i]->inorder2){
//					printf("combo found in key already\n");
					found=1;
					mc2_num_outsections[j]++;
					break;
				}
				j++;
			}
			if(found==0){	//need to record sections into mc2_sectionkey2SAP
//				printf("combo not in key yet. adding: j=%d, sec1=%d, sec2=%lu\n", j, i, mc2_current_hinge_span[i]->inorder2);
				keycounter++;
				mc2_sectionkey2SAP[j][0] = i;
				mc2_sectionkey2SAP[j][1] = mc2_current_hinge_span[i]->inorder2;
				mc2_num_outsections[j]++;
			}	
		}
	}

	printf("Total of %d 2SAP sections\n", keycounter);

/*	for(i=1; i<=max_keynum; i++){
		printf("keynum=%d, sections: %lu,%lu (OG: %lu,%lu)\n", i, mc2_sectionkey2SAP[i][0], mc2_sectionkey2SAP[i][1], mc2_sectionkey[mc2_sectionkey2SAP[i][0]], mc2_sectionkey[mc2_sectionkey2SAP[i][1]]);
		mc2_printsection(mc2_sectionkey[mc2_sectionkey2SAP[i][0]]);
		mc2_printsection(mc2_sectionkey[mc2_sectionkey2SAP[i][1]]);
		printf("\n");
	}
*/

	for(i=1; i<=max_keynum; i++){
		mc2_t_outsection[i] = mc2_unsgnlong_vecalloc(1, mc2_num_outsections[i]); //allocate memory
		mc2_t_num_walks[i] = mc2_unsgn_vecalloc(1, mc2_num_outsections[i]);
		mc2_t_num_walks2[i] = mc2_unsgn_vecalloc(1, mc2_num_outsections[i]);
		mc2_t_nrr[i] = mc2_unsgnlong_vecalloc(1, mc2_num_outsections[i]);
	}



	unsigned long int section_num;
	unsigned long int sec1key;
	unsigned long int sec2key;
	unsigned int curArrayEnt[MAX_KEYNUM_ARR+1];
	for(i=1; i<=max_keynum; i++){
		curArrayEnt[i]=1;
	}

	struct hinge_span *hinge_to_free;


	//fill mc2_t_outsection and mc2_t_num_walks and mc2_t_num_walks2
	for(section_num = 1; section_num <= max_sections; section_num++) {
		mc2_current_hinge_span[section_num] = mc2_first_hinge_span[section_num];
		while ((*mc2_current_hinge_span[section_num]).nexthinge != NULL) {
			mc2_current_hinge_span[section_num] = (*mc2_current_hinge_span[section_num]).nexthinge;

//			printf("this hinge is made up of sections %lu,%lu -> %lu,%lu\n", section_num, (*mc2_current_hinge_span[section_num]).inorder2, (*mc2_current_hinge_span[section_num]).outorder, (*mc2_current_hinge_span[section_num]).outorder2);

			//get keynums
			sec1key=0;
			for(i=1; i<=max_keynum; i++){
//				printf("keynum=%d, match: %lu,%lu ?\n", i, mc2_sectionkey2SAP[i][0], mc2_sectionkey2SAP[i][1]);
				if(mc2_sectionkey2SAP[i][0]==section_num && mc2_sectionkey2SAP[i][1]==(*mc2_current_hinge_span[section_num]).inorder2){
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
				if(mc2_sectionkey2SAP[i][0]==(*mc2_current_hinge_span[section_num]).outorder && mc2_sectionkey2SAP[i][1]==(*mc2_current_hinge_span[section_num]).outorder2){
//					printf("sec2key=%d\n", i);
					sec2key=i;
					break;
				}
			}
			if(sec2key==0){
				printf("error: didn't find sec2key. Sections were %lu,%lu. Exitting\n", (*mc2_current_hinge_span[section_num]).outorder, (*mc2_current_hinge_span[section_num]).outorder2);
				exit(1);
			}
//			printf("sec1key=%lu, sec2key=%lu\n", sec1key, sec2key);
			mc2_t_outsection[sec1key][curArrayEnt[sec1key]] = sec2key;
			mc2_t_num_walks[sec1key][curArrayEnt[sec1key]] = (*mc2_current_hinge_span[section_num]).num_of_walks;
			mc2_t_num_walks2[sec1key][curArrayEnt[sec1key]] = (*mc2_current_hinge_span[section_num]).num_of_walks2;
			curArrayEnt[sec1key]++;
		}
	}



	//allocate space for mc2_t_start, mc2_t_end, and mc2_t_walks -- AND -- mc2_t_start2, mc2_t_end2, and mc2_t_walks2
	printf("allocating space for mc2_t_start, mc2_t_end, and mc2_t_walks AND mc2_t_start2, mc2_t_end2, and mc2_t_walks2\n");

	for(section_num=1; section_num<=max_keynum; section_num++){
		mc2_t_start[section_num] = mc_2sap_alloc_int3_table(mc2_num_outsections[section_num], 3, mc2_t_num_walks[section_num], "2SAP transition starts");
		mc2_t_start2[section_num] = mc_2sap_alloc_int3_table(mc2_num_outsections[section_num], 3, mc2_t_num_walks2[section_num], "2SAP transition starts2");
		mc2_t_end[section_num] = mc_2sap_alloc_int3_table(mc2_num_outsections[section_num], 3, mc2_t_num_walks[section_num], "2SAP transition ends");
		mc2_t_end2[section_num] = mc_2sap_alloc_int3_table(mc2_num_outsections[section_num], 3, mc2_t_num_walks2[section_num], "2SAP transition ends2");
		mc2_t_walks[section_num] = mc_2sap_alloc_int3_table(mc2_num_outsections[section_num], vM * vL + 1, mc2_t_num_walks[section_num], "2SAP transition walks");
		mc2_t_walks2[section_num] = mc_2sap_alloc_int3_table(mc2_num_outsections[section_num], vM * vL + 1, mc2_t_num_walks2[section_num], "2SAP transition walks2");
	}
	printf("Done allocating space for mc2_t_start, mc2_t_end, and mc2_t_walks AND mc2_t_start2, mc2_t_end2, and mc2_t_walks2\n");



	///Fill mc2_t_start, mc2_t_end, mc2_t_walks -- AND -- mc2_t_start2, mc2_t_end2, mc2_t_walks2
	printf("Filling mc2_t_start, mc2_t_end, and mc2_t_walks AND mc2_t_start2, mc2_t_end2, and mc2_t_walks2\n");
	for(i=1; i<=max_keynum; i++){
		curArrayEnt[i]=1;
	}
	for (section_num = 1; section_num <= max_sections; section_num++) {
//		printf("section_num=%lu\n", section_num);
		mc2_current_hinge_span[section_num] = mc2_first_hinge_span[section_num];
		while ((*mc2_current_hinge_span[section_num]).nexthinge != NULL) {
			hinge_to_free = mc2_current_hinge_span[section_num];
			mc2_current_hinge_span[section_num] = (*mc2_current_hinge_span[section_num]).nexthinge;
			free(hinge_to_free);

			//get keynums
			sec1key=0;
			for(i=1; i<=max_keynum; i++){
//				printf("keynum=%d, match: %lu,%lu ?\n", i, mc2_sectionkey2SAP[i][0], mc2_sectionkey2SAP[i][1]);
				if(mc2_sectionkey2SAP[i][0]==section_num && mc2_sectionkey2SAP[i][1]==(*mc2_current_hinge_span[section_num]).inorder2){
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
				if(mc2_sectionkey2SAP[i][0]==(*mc2_current_hinge_span[section_num]).outorder && mc2_sectionkey2SAP[i][1]==(*mc2_current_hinge_span[section_num]).outorder2){
//					printf("sec2key=%d\n", i);
					sec2key=i;
					break;
				}
			}
			if(sec2key==0){
				printf("error: didn't find sec2key. Sections were %lu,%lu. Exitting\n", (*mc2_current_hinge_span[section_num]).outorder, (*mc2_current_hinge_span[section_num]).outorder2);
				exit(1);
			}
//			printf("\nsec1key=%lu, sec2key=%lu\n", sec1key, sec2key);



			for(i=0; i<=(mc2_current_hinge_span[section_num]->num_of_walks)-1; i++){
				for(j=0; j<=2; j++){
//					printf("recording %d-th coordinate of %d-th walk\n", j, i);
//					printf("keynum=%lu, arr_ent=%lu, j=%d, i=%d\n", sec1key, curArrayEnt[sec1key], j, i);
					mc2_t_start[sec1key][curArrayEnt[sec1key]][j][i] = mc2_current_hinge_span[section_num]->start[j][i];
//					printf("mc2_t_start[%lu][%lu][%d][%d] = %d\n", sec1key, curArrayEnt[sec1key], j, i, mc2_t_start[sec1key][curArrayEnt[sec1key]][j][i]);
					mc2_t_end[sec1key][curArrayEnt[sec1key]][j][i] = mc2_current_hinge_span[section_num]->end[j][i];
//					printf("mc2_t_end[%lu][%lu][%d][%d] = %d\n", sec1key, curArrayEnt[sec1key], j, i, mc2_t_end[sec1key][curArrayEnt[sec1key]][j][i]);
				}
				for(j=0; j<=vM*vL; j++){
					mc2_t_walks[sec1key][curArrayEnt[sec1key]][j][i] = mc2_current_hinge_span[section_num]->walks[j][i];
				}
			}

			for(i=0; i<=(mc2_current_hinge_span[section_num]->num_of_walks2)-1; i++){
				for(j=0; j<=2; j++){
//					printf("recording %d-th coordinate\n", j);
//					printf("keynum=%lu, arr_ent=%lu, j=%d, i=%d\n", sec1key, curArrayEnt[sec1key], j, i);
					mc2_t_start2[sec1key][curArrayEnt[sec1key]][j][i] = mc2_current_hinge_span[section_num]->start2[j][i];
//					printf("mc2_t_start2[%lu][%lu][%d][%d] = %d\n", sec1key, curArrayEnt[sec1key], j, i, mc2_t_start2[sec1key][curArrayEnt[sec1key]][j][i]);
					mc2_t_end2[sec1key][curArrayEnt[sec1key]][j][i] = mc2_current_hinge_span[section_num]->end2[j][i];
//					printf("mc2_t_end2[%lu][%lu][%d][%d] = %d\n", sec1key, curArrayEnt[sec1key], j, i, mc2_t_end2[sec1key][curArrayEnt[sec1key]][j][i]);
				}
				for(j=0; j<=vM*vL; j++){
					mc2_t_walks2[sec1key][curArrayEnt[sec1key]][j][i] = mc2_current_hinge_span[section_num]->walks2[j][i];
				}
			}
			curArrayEnt[sec1key]++;

		}
	}
	printf("Done filling mc2_t_start, mc2_t_end, and mc2_t_walks\n");





	//fill mc2_t_nrr
	unsigned long int tspan_num=0;
	unsigned long int nth_outsection;
	for (section_num=1; section_num<=max_keynum;section_num++){
		for (nth_outsection=1;nth_outsection<=mc2_num_outsections[section_num];nth_outsection++){
			mc2_t_nrr[section_num][nth_outsection] = ++tspan_num;
		}
	}

/*
	int counter=1;
	for (section_num=1; section_num<=max_keynum;section_num++){
		for (nth_outsection=1;nth_outsection<=mc2_num_outsections[section_num];nth_outsection++){
			printf("This tspan (%d=%lu) consists of keys %lu->%lu, which are pairs of section keys (%lu,%lu)->(%lu,%lu), which have OG section_nums (%lu, %lu)->(%lu, %lu)\n", counter, tspans_nrr[section_num][nth_outsection], section_num, tspans_outsection[section_num][nth_outsection], mc2_sectionkey2SAP[0][section_num], mc2_sectionkey2SAP[1][section_num], mc2_sectionkey2SAP[0][tspans_outsection[section_num][nth_outsection]], mc2_sectionkey2SAP[1][tspans_outsection[section_num][nth_outsection]], mc2_sectionkey[mc2_sectionkey2SAP[0][section_num]], mc2_sectionkey[mc2_sectionkey2SAP[1][section_num]], mc2_sectionkey[mc2_sectionkey2SAP[0][tspans_outsection[section_num][nth_outsection]]], mc2_sectionkey[mc2_sectionkey2SAP[1][tspans_outsection[section_num][nth_outsection]]]);
			counter++;
		}
	}
	exit(1);
*/

	





	return;
}




/***** start of function mc2_conv_to_array *************************************/
/* called only by main */


void mc2_conv_endhinges_to_array(void)
{



	int i, j;

	unsigned long int section_num;
	unsigned long int sec1key;
	struct endhinge *hinge_to_free;

	for(i=1; i<=max_keynum; i++){
		mc2_num_left_endhinges[i]=0;
		mc2_num_right_endhinges[i]=0;
	}

	for (section_num = 1; section_num <= max_sections; section_num++) {
		mc2_currentendhinge[section_num] = mc2_firstendhinge[section_num];
		while ((*mc2_currentendhinge[section_num]).nextendhinge != NULL) {
			mc2_currentendhinge[section_num] = (*mc2_currentendhinge[section_num]).nextendhinge;

			//get keynum
			sec1key=0;
			for(i=1; i<=max_keynum; i++){
//				printf("keynum=%d, match: %lu,%lu ?\n", i, mc2_sectionkey2SAP[i][0], mc2_sectionkey2SAP[i][1]);
				if(mc2_sectionkey2SAP[i][0]==section_num && mc2_sectionkey2SAP[i][1]==(*mc2_currentendhinge[section_num]).sec2){
//					printf("sec1key=%d\n", i);
					sec1key=i;
					break;
				}
			}
			if(sec1key==0){
				printf("error: didn't find sec1key. Exitting\n");
				exit(1);
			}

			if(mc2_currentendhinge[section_num]->side==0){
				mc2_num_left_endhinges[sec1key]++;
			}
			else{
				mc2_num_right_endhinges[sec1key]++;
			}
		}
	}

	for(section_num=1; section_num<=max_keynum; section_num++){
//		printf("mc2_num_left_endhinges[%lu]=%lu\n", section_num, mc2_num_left_endhinges[section_num]);
//		printf("mc2_num_right_endhinges[%lu]=%lu\n", section_num, mc2_num_right_endhinges[section_num]);
		mc2_Lend_num_walks[section_num] = mc2_unsgn_vecalloc(1, mc2_num_left_endhinges[section_num]);
		mc2_Rend_num_walks[section_num] = mc2_unsgn_vecalloc(1, mc2_num_right_endhinges[section_num]);
		mc2_Lend_num_walks2[section_num] = mc2_unsgn_vecalloc(1, mc2_num_left_endhinges[section_num]);
		mc2_Rend_num_walks2[section_num] = mc2_unsgn_vecalloc(1, mc2_num_right_endhinges[section_num]);
	}





	//allocate space for mc2_Lend_start, mc2_Lend_end, and mc2_Lend_walks AND mc2_Lend_start2, mc2_Lend_end2, and mc2_Lend_walks2 (and Rend stuff)

	for(section_num=1; section_num<=max_keynum; section_num++){
		mc2_Lend_start[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_left_endhinges[section_num], 3, MAX_vMvL + 1, "2SAP left end starts");
		mc2_Lend_start2[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_left_endhinges[section_num], 3, MAX_vMvL + 1, "2SAP left end starts2");
		mc2_Lend_end[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_left_endhinges[section_num], 3, MAX_vMvL + 1, "2SAP left end ends");
		mc2_Lend_end2[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_left_endhinges[section_num], 3, MAX_vMvL + 1, "2SAP left end ends2");
		mc2_Lend_walks[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_left_endhinges[section_num], vM * vL + 1, MAX_vMvL + 1, "2SAP left end walks");
		mc2_Lend_walks2[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_left_endhinges[section_num], vM * vL + 1, MAX_vMvL + 1, "2SAP left end walks2");
		mc2_Rend_start[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_right_endhinges[section_num], 3, MAX_vMvL + 1, "2SAP right end starts");
		mc2_Rend_start2[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_right_endhinges[section_num], 3, MAX_vMvL + 1, "2SAP right end starts2");
		mc2_Rend_end[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_right_endhinges[section_num], 3, MAX_vMvL + 1, "2SAP right end ends");
		mc2_Rend_end2[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_right_endhinges[section_num], 3, MAX_vMvL + 1, "2SAP right end ends2");
		mc2_Rend_walks[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_right_endhinges[section_num], vM * vL + 1, MAX_vMvL + 1, "2SAP right end walks");
		mc2_Rend_walks2[section_num] = mc_2sap_alloc_int3_fixed(mc2_num_right_endhinges[section_num], vM * vL + 1, MAX_vMvL + 1, "2SAP right end walks2");
	}




	///Fill L/mc2_Rend_start, L/mc2_Rend_end, L/mc2_Rend_walks
	unsigned int L_curArrayEnt[MAX_KEYNUM_ARR+1];
	unsigned int R_curArrayEnt[MAX_KEYNUM_ARR+1];
	for(i=1; i<=max_keynum; i++){
		L_curArrayEnt[i]=1;
		R_curArrayEnt[i]=1;
	}
	for (section_num = 1; section_num <= max_sections; section_num++) {
		mc2_currentendhinge[section_num] = mc2_firstendhinge[section_num];
		while ((*mc2_currentendhinge[section_num]).nextendhinge != NULL) {
			hinge_to_free = mc2_currentendhinge[section_num];
			mc2_currentendhinge[section_num] = (*mc2_currentendhinge[section_num]).nextendhinge;
			free(hinge_to_free);

			//get keynum
			sec1key=0;
			for(i=1; i<=max_keynum; i++){
//				printf("keynum=%d, match: %lu,%lu ?\n", i, mc2_sectionkey2SAP[i][0], mc2_sectionkey2SAP[i][1]);
				if(mc2_sectionkey2SAP[i][0]==section_num && mc2_sectionkey2SAP[i][1]==(*mc2_currentendhinge[section_num]).sec2){
//					printf("sec1key=%d\n", i);
					sec1key=i;
					break;
				}
			}
			if(sec1key==0){
				printf("error: didn't find sec1key. Exitting\n");
				exit(1);
			}

			if(mc2_currentendhinge[section_num]->side==0){	//leftendhinge
				mc2_Lend_num_walks[sec1key][L_curArrayEnt[sec1key]] = mc2_currentendhinge[section_num]->num_of_walks;
				mc2_Lend_num_walks2[sec1key][L_curArrayEnt[sec1key]] = mc2_currentendhinge[section_num]->num_of_walks2;
				for(i=0; i<=(mc2_currentendhinge[section_num]->num_of_walks)-1; i++){
					for(j=0; j<=2; j++){
						mc2_Lend_start[sec1key][L_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->start[j][i];
						mc2_Lend_end[sec1key][L_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->end[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						mc2_Lend_walks[sec1key][L_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->walks[j][i];
					}
				}
				for(i=0; i<=(mc2_currentendhinge[section_num]->num_of_walks2)-1; i++){
					for(j=0; j<=2; j++){
						mc2_Lend_start2[sec1key][L_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->start2[j][i];
						mc2_Lend_end2[sec1key][L_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->end2[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						mc2_Lend_walks2[sec1key][L_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->walks2[j][i];
					}
				}
				L_curArrayEnt[sec1key]++;
			}
			else{ //rightendhinge
				mc2_Rend_num_walks[sec1key][R_curArrayEnt[sec1key]] = mc2_currentendhinge[section_num]->num_of_walks;
				mc2_Rend_num_walks2[sec1key][R_curArrayEnt[sec1key]] = mc2_currentendhinge[section_num]->num_of_walks2;
				for(i=0; i<=(mc2_currentendhinge[section_num]->num_of_walks)-1; i++){
					for(j=0; j<=2; j++){
						mc2_Rend_start[sec1key][R_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->start[j][i];
						mc2_Rend_end[sec1key][R_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->end[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						mc2_Rend_walks[sec1key][R_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->walks[j][i];
					}
				}
				for(i=0; i<=(mc2_currentendhinge[section_num]->num_of_walks2)-1; i++){
					for(j=0; j<=2; j++){
						mc2_Rend_start2[sec1key][R_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->start2[j][i];
						mc2_Rend_end2[sec1key][R_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->end2[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						mc2_Rend_walks2[sec1key][R_curArrayEnt[sec1key]][j][i] = mc2_currentendhinge[section_num]->walks2[j][i];
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
mc2_enterhinge(int i, int j, int side, int (*pointordNum)[6], int curlength)
{
	if( !(mc2_alreadyentered[i][j]==1 && side==0) ){
	//	printf("entered hinge (i=%d, j=%d, side=%d\n", i, j, side);

		curlength++;
		mc2_num_walks++;
		if(side==0){
			mc2_curstart[0][mc2_num_walks-1]=0;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 1;
		}
		else{ //side==1
			mc2_curstart[0][mc2_num_walks-1]=2;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 2;
		}
		mc2_curstart[1][mc2_num_walks-1]=j;
		mc2_curstart[2][mc2_num_walks-1]=i;


		mc2_hingestatus[i][j]=1;
		mc2_ordertemplate[side][i][j] = (*pointordNum)[side];
			/* section edge through which the SAW enters the hinge is given a number */
			/* this number represents the order in which it occurs in the section    */

		(*pointordNum)[side]++;
			/* the edge number through which the hinge is entered */
			/* is incremented so that the next edge gets a higher number */

		/* the SAW can now do one of three things */
		/* it may leave the hinge in which case mc2_leavehinge(...) is called */
		/* it may explore edges in the horizontal direction => mc2_rowedges (...) is called */
		/* it may explore edges in the vertical direction => coleges (...) is called    */
		/* in each case the currrent state of the hinge is passed via the parameters    */

		mc2_leavehinge(i, j, side ^ 1, pointordNum, curlength);
			/* the SAW leaves the hinge out the opposite side, hence (side ^ 1) is passed */
			/* if side = 1 then side ^ 1 = 0 */
			/* if side = 0 then side ^ 1 = 1 */
		mc2_rowedges(i, j, pointordNum, curlength);
			/* horizontal edges are explored */
		mc2_coledges(i, j, pointordNum, curlength);
			/* vertical edges are explored */

		(*pointordNum)[side]--;
			/* the edge number is decremented since, control is being passed back to  */
			/* the function that called mc2_enterhinge. That is to say, the program will  */
			/* now explore alternate ways of entering the hinge (i.e. from different  */
			/* locations). In essense the intial incremention above is nullified so	  */
			/* when mc2_enterhinge is called again the entering edge is numbered properly */
		
		mc2_ordertemplate[side][i][j] = 0;
			/* the edge that was recorded must now be erased as this edge may not exist */
			/* as part of the SAW when alternate possibilities are explored. In the case */
			/* that it does not exist it needs to be numbered as "0" (zero).             */
		mc2_hingestatus[i][j]=0;

		mc2_curstart[0][mc2_num_walks-1]=-1;
		mc2_curstart[1][mc2_num_walks-1]=-1;
		mc2_curstart[2][mc2_num_walks-1]=-1;
		mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;

		curlength--;
		mc2_num_walks--;
	//	printf("finished mc2_enterhinge\n");
	}
	return;

}				/* end of function mc2_enterhinge */

/***************************************************************************/

void
mc2_leavehinge(int i, int j, int side, int (*pointordNum)[6], int curlength)
/* called by: mc2_rowedges,mc2_coledges,mc2_enterhinge */
{
	if( !(mc2_alreadyentered[i][j]==1 && side==0) ){
	//	printf("left hinge (i=%d, j=%d, side=%d\n", i, j, side);

		curlength++;
		if(side==0){
			mc2_curend[0][mc2_num_walks-1]=0;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 2;
		}
		else{ //side==1
			mc2_curend[0][mc2_num_walks-1]=2;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 1;
		}
		mc2_curend[1][mc2_num_walks-1]=j;
		mc2_curend[2][mc2_num_walks-1]=i;


		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		mc2_ordertemplate[side][i][j] = (*pointordNum)[side];	/* record the section
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
						if (mc2_hingestatus[ii][jj] == 0) {
							// if the vertex is not occupied then explore entering the hinge
							// at this vertex with a second polygon
							mc2_enterhinge2(ii, jj, side, pointordNum, 0);
						}
						mc2_alreadyentered2[ii][jj]=1;
					}
				}
				for (ii = 0; ii <= M; ii++) {
					for (jj = 0; jj <= L; jj++) {
						mc2_alreadyentered2[ii][jj]=0;
					}
				}
			}
		}
		/* since a call to mc2_leavehinge means the "headernode" is no longer in the two span     */
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
				if (mc2_hingestatus[ii][jj] == 0) {
				/* if the vertex is not occupied then explore entering the hinge at this vertex */
				/* this maintains the self avoiding nature of the walk */
					mc2_enterhinge(ii, jj, side, pointordNum,0);
				}
			}
		}
		(*pointordNum)[side]--;
			/* the edge number is decremented since, control is being passed back to  */
			/* the function that called mc2_leavehinge. That is to say, the program will  */
			/* now explore alternate ways of leaving the hinge (i.e. via different    */
			/* paths. this may be as simple as going out the other side). In essense  */
			/* the intial incremention above is nullified so when mc2_leavehinge is 	  */
			/* called again the leaving edge is numbered properly  			  */
		
		mc2_ordertemplate[side][i][j] = 0;
			/* the edge that was recorded must now be erased as this edge may not exist  */
			/* as part of the SAW when alternate possibilities are explored. In the case */
			/* that it does not exist it needs to be numbered as "0" (zero).             */
	//		printf("finished mc2_leavehinge\n");


		mc2_curend[0][mc2_num_walks-1]=-1;
		mc2_curend[1][mc2_num_walks-1]=-1;
		mc2_curend[2][mc2_num_walks-1]=-1;
		mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;
		curlength--;
	}
	return;
}
/***************************************************************************/

void mc2_rowedges(int i, int j, int (*pointordNum)[6], int curlength)  //y-direction
{
//	printf("mc2_rowedges called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(mc2_hingestatus[i][j-1]==0){
			mc2_hingestatus[i][j-1]=1;

			curlength++;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 4;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			mc2_rowhingeedges[i][j-1]=1;
			mc2_leavehinge(i, j-1, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			mc2_leavehinge(i, j-1, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			mc2_rowedges(i, j-1, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			mc2_coledges(i, j-1, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			mc2_hingestatus[i][j-1]=0;
			mc2_rowhingeedges[i][j-1]=0;

			mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;
			curlength--;
		}
	}
	if(j<L){
		if(mc2_hingestatus[i][j+1]==0){
			mc2_hingestatus[i][j+1]=1;

			curlength++;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 3;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			mc2_rowhingeedges[i][j]=1;
			mc2_leavehinge(i, j+1, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			mc2_leavehinge(i, j+1, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			mc2_rowedges(i, j+1, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			mc2_coledges(i, j+1, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			mc2_hingestatus[i][j+1]=0;
			mc2_rowhingeedges[i][j]=0;

			mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;
			curlength--;
		}
	}
//	printf("finished mc2_rowedges\n");
	return;
}				/* end of function mc2_rowedges */

/***************************************************************************/

void mc2_coledges(int i, int j, int (*pointordNum)[6], int curlength)  //z-direction
/* called by mc2_rowedges,mc2_coledges */
{
//	printf("mc2_coledges called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(mc2_hingestatus[i-1][j]==0){
			mc2_hingestatus[i-1][j]=1;

			curlength++;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 6;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			mc2_colhingeedges[i-1][j]=1;
			mc2_leavehinge(i-1, j, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			mc2_leavehinge(i-1, j, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			mc2_rowedges(i-1, j, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			mc2_coledges(i-1, j, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			mc2_hingestatus[i-1][j]=0;
			mc2_colhingeedges[i-1][j]=0;

			mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;
			curlength--;
		}
	}
	if(i<M){
		if(mc2_hingestatus[i+1][j]==0){
			mc2_hingestatus[i+1][j]=1;

			curlength++;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 5;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			mc2_colhingeedges[i][j]=1;
			mc2_leavehinge(i+1, j, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			mc2_leavehinge(i+1, j, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			mc2_rowedges(i+1, j, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			mc2_coledges(i+1, j, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			mc2_hingestatus[i+1][j]=0;
			mc2_colhingeedges[i][j]=0;

			mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;
			curlength--;
		}
	}
//	printf("finished mc2_coledges\n");
	return;
}				/* end of function mc2_coledges */

/***************************************************************************/







/***************************************************************************/

void
mc2_enterhinge2(int i, int j, int side, int (*pointordNum)[6], int curlength2)
{
	if( !(mc2_alreadyentered2[i][j]==1 && side==0) ){
	//	printf("entered hinge2 (i=%d, j=%d, side=%d\n", i, j, side);

		curlength2++;
		mc2_num_walks2++;
		if(side==0){
			mc2_curstart2[0][mc2_num_walks2-1]=0;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 1;
		}
		else{ //side==1
			mc2_curstart2[0][mc2_num_walks2-1]=2;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 2;
		}
		mc2_curstart2[1][mc2_num_walks2-1]=j;
		mc2_curstart2[2][mc2_num_walks2-1]=i;


		mc2_hingestatus[i][j]=1;
		mc2_ordertemplate2[side][i][j] = (*pointordNum)[side+3];
			/* section edge through which the SAW enters the hinge is given a number */
			/* this number represents the order in which it occurs in the section    */

		(*pointordNum)[side+3]++;
			/* the edge number through which the hinge is entered */
			/* is incremented so that the next edge gets a higher number */

		/* the SAW can now do one of three things */
		/* it may leave the hinge in which case mc2_leavehinge(...) is called */
		/* it may explore edges in the horizontal direction => mc2_rowedges (...) is called */
		/* it may explore edges in the vertical direction => coleges (...) is called    */
		/* in each case the currrent state of the hinge is passed via the parameters    */

		mc2_leavehinge2(i, j, side ^ 1, pointordNum, curlength2);
			/* the SAW leaves the hinge out the opposite side, hence (side ^ 1) is passed */
			/* if side = 1 then side ^ 1 = 0 */
			/* if side = 0 then side ^ 1 = 1 */
		mc2_rowedges2(i, j, pointordNum, curlength2);
			/* horizontal edges are explored */
		mc2_coledges2(i, j, pointordNum, curlength2);
			/* vertical edges are explored */

		(*pointordNum)[side+3]--;
			/* the edge number is decremented since, control is being passed back to  */
			/* the function that called mc2_enterhinge. That is to say, the program will  */
			/* now explore alternate ways of entering the hinge (i.e. from different  */
			/* locations). In essense the intial incremention above is nullified so	  */
			/* when mc2_enterhinge is called again the entering edge is numbered properly */
		
		mc2_ordertemplate2[side][i][j] = 0;
			/* the edge that was recorded must now be erased as this edge may not exist */
			/* as part of the SAW when alternate possibilities are explored. In the case */
			/* that it does not exist it needs to be numbered as "0" (zero).             */
		mc2_hingestatus[i][j]=0;

		mc2_curstart2[0][mc2_num_walks2-1]=-1;
		mc2_curstart2[1][mc2_num_walks2-1]=-1;
		mc2_curstart2[2][mc2_num_walks2-1]=-1;
		mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;

		curlength2--;
		mc2_num_walks2--;
	//	printf("finished mc2_enterhinge\n");
	}
	return;

}				/* end of function mc2_enterhinge */

/***************************************************************************/

void
mc2_leavehinge2(int i, int j, int side, int (*pointordNum)[6], int curlength2)
/* called by: mc2_rowedges,mc2_coledges,mc2_enterhinge */
{
	if( !(mc2_alreadyentered2[i][j]==1 && side==0) ){
	//	printf("left hinge (i=%d, j=%d, side=%d\n", i, j, side);

		curlength2++;
		if(side==0){
			mc2_curend2[0][mc2_num_walks2-1]=0;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 2;
		}
		else{ //side==1
			mc2_curend2[0][mc2_num_walks2-1]=2;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 1;
		}
		mc2_curend2[1][mc2_num_walks2-1]=j;
		mc2_curend2[2][mc2_num_walks2-1]=i;


		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		mc2_ordertemplate2[side][i][j] = (*pointordNum)[side+3];	/* record the section
									 * edge (i.e. the
									 * leaving edge) */
		(*pointordNum)[side+3]++;
			/* the edge number is incremented so that the next edge gets a higher number */

		if (side == 0 && (*pointordNum)[4] > 1) {
	 	/* if the leaving edge is on the left and there is at least one edge on the right section */
		/* this may be a valid two span */
			mc2_num_2_spans++; /*Keep track of all 2 spans, valid and non-valid, just for interest sake*/
			if (mc2_LFlag2(pointordNum)){
			/* if the two span connects to phi on the left and on the right then do the following*/

				mc2_valid_2_spans++;                //This is a valid 2 span so it is counted
				//printf("Pattern number: %d \n", mc2_valid_2_spans);
				mc2_fillreordertemplate((*pointordNum)[0]-1, (*pointordNum)[1]-1);
				mc2_fillreordertemplate2((*pointordNum)[3]-1, (*pointordNum)[4]-1);
				mc2_recordtemplate(pointordNum);   //the information contained in mc2_ordertemplate is recorded
			}

		}
		/* since a call to mc2_leavehinge means the "headernode" is no longer in the two span     */
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
				if (mc2_hingestatus[ii][jj] == 0) {
				/* if the vertex is not occupied then explore entering the hinge at this vertex */
				/* this maintains the self avoiding nature of the walk */
					mc2_enterhinge2(ii, jj, side, pointordNum,0);
				}
			}
		}
		(*pointordNum)[side+3]--;
			/* the edge number is decremented since, control is being passed back to  */
			/* the function that called mc2_leavehinge. That is to say, the program will  */
			/* now explore alternate ways of leaving the hinge (i.e. via different    */
			/* paths. this may be as simple as going out the other side). In essense  */
			/* the intial incremention above is nullified so when mc2_leavehinge is 	  */
			/* called again the leaving edge is numbered properly  			  */
		
		mc2_ordertemplate2[side][i][j] = 0;
			/* the edge that was recorded must now be erased as this edge may not exist  */
			/* as part of the SAW when alternate possibilities are explored. In the case */
			/* that it does not exist it needs to be numbered as "0" (zero).             */
	//		printf("finished mc2_leavehinge\n");


		mc2_curend2[0][mc2_num_walks2-1]=-1;
		mc2_curend2[1][mc2_num_walks2-1]=-1;
		mc2_curend2[2][mc2_num_walks2-1]=-1;
		mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;
		curlength2--;
	}
	return;
}
/***************************************************************************/

void mc2_rowedges2(int i, int j, int (*pointordNum)[6], int curlength2)  //y-direction
{
//	printf("mc2_rowedges2 called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(mc2_hingestatus[i][j-1]==0){
			mc2_hingestatus[i][j-1]=1;

			curlength2++;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 4;

			(*pointordNum)[5]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			mc2_rowhingeedges[i][j-1]=1;
			mc2_leavehinge2(i, j-1, 0, pointordNum, curlength2);
				/* exit the hinge on side 0 with j moved*/
			mc2_leavehinge2(i, j-1, 1, pointordNum, curlength2);
				/* exit the hinge on side 1 with j moved*/
			mc2_rowedges2(i, j-1, pointordNum, curlength2);
				/* explore horizontal moves within the hinge */
			mc2_coledges2(i, j-1, pointordNum, curlength2);
				/* explore vertical moves within the hinge */
			(*pointordNum)[5]--;
				/* decrement the number of edges in the hinge */
			mc2_hingestatus[i][j-1]=0;
			mc2_rowhingeedges[i][j-1]=0;

			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;
			curlength2--;
		}
	}
	if(j<L){
		if(mc2_hingestatus[i][j+1]==0){
			mc2_hingestatus[i][j+1]=1;

			curlength2++;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 3;

			(*pointordNum)[5]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			mc2_rowhingeedges[i][j]=1;
			mc2_leavehinge2(i, j+1, 0, pointordNum, curlength2);
				/* exit the hinge on side 0 with j moved*/
			mc2_leavehinge2(i, j+1, 1, pointordNum, curlength2);
				/* exit the hinge on side 1 with j moved*/
			mc2_rowedges2(i, j+1, pointordNum, curlength2);
				/* explore horizontal moves within the hinge */
			mc2_coledges2(i, j+1, pointordNum, curlength2);
				/* explore vertical moves within the hinge */
			(*pointordNum)[5]--;
				/* decrement the number of edges in the hinge */
			mc2_hingestatus[i][j+1]=0;
			mc2_rowhingeedges[i][j]=0;

			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;
			curlength2--;
		}
	}
//	printf("finished mc2_rowedges2\n");
	return;
}				/* end of function mc2_rowedges */

/***************************************************************************/

void mc2_coledges2(int i, int j, int (*pointordNum)[6], int curlength2)  //z-direction
/* called by mc2_rowedges,mc2_coledges */
{
//	printf("mc2_coledges2 called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(mc2_hingestatus[i-1][j]==0){
			mc2_hingestatus[i-1][j]=1;

			curlength2++;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 6;

			(*pointordNum)[5]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			mc2_colhingeedges[i-1][j]=1;
			mc2_leavehinge2(i-1, j, 0, pointordNum, curlength2);
				/* exit the hinge on side 0 with j moved*/
			mc2_leavehinge2(i-1, j, 1, pointordNum, curlength2);
				/* exit the hinge on side 1 with j moved*/
			mc2_rowedges2(i-1, j, pointordNum, curlength2);
				/* explore horizontal moves within the hinge */
			mc2_coledges2(i-1, j, pointordNum, curlength2);
				/* explore vertical moves within the hinge */
			(*pointordNum)[5]--;
				/* decrement the number of edges in the hinge */
			mc2_hingestatus[i-1][j]=0;
			mc2_colhingeedges[i-1][j]=0;

			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;
			curlength2--;
		}
	}
	if(i<M){
		if(mc2_hingestatus[i+1][j]==0){
			mc2_hingestatus[i+1][j]=1;

			curlength2++;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 5;

			(*pointordNum)[5]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			mc2_colhingeedges[i][j]=1;
			mc2_leavehinge2(i+1, j, 0, pointordNum, curlength2);
				/* exit the hinge on side 0 with j moved*/
			mc2_leavehinge2(i+1, j, 1, pointordNum, curlength2);
				/* exit the hinge on side 1 with j moved*/
			mc2_rowedges2(i+1, j, pointordNum, curlength2);
				/* explore horizontal moves within the hinge */
			mc2_coledges2(i+1, j, pointordNum, curlength2);
				/* explore vertical moves within the hinge */
			(*pointordNum)[5]--;
				/* decrement the number of edges in the hinge */
			mc2_hingestatus[i+1][j]=0;
			mc2_colhingeedges[i][j]=0;

			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;
			curlength2--;
		}
	}
//	printf("finished mc2_coledges2\n");
	return;
}				/* end of function mc2_coledges */

/***************************************************************************/



















/***************************************************************************/
void
mc2_recordtemplate( int (*pointordNum)[6])
/* This function records the pertinent information about a two-span */
/* which has been verified as being valid into a linked list */
{

	(void)pointordNum;

	int i,j;

	unsigned long int	inNum = mc2_num_section_norder(0);
//	printf("OGinNum=%d\n", inNum);
	//Change inNum to it's key number.
	i=1;
	int found=0;
	while(mc2_sectionkey[i]!=0){
		if(mc2_sectionkey[i]==inNum){
			found=1;
			inNum=i;
			break;
		}
		i++;
	}
	if(found==0){	//need to record inNum in mc2_sectionkey
		mc2_sectionkey[i]=inNum;
		inNum=i;
	}
//	printf("NEWinNum=%d\n", inNum);

	unsigned long int	outNum = mc2_num_section_norder(1);
//	printf("OGoutNum=%d\n", outNum);
	//reset outNum to it's key number.
	i=1;
	found=0;
	while(mc2_sectionkey[i]!=0){
		if(mc2_sectionkey[i]==outNum){
			found=1;
			outNum=i;
			break;
		}
		i++;
	}
	if(found==0){	//need to record inNum in mc2_sectionkey
		mc2_sectionkey[i]=outNum;
		outNum=i;
	}
//	printf("NEWoutNum=%d\n", outNum);


	unsigned long int	inNum2 = mc2_num_section_norder2(0);
//	printf("OGinNum2=%d\n", inNum2);
	//Change inNum2 to it's key number.
	i=1;
	found=0;
	while(mc2_sectionkey[i]!=0){
		if(mc2_sectionkey[i]==inNum2){
			found=1;
			inNum2=i;
			break;
		}
		i++;
	}
	if(found==0){	//need to record inNum in mc2_sectionkey
		mc2_sectionkey[i]=inNum2;
		inNum2=i;
	}
//	printf("NEWinNum2=%d\n", inNum);

	unsigned long int	outNum2 = mc2_num_section_norder2(1);
//	printf("OGoutNum2=%d\n", outNum2);
	//reset outNum2 to it's key number.
	i=1;
	found=0;
	while(mc2_sectionkey[i]!=0){
		if(mc2_sectionkey[i]==outNum2){
			found=1;
			outNum2=i;
			break;
		}
		i++;
	}
	if(found==0){	//need to record inNum in mc2_sectionkey
		mc2_sectionkey[i]=outNum2;
		outNum2=i;
	}
//	printf("NEWoutNum2=%d\n", outNum2);


	int temp_hedges[MAX_vMvL];
	int index=0;
	int edgenum=1;


	for(i=0; i<=vM*vL-1; i++){
		temp_hedges[i]=0;
	}

	for(i=0; i<=vM-1; i++){
		for(j=0; j<=L-1; j++){
			if(mc2_rowhingeedges[i][j]==1){
				temp_hedges[index]=edgenum;
				index++;
			}
			edgenum++;
		}
	}
	for(i=0; i<=M-1; i++){
		for(j=0; j<=vL-1; j++){
			if(mc2_colhingeedges[i][j]==1){
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

	struct hinge_span* test = mc2_first_hinge_span[inNum];
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
					mc2_dupcounter++;
					break;
				}
			}
		}
	}

//	printf("duplicate=%d\n", duplicate);
	if(duplicate==0){
		mc2_num_tspans++;

//		printf("actually recording\n");
		(*mc2_current_hinge_span[inNum]).nexthinge = mc2_newhinge(mc2_num_walks, mc2_num_walks2); /* add a mc2_newhinge to the linked list */
	
		mc2_current_hinge_span[inNum] = (*mc2_current_hinge_span[inNum]).nexthinge; /* let mc2_current_hinge_span[inNum] point to mc2_newhinge */
		(*mc2_current_hinge_span[inNum]).inorder2 = inNum2;	/* record the ordering on the out side */
		(*mc2_current_hinge_span[inNum]).outorder = outNum;	/* record the ordering on the out side */
		(*mc2_current_hinge_span[inNum]).outorder2 = outNum2;	/* record the ordering on the out side */
		(*mc2_current_hinge_span[inNum]).num_of_walks = mc2_num_walks;
		(*mc2_current_hinge_span[inNum]).num_of_walks2 = mc2_num_walks2;
		for(i=0; i<=vM*vL-1; i++){
			(*mc2_current_hinge_span[inNum]).hedges[i]=temp_hedges[i];
		}

		for(i=0; i<=mc2_num_walks-1; i++){
			for(j=0; j<=2; j++){
				mc2_current_hinge_span[inNum]->start[j][i] = mc2_curstart[j][i];
				mc2_current_hinge_span[inNum]->end[j][i] = mc2_curend[j][i];
			}
			for(j=0; j<=vM*vL; j++){
				mc2_current_hinge_span[inNum]->walks[j][i] = mc2_curwalks[j][i];
			}
		}

		for(i=0; i<=mc2_num_walks2-1; i++){
			for(j=0; j<=2; j++){
				mc2_current_hinge_span[inNum]->start2[j][i] = mc2_curstart2[j][i];
				mc2_current_hinge_span[inNum]->end2[j][i] = mc2_curend2[j][i];
			}
			for(j=0; j<=vM*vL; j++){
				mc2_current_hinge_span[inNum]->walks2[j][i] = mc2_curwalks2[j][i];
			}
		}

/*		printf("Recorded 2-span. Sections %d,%d to Sections %d,%d with hedges: ", inNum, inNum2, outNum, outNum2);
		for(i=0; i<=vM*vL-1; i++){
			printf("%d ", (*mc2_current_hinge_span[inNum]).hedges[i]);
		}
		printf("\n");
*/	}
	else{
		mc2_num_duplicate_tspans++;
	}
//	printf("finished recording\n");
	return;
}




///////////////ENDHINGE FUNCTIONS///////////////////////////////////////////////////////////////////


void
mc2_enterendhinge(int i, int j, int curlength)
{
	if( mc2_alreadyentered[i][j]!=1 ){
	//	printf("entered endhinge (i=%d, j=%d\n", i, j);

		curlength++;
		mc2_num_walks++;
		mc2_curstart[0][mc2_num_walks-1]=0;
		mc2_curstart[1][mc2_num_walks-1]=j;
		mc2_curstart[2][mc2_num_walks-1]=i;

		mc2_curwalks[curlength-1][mc2_num_walks-1] = 1;

		mc2_hingestatus[i][j]=1;
		mc2_endtemplate[i][j] = mc2_EndOrdNum[0];
		mc2_EndOrdNum[0]++;


		mc2_endrowedges(i, j, curlength);
			/* horizontal edges are explored */
		mc2_endcoledges(i, j, curlength);
			/* vertical edges are explored */

		//reset stuff
		mc2_EndOrdNum[0]--;
		mc2_endtemplate[i][j] = 0;
		mc2_hingestatus[i][j]=0;

		mc2_curstart[0][mc2_num_walks-1]=-1;
		mc2_curstart[1][mc2_num_walks-1]=-1;
		mc2_curstart[2][mc2_num_walks-1]=-1;
		mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;

		curlength--;
		mc2_num_walks--;
	//	printf("finished mc2_enterhinge\n");
	}
	return;

}				/* end of function mc2_enterhinge */

/***************************************************************************/

void
mc2_leaveendhinge(int i, int j, int curlength)
/* called by: mc2_rowedges,mc2_coledges,mc2_enterhinge */
{
	if( mc2_alreadyentered[i][j]!=1 ){
	//	printf("left hinge (i=%d, j=%d\n", i, j);

		curlength++;
		mc2_curend[0][mc2_num_walks-1]=0;
		mc2_curend[1][mc2_num_walks-1]=j;
		mc2_curend[2][mc2_num_walks-1]=i;
		mc2_curwalks[curlength-1][mc2_num_walks-1] = 2;

		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		mc2_endtemplate[i][j] = mc2_EndOrdNum[0];
		mc2_EndOrdNum[0]++;

		if ( mc2_EndOrdNum[0]-1 <= vM*vL-2 ){	//make sure there's room for the 2nd poly to enter+exit
			for (ii = 0; ii <= M; ii++) { // these are the "vertical" locations of the vertices being considered
				for (jj = 0; jj <= L; jj++) { // these are the horizontal locations of the vertices being considered
					if (mc2_hingestatus[ii][jj] == 0) {
						// if the vertex is not occupied then explore entering the hinge
						// at this vertex with a second polygon
						mc2_enterendhinge2(ii, jj, 0);
					}
					mc2_alreadyentered2[ii][jj]=1;
				}
			}
			for (ii = 0; ii <= M; ii++) {
				for (jj = 0; jj <= L; jj++) {
					mc2_alreadyentered2[ii][jj]=0;
				}
			}
		}




		for (ii = 0; ii <= M; ii++) {
		/* these are the "vertical" locations of the vertices being considered */
			for (jj = 0; jj <= L; jj++) {
			/* these are the horizontal locations of the vertices being considered */
				if (mc2_hingestatus[ii][jj] == 0) {
				/* if the vertex is not occupied then explore entering the hinge at this vertex */
				/* this maintains the self avoiding nature of the walk */
					mc2_enterendhinge(ii, jj, 0);
				}
			}
		}

		mc2_EndOrdNum[0]--;
		mc2_endtemplate[i][j] = 0;

		mc2_curend[0][mc2_num_walks-1]=-1;
		mc2_curend[1][mc2_num_walks-1]=-1;
		mc2_curend[2][mc2_num_walks-1]=-1;
		mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;
		curlength--;
	}
	return;
}
/***************************************************************************/

void mc2_endrowedges(int i, int j, int curlength)  //y-direction
{
//	printf("mc2_endrowedges called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(mc2_hingestatus[i][j-1]==0){
			mc2_hingestatus[i][j-1]=1;

			curlength++;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 4;

			mc2_EndOrdNum[1]++;
			mc2_rowhingeedges[i][j-1]=1;

			mc2_leaveendhinge(i, j-1, curlength);
			mc2_endrowedges(i, j-1, curlength);
			mc2_endcoledges(i, j-1, curlength);

			mc2_rowhingeedges[i][j-1]=0;
			mc2_EndOrdNum[1]--;
			mc2_hingestatus[i][j-1]=0;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;
			curlength--;
		}
	}
	if(j<L){
		if(mc2_hingestatus[i][j+1]==0){
			mc2_hingestatus[i][j+1]=1;

			curlength++;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 3;

			mc2_EndOrdNum[1]++;
			mc2_rowhingeedges[i][j]=1;

			mc2_leaveendhinge(i, j+1, curlength);
			mc2_endrowedges(i, j+1, curlength);
			mc2_endcoledges(i, j+1, curlength);

			mc2_rowhingeedges[i][j]=0;
			mc2_EndOrdNum[1]--;
			mc2_hingestatus[i][j+1]=0;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;
			curlength--;
		}
	}
//	printf("finished mc2_endrowedges\n");
	return;
}				/* end of function mc2_rowedges */

/***************************************************************************/

void mc2_endcoledges(int i, int j, int curlength)  //z-direction
{
//	printf("mc2_endcoledges called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(mc2_hingestatus[i-1][j]==0){
			mc2_hingestatus[i-1][j]=1;

			curlength++;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 6;

			mc2_EndOrdNum[1]++;
			mc2_colhingeedges[i-1][j]=1;

			mc2_leaveendhinge(i-1, j, curlength);
			mc2_endrowedges(i-1, j, curlength);
			mc2_endcoledges(i-1, j, curlength);

			mc2_colhingeedges[i-1][j]=0;
			mc2_EndOrdNum[1]--;
			mc2_hingestatus[i-1][j]=0;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;
			curlength--;
		}
	}
	if(i<M){
		if(mc2_hingestatus[i+1][j]==0){
			mc2_hingestatus[i+1][j]=1;

			curlength++;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 5;

			mc2_EndOrdNum[1]++;
			mc2_colhingeedges[i][j]=1;

			mc2_leaveendhinge(i+1, j, curlength);
			mc2_endrowedges(i+1, j, curlength);
			mc2_endcoledges(i+1, j, curlength);

			mc2_colhingeedges[i][j]=0;
			mc2_EndOrdNum[1]--;
			mc2_hingestatus[i+1][j]=0;
			mc2_curwalks[curlength-1][mc2_num_walks-1] = 0;
			curlength--;
		}
	}
//	printf("finished mc2_endcoledges\n");
	return;
}				/* end of function mc2_coledges */

/***************************************************************************/










void
mc2_enterendhinge2(int i, int j, int curlength2)
{
	if( mc2_alreadyentered2[i][j]!=1 ){
	//	printf("entered endhinge2 (i=%d, j=%d\n", i, j);

		curlength2++;
		mc2_num_walks2++;
		mc2_curstart2[0][mc2_num_walks2-1]=0;
		mc2_curstart2[1][mc2_num_walks2-1]=j;
		mc2_curstart2[2][mc2_num_walks2-1]=i;

		mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 1;

		mc2_hingestatus[i][j]=2;
		mc2_endtemplate3[i][j] = mc2_EndOrdNum2[0];
		mc2_EndOrdNum2[0]++;


		mc2_endrowedges2(i, j, curlength2);
			/* horizontal edges are explored */
		mc2_endcoledges2(i, j, curlength2);
			/* vertical edges are explored */

		//reset stuff
		mc2_EndOrdNum2[0]--;
		mc2_endtemplate3[i][j] = 0;
		mc2_hingestatus[i][j]=0;

		mc2_curstart2[0][mc2_num_walks2-1]=-1;
		mc2_curstart2[1][mc2_num_walks2-1]=-1;
		mc2_curstart2[2][mc2_num_walks2-1]=-1;
		mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;

		curlength2--;
		mc2_num_walks2--;
	//	printf("finished mc2_enterhinge2\n");
	}
	return;

}				/* end of function mc2_enterhinge */

/***************************************************************************/

void
mc2_leaveendhinge2(int i, int j, int curlength2)
/* called by: mc2_rowedges,mc2_coledges,mc2_enterhinge */
{
	if( mc2_alreadyentered2[i][j]!=1 ){
	//	printf("left hinge (i=%d, j=%d\n", i, j);

		curlength2++;
		mc2_curend2[0][mc2_num_walks2-1]=0;
		mc2_curend2[1][mc2_num_walks2-1]=j;
		mc2_curend2[2][mc2_num_walks2-1]=i;
		mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 2;

		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		mc2_endtemplate3[i][j] = mc2_EndOrdNum2[0];
		mc2_EndOrdNum2[0]++;

		mc2_num_potential_endhinges++;
		if (mc2_LFlag_endhinge2()){
			mc2_num_valid_endhinges++;
			mc2_fillrendtemplate();
			mc2_fillrendtemplate3();
			mc2_recordendtemplate(); /*the information contained in mc2_endtemplate is recorded*/
			//duplicate check is performed in mc2_recordtemplate.
		}



		for (ii = 0; ii <= M; ii++) {
		/* these are the "vertical" locations of the vertices being considered */
			for (jj = 0; jj <= L; jj++) {
			/* these are the horizontal locations of the vertices being considered */
				if (mc2_hingestatus[ii][jj] == 0) {
				/* if the vertex is not occupied then explore entering the hinge at this vertex */
				/* this maintains the self avoiding nature of the walk */
					mc2_enterendhinge2(ii, jj, 0);
				}
			}
		}

		mc2_EndOrdNum2[0]--;
		mc2_endtemplate3[i][j] = 0;

		mc2_curend2[0][mc2_num_walks2-1]=-1;
		mc2_curend2[1][mc2_num_walks2-1]=-1;
		mc2_curend2[2][mc2_num_walks2-1]=-1;
		mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;
		curlength2--;
	}
	return;
}
/***************************************************************************/

void mc2_endrowedges2(int i, int j, int curlength2)  //y-direction
{
//	printf("mc2_endrowedges2 called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(mc2_hingestatus[i][j-1]==0){
			mc2_hingestatus[i][j-1]=2;

			curlength2++;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 4;

			mc2_EndOrdNum2[1]++;
			mc2_rowhingeedges[i][j-1]=1;

			mc2_leaveendhinge2(i, j-1, curlength2);
			mc2_endrowedges2(i, j-1, curlength2);
			mc2_endcoledges2(i, j-1, curlength2);

			mc2_rowhingeedges[i][j-1]=0;
			mc2_EndOrdNum2[1]--;
			mc2_hingestatus[i][j-1]=0;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;
			curlength2--;
		}
	}
	if(j<L){
		if(mc2_hingestatus[i][j+1]==0){
			mc2_hingestatus[i][j+1]=2;

			curlength2++;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 3;

			mc2_EndOrdNum2[1]++;
			mc2_rowhingeedges[i][j]=1;

			mc2_leaveendhinge2(i, j+1, curlength2);
			mc2_endrowedges2(i, j+1, curlength2);
			mc2_endcoledges2(i, j+1, curlength2);

			mc2_rowhingeedges[i][j]=0;
			mc2_EndOrdNum2[1]--;
			mc2_hingestatus[i][j+1]=0;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;
			curlength2--;
		}
	}
//	printf("finished mc2_endrowedges2\n");
	return;
}				/* end of function mc2_rowedges */

/***************************************************************************/

void mc2_endcoledges2(int i, int j, int curlength2)  //z-direction
{
//	printf("mc2_endcoledges2 called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(mc2_hingestatus[i-1][j]==0){
			mc2_hingestatus[i-1][j]=2;

			curlength2++;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 6;

			mc2_EndOrdNum2[1]++;
			mc2_colhingeedges[i-1][j]=1;

			mc2_leaveendhinge2(i-1, j, curlength2);
			mc2_endrowedges2(i-1, j, curlength2);
			mc2_endcoledges2(i-1, j, curlength2);

			mc2_colhingeedges[i-1][j]=0;
			mc2_EndOrdNum2[1]--;
			mc2_hingestatus[i-1][j]=0;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;
			curlength2--;
		}
	}
	if(i<M){
		if(mc2_hingestatus[i+1][j]==0){
			mc2_hingestatus[i+1][j]=2;

			curlength2++;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 5;

			mc2_EndOrdNum2[1]++;
			mc2_colhingeedges[i][j]=1;

			mc2_leaveendhinge2(i+1, j, curlength2);
			mc2_endrowedges2(i+1, j, curlength2);
			mc2_endcoledges2(i+1, j, curlength2);

			mc2_colhingeedges[i][j]=0;
			mc2_EndOrdNum2[1]--;
			mc2_hingestatus[i+1][j]=0;
			mc2_curwalks2[curlength2-1][mc2_num_walks2-1] = 0;
			curlength2--;
		}
	}
//	printf("finished mc2_endcoledges2\n");
	return;
}				/* end of function mc2_coledges */

/***************************************************************************/












/***************************************************************************/
void
mc2_recordendtemplate()
{
//	printf("recording an mc2_endtemplate\n");
	int i,j;

	unsigned long int	secnum = mc2_num_section_endhinge_norder();
//	printf("OG secnum=%d\n", secnum);
	//Change secnum to it's key number.
	i=1;
	int found=0;
	while(mc2_sectionkey[i]!=0){
		if(mc2_sectionkey[i]==secnum){
			found=1;
			secnum=i;
			break;
		}
		i++;
	}
	if(found==0){	//This shouldn't happen, since all valid section should have been generated during 2-span geneartion.
		printf("another section was generated... problem. Exitting.\n");
		exit(1);
//		mc2_sectionkey[i]=secnum;
//		secnum=i;
	}
//	printf("NEW secnum=%lu\n", secnum);


	unsigned long int	secnum3 = mc2_num_section_endhinge_norder3();
//	printf("OG secnum3=%d\n", secnum3);
	//Change secnum to it's key number.
	i=1;
	found=0;
	while(mc2_sectionkey[i]!=0){
		if(mc2_sectionkey[i]==secnum3){
			found=1;
			secnum3=i;
			break;
		}
		i++;
	}
	if(found==0){	//This shouldn't happen, since all valid section should have been generated during 2-span geneartion.
		printf("another section was generated... problem. Exitting.\n");
		exit(1);
//		mc2_sectionkey[i]=secnum3;
//		secnum3=i;
	}
//	printf("NEW secnum3=%lu\n", secnum3);




	int temp_hedges[MAX_vMvL];
	int index=0;
	int edgenum=1;

	for(i=0; i<=vM*vL-1; i++){
		temp_hedges[i]=0;
	}

	for(i=0; i<=vM-1; i++){
		for(j=0; j<=L-1; j++){
			if(mc2_rowhingeedges[i][j]==1){
				temp_hedges[index]=edgenum;
				index++;
			}
			edgenum++;
		}
	}
	for(i=0; i<=M-1; i++){
		for(j=0; j<=vL-1; j++){
			if(mc2_colhingeedges[i][j]==1){
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

	struct endhinge* test = mc2_firstendhinge[secnum];
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
		mc2_tot_right_endhinges++;
		(*mc2_currentendhinge[secnum]).nextendhinge = mc2_newendhinge(mc2_num_walks, mc2_num_walks2); /* add a mc2_newendhinge to the linked list */
	
		mc2_currentendhinge[secnum] = (*mc2_currentendhinge[secnum]).nextendhinge;
		(*mc2_currentendhinge[secnum]).sec2 = secnum3;
		(*mc2_currentendhinge[secnum]).side = 1;
		(*mc2_currentendhinge[secnum]).num_of_walks = mc2_num_walks;
		(*mc2_currentendhinge[secnum]).num_of_walks2 = mc2_num_walks2;
		for(i=0; i<=vM*vL-1; i++){
			(*mc2_currentendhinge[secnum]).hedges[i]=temp_hedges[i];
		}

		for(i=0; i<=mc2_num_walks-1; i++){
			for(j=0; j<=2; j++){
				mc2_currentendhinge[secnum]->start[j][i] = mc2_curstart[j][i];
				mc2_currentendhinge[secnum]->end[j][i] = mc2_curend[j][i];
			}
			for(j=0; j<=vM*vL; j++){
				mc2_currentendhinge[secnum]->walks[j][i] = mc2_curwalks[j][i];
			}
		}

		for(i=0; i<=mc2_num_walks2-1; i++){
			for(j=0; j<=2; j++){
				mc2_currentendhinge[secnum]->start2[j][i] = mc2_curstart2[j][i];
				mc2_currentendhinge[secnum]->end2[j][i] = mc2_curend2[j][i];
			}
			for(j=0; j<=vM*vL; j++){
				mc2_currentendhinge[secnum]->walks2[j][i] = mc2_curwalks2[j][i];
			}
		}

	

		//now record the mirror image left endhinge version
		//first check if the left version will follow a lexicographical ordering on left endhinges: make poly1 be the first in the lexicographical ordering (of vertices)
		//lexi_left_endhinges counts number of left endhinges which have poly1 as the first (lexicographically) poly.

		int recordleft=0;	//0 if haven't found an occupied vertex yet. 1 if first vertex is from poly1. 2 if first vertex if from poly2

		for(i=0; i<=M; i++){
			for(j=0; j<=L; j++){
				if(recordleft==0){
					if(mc2_hingestatus[i][j]==1){
						recordleft=1;
					}
					else if(mc2_hingestatus[i][j]==2){
						recordleft=2;
					}
				}
				else{
					break;
				}
			}
		}
		

		if(recordleft==1){

			mc2_fillendtemplate2();
			mc2_fillrendtemplate2();
			mc2_fillendtemplate4();
			mc2_fillrendtemplate4();


			unsigned long int secnum2 = mc2_num_section_endhinge_norder2();
	//		printf("OG secnum2=%d\n", secnum2);
			//Change secnum2 to it's key number.
			i=1;
			found=0;
			while(mc2_sectionkey[i]!=0){
				if(mc2_sectionkey[i]==secnum2){
					found=1;
					secnum2=i;
					break;
				}
				i++;
			}
			if(found==0){	//This shouldn't happen, since all valid section should have been generated during 2-span geneartion.
				printf("another section was generated... problem. Exitting.\n");
				exit(1);
	//			mc2_sectionkey[i]=secnum2;
	//			secnum2=i;
			}
//			printf("NEW secnum2=%lu\n", secnum2);



			unsigned long int secnum4 = mc2_num_section_endhinge_norder4();
	//		printf("OG secnum4=%d\n", secnum4);
			//Change secnum4 to it's key number.
			i=1;
			found=0;
			while(mc2_sectionkey[i]!=0){
				if(mc2_sectionkey[i]==secnum4){
					found=1;
					secnum4=i;
					break;
				}
				i++;
			}
			if(found==0){	//This shouldn't happen, since all valid section should have been generated during 2-span geneartion.
				printf("another section was generated... problem. Exitting.\n");
				exit(1);
	//			mc2_sectionkey[i]=secnum4;
	//			secnum4=i;
			}
//			printf("NEW secnum4=%lu\n", secnum4);




			//Need to check if this left hinge will be a duplicate.

			//check if this endhinge configuration has already been recorded
			duplicate=0;
			test = mc2_firstendhinge[secnum2];
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
				mc2_tot_left_endhinges++;
//				printf("recording left endhinge\n");
				(*mc2_currentendhinge[secnum2]).nextendhinge = mc2_newendhinge(mc2_num_walks, mc2_num_walks2); /* add a mc2_newendhinge to the linked list */

				mc2_currentendhinge[secnum2] = (*mc2_currentendhinge[secnum2]).nextendhinge;
				(*mc2_currentendhinge[secnum2]).sec2 = secnum4;
				(*mc2_currentendhinge[secnum2]).side = 0;
				(*mc2_currentendhinge[secnum2]).num_of_walks = mc2_num_walks;
				(*mc2_currentendhinge[secnum2]).num_of_walks2 = mc2_num_walks2;

				for(i=0; i<=vM*vL-1; i++){
					(*mc2_currentendhinge[secnum2]).hedges[i]=temp_hedges[i];
				}

				for(i=0; i<=mc2_num_walks-1; i++){
					mc2_currentendhinge[secnum2]->start[0][i] = 2;	//change
					mc2_currentendhinge[secnum2]->end[0][i] = 2;		//change
					for(j=1; j<=2; j++){
						mc2_currentendhinge[secnum2]->start[j][i] = mc2_curstart[j][i];
						mc2_currentendhinge[secnum2]->end[j][i] = mc2_curend[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						if(mc2_curwalks[j][i]==1){	//change
							mc2_currentendhinge[secnum2]->walks[j][i]=2;
						}
						else if(mc2_curwalks[j][i]==2){	//change
							mc2_currentendhinge[secnum2]->walks[j][i]=1;
						}
						else{
							mc2_currentendhinge[secnum2]->walks[j][i] = mc2_curwalks[j][i];
						}
					}
				}


				for(i=0; i<=mc2_num_walks2-1; i++){
					mc2_currentendhinge[secnum2]->start2[0][i] = 2;	//change
					mc2_currentendhinge[secnum2]->end2[0][i] = 2;		//change
					for(j=1; j<=2; j++){
						mc2_currentendhinge[secnum2]->start2[j][i] = mc2_curstart2[j][i];
						mc2_currentendhinge[secnum2]->end2[j][i] = mc2_curend2[j][i];
					}
					for(j=0; j<=vM*vL; j++){
						if(mc2_curwalks2[j][i]==1){	//change
							mc2_currentendhinge[secnum2]->walks2[j][i]=2;
						}
						else if(mc2_curwalks2[j][i]==2){	//change
							mc2_currentendhinge[secnum2]->walks2[j][i]=1;
						}
						else{
							mc2_currentendhinge[secnum2]->walks2[j][i] = mc2_curwalks2[j][i];
						}
					}
				}
			}
			else{
				mc2_num_duplicate_left_endhinges++;
			}
		}
		else{	//poly2 was first poly accoridng to lexi ordering
			mc2_num_nonlexi_endhinges++;
		}		
	}
	else{
		mc2_num_duplicate_right_endhinges++;
	}
//	printf("finished recording\n");
	return;
}



///////////////////END OF ENDHINGE FUNCTIONS////////////////////////////////////////////////////////



/***************************************************************************/
struct hinge_span* mc2_newhinge(unsigned int numberofwalks, unsigned int numberofwalks2)		/* returns a pointer to a hinge_span structure */
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
		nextnewhinge->start[i] = mc2_int_vecalloc(0,numberofwalks-1);
		nextnewhinge->end[i] = mc2_int_vecalloc(0,numberofwalks-1);
		nextnewhinge->start2[i] = mc2_int_vecalloc(0,numberofwalks2-1);
		nextnewhinge->end2[i] = mc2_int_vecalloc(0,numberofwalks2-1);
	}
	for(i=0; i<=vM*vL; i++){
		nextnewhinge->walks[i] = mc2_int_vecalloc(0,numberofwalks-1);
		nextnewhinge->walks2[i] = mc2_int_vecalloc(0,numberofwalks2-1);
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
struct endhinge* mc2_newendhinge(unsigned int numberofwalks, unsigned int numberofwalks2)		/* returns a pointer to an endhinge structure */
/* this is a memory allocation function */
{
	struct endhinge *nextnewendhinge;
	unsigned int max_walks = numberofwalks > numberofwalks2 ? numberofwalks : numberofwalks2;

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
		nextnewendhinge->start[i] = mc2_int_vecalloc(0,numberofwalks-1);
		nextnewendhinge->end[i] = mc2_int_vecalloc(0,numberofwalks-1);
		nextnewendhinge->start2[i] = mc2_int_vecalloc(0,max_walks-1);
		nextnewendhinge->end2[i] = mc2_int_vecalloc(0,max_walks-1);
	}
	for(i=0; i<=vM*vL; i++){
		nextnewendhinge->walks[i] = mc2_int_vecalloc(0,numberofwalks-1);
		nextnewendhinge->walks2[i] = mc2_int_vecalloc(0,max_walks-1);
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




void mc2_fillreordertemplate(int ledges, int redges){ //fills mc2_reordertemplate appropriately based on mc2_ordertemplate
	int i, j, ii, jj, s;
	int firstentry, secondentry;
	int connectingedge;

//	printf("mc2_fillreordertemplate called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			mc2_reordertemplate[0][i][j]=0;
			mc2_reordertemplate[1][i][j]=0;
		}
	}

//	printf("initialized\n");
//	printf("ledges=%d, redges=%d\n", ledges, redges);
	for(s=0; s<=1; s++){	//for each side
		firstentry=1;
		for(i=0; i<=M; i++){	
			for(j=0; j<=L; j++){
				if(mc2_ordertemplate[s][i][j] % 2 != 0 && mc2_reordertemplate[s][i][j]==0){	//if odd
//					printf("odd edge found in (%d,%d)\n", j, i);
					mc2_reordertemplate[s][i][j] = firstentry;
//					printf("recorded edge %d in reorder\n", firstentry);
					if(mc2_ordertemplate[s][i][j]==1){	//if 1, then have to find the location of ledges or redges.
						if(s==0){
							connectingedge=ledges;
						}
						else{
							connectingedge=redges;
						}
					}
					else{	//if not 1, then find location the edgenum-1
						connectingedge=mc2_ordertemplate[s][i][j]-1;
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
							if(mc2_ordertemplate[s][ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
								mc2_reordertemplate[s][ii][jj] = secondentry;
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
				else if(mc2_ordertemplate[s][i][j]>0 && mc2_reordertemplate[s][i][j]==0){	// if even
//					printf("even edge found in (%d,%d)\n", j, i);
					mc2_reordertemplate[s][i][j] = firstentry;
//					printf("recorded edge %d in reorder\n", firstentry);
					if(s==0 && mc2_ordertemplate[s][i][j]==ledges){	//if max, then have to find the location of edge 1.
						connectingedge=1;
					}
					else if(s==1 && mc2_ordertemplate[s][i][j]==redges){
						connectingedge=1;
					}
					else{	//if not max, then find location the edgenum+1
						connectingedge=mc2_ordertemplate[s][i][j]+1;
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
							if(mc2_ordertemplate[s][ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
								mc2_reordertemplate[s][ii][jj] = secondentry;
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

void mc2_fillreordertemplate2(int ledges, int redges){ //fills mc2_reordertemplate2 appropriately based on mc2_ordertemplate
	int i, j, ii, jj, s;
	int firstentry, secondentry;
	int connectingedge;

//	printf("mc2_fillreordertemplate2 called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			mc2_reordertemplate2[0][i][j]=0;
			mc2_reordertemplate2[1][i][j]=0;
		}
	}

//	printf("initialized\n");
//	printf("ledges=%d, redges=%d\n", ledges, redges);
	for(s=0; s<=1; s++){	//for each side
		firstentry=1;
		for(i=0; i<=M; i++){	
			for(j=0; j<=L; j++){
				if(mc2_ordertemplate2[s][i][j] % 2 != 0 && mc2_reordertemplate2[s][i][j]==0){	//if odd
//					printf("odd edge found in (%d,%d)\n", j, i);
					mc2_reordertemplate2[s][i][j] = firstentry;
//					printf("recorded edge %d in reorder\n", firstentry);
					if(mc2_ordertemplate2[s][i][j]==1){	//if 1, then have to find the location of ledges or redges.
						if(s==0){
							connectingedge=ledges;
						}
						else{
							connectingedge=redges;
						}
					}
					else{	//if not 1, then find location the edgenum-1
						connectingedge=mc2_ordertemplate2[s][i][j]-1;
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
							if(mc2_ordertemplate2[s][ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
								mc2_reordertemplate2[s][ii][jj] = secondentry;
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
				else if(mc2_ordertemplate2[s][i][j]>0 && mc2_reordertemplate2[s][i][j]==0){	// if even
//					printf("even edge found in (%d,%d)\n", j, i);
					mc2_reordertemplate2[s][i][j] = firstentry;
//					printf("recorded edge %d in reorder\n", firstentry);
					if(s==0 && mc2_ordertemplate2[s][i][j]==ledges){	//if max, then have to find the location of edge 1.
						connectingedge=1;
					}
					else if(s==1 && mc2_ordertemplate2[s][i][j]==redges){
						connectingedge=1;
					}
					else{	//if not max, then find location the edgenum+1
						connectingedge=mc2_ordertemplate2[s][i][j]+1;
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
							if(mc2_ordertemplate2[s][ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
								mc2_reordertemplate2[s][ii][jj] = secondentry;
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

void mc2_printordtemp(){
	int i, j;
	printf("printing mc2_ordertemplate. Side 0:\n");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("%d ", mc2_ordertemplate[0][i][j]);
		}
		printf("\n");
	}

	printf("mc2_ordertemplate side1:\n");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("%d ", mc2_ordertemplate[1][i][j]);
		}
		printf("\n");
	}
}

void mc2_printreordtemp(){
	int i, j;
	printf("printing mc2_reordertemplate. Side 0:\n");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("%d ", mc2_reordertemplate[0][i][j]);
		}
		printf("\n");
	}

	printf("mc2_reordertemplate side1:\n");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("%d ", mc2_reordertemplate[1][i][j]);
		}
		printf("\n");
	}
}



void mc2_fillrendtemplate(){ //fills mc2_rendtemplate appropriately based on mc2_endtemplate
	int i, j, ii, jj;
	int firstentry, secondentry;
	int connectingedge;

//	printf("mc2_fillrendtemplate called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			mc2_rendtemplate[i][j]=0;
			mc2_rendtemplate[i][j]=0;
		}
	}

//	printf("initialized\n");
	firstentry=1;
	for(i=0; i<=M; i++){	
		for(j=0; j<=L; j++){
			if(mc2_endtemplate[i][j] % 2 != 0 && mc2_rendtemplate[i][j]==0){	//if odd
//				printf("odd edge found in (%d,%d)\n", j, i);
				mc2_rendtemplate[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(mc2_endtemplate[i][j]==1){	//if 1, then have to find the location of ledges or redges.
					connectingedge=mc2_EndOrdNum[0]-1;
				}
				else{	//if not 1, then find location the edgenum-1
					connectingedge=mc2_endtemplate[i][j]-1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=mc2_EndOrdNum[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(mc2_endtemplate[ii][jj]==connectingedge){
//							printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							mc2_rendtemplate[ii][jj] = secondentry;
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
			else if(mc2_endtemplate[i][j]>0 && mc2_rendtemplate[i][j]==0){	// if even
//				printf("even edge found in (%d,%d)\n", j, i);
				mc2_rendtemplate[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(mc2_endtemplate[i][j]==mc2_EndOrdNum[0]-1){	//if max, then have to find the location of edge 1.
					connectingedge=1;
				}
				else{	//if not max, then find location the edgenum+1
					connectingedge=mc2_endtemplate[i][j]+1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=mc2_EndOrdNum[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(mc2_endtemplate[ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							mc2_rendtemplate[ii][jj] = secondentry;
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




void mc2_fillrendtemplate3(){ //fills mc2_rendtemplate3 appropriately based on mc2_endtemplate3
	int i, j, ii, jj;
	int firstentry, secondentry;
	int connectingedge;

//	printf("mc2_fillrendtemplate called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			mc2_rendtemplate3[i][j]=0;
			mc2_rendtemplate3[i][j]=0;
		}
	}

//	printf("initialized\n");
	firstentry=1;
	for(i=0; i<=M; i++){	
		for(j=0; j<=L; j++){
			if(mc2_endtemplate3[i][j] % 2 != 0 && mc2_rendtemplate3[i][j]==0){	//if odd
//				printf("odd edge found in (%d,%d)\n", j, i);
				mc2_rendtemplate3[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(mc2_endtemplate3[i][j]==1){	//if 1, then have to find the location of ledges or redges.
					connectingedge=mc2_EndOrdNum2[0]-1;
				}
				else{	//if not 1, then find location the edgenum-1
					connectingedge=mc2_endtemplate3[i][j]-1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=mc2_EndOrdNum2[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(mc2_endtemplate3[ii][jj]==connectingedge){
//							printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							mc2_rendtemplate3[ii][jj] = secondentry;
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
			else if(mc2_endtemplate3[i][j]>0 && mc2_rendtemplate3[i][j]==0){	// if even
//				printf("even edge found in (%d,%d)\n", j, i);
				mc2_rendtemplate3[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(mc2_endtemplate3[i][j]==mc2_EndOrdNum2[0]-1){	//if max, then have to find the location of edge 1.
					connectingedge=1;
				}
				else{	//if not max, then find location the edgenum+1
					connectingedge=mc2_endtemplate3[i][j]+1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=mc2_EndOrdNum2[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(mc2_endtemplate3[ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							mc2_rendtemplate3[ii][jj] = secondentry;
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




void mc2_fillendtemplate2(){ //fills mc2_endtemplate2 (left endhinge) appropriately based on mc2_endtemplate (right endhinge)
	int i, j;
	int max=0;
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			mc2_endtemplate2[i][j] = 0;
			if(mc2_endtemplate[i][j]>max){
				max = mc2_endtemplate[i][j];
			}
		}
	}
	
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			if(mc2_endtemplate[i][j]>0 && mc2_endtemplate[i][j] < max){
				mc2_endtemplate2[i][j] = mc2_endtemplate[i][j] + 1;
			}
			else if(mc2_endtemplate[i][j]==max){
				mc2_endtemplate2[i][j]=1;
			}
		}
	}
}






void mc2_fillendtemplate4(){ //fills mc2_endtemplate4 (left endhinge) appropriately based on mc2_endtemplate3 (right endhinge)
	int i, j;
	int max=0;
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			mc2_endtemplate4[i][j] = 0;
			if(mc2_endtemplate3[i][j]>max){
				max = mc2_endtemplate3[i][j];
			}
		}
	}
	
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			if(mc2_endtemplate3[i][j]>0 && mc2_endtemplate3[i][j] < max){
				mc2_endtemplate4[i][j] = mc2_endtemplate3[i][j] + 1;
			}
			else if(mc2_endtemplate3[i][j]==max){
				mc2_endtemplate4[i][j]=1;
			}
		}
	}
}




void mc2_fillrendtemplate2(){ //fills mc2_rendtemplate2 appropriately based on mc2_endtemplate2
	int i, j, ii, jj;
	int firstentry, secondentry;
	int connectingedge;

//	printf("mc2_fillrendtemplate2 called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			mc2_rendtemplate2[i][j]=0;
			mc2_rendtemplate2[i][j]=0;
		}
	}

//	printf("initialized\n");
	firstentry=1;
	for(i=0; i<=M; i++){	
		for(j=0; j<=L; j++){
			if(mc2_endtemplate2[i][j] % 2 != 0 && mc2_rendtemplate2[i][j]==0){	//if odd
//				printf("odd edge found in (%d,%d)\n", j, i);
				mc2_rendtemplate2[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(mc2_endtemplate2[i][j]==1){	//if 1, then have to find the location of ledges or redges.
					connectingedge=mc2_EndOrdNum[0]-1;
				}
				else{	//if not 1, then find location the edgenum-1
					connectingedge=mc2_endtemplate2[i][j]-1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=mc2_EndOrdNum[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(mc2_endtemplate2[ii][jj]==connectingedge){
//							printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							mc2_rendtemplate2[ii][jj] = secondentry;
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
			else if(mc2_endtemplate2[i][j]>0 && mc2_rendtemplate2[i][j]==0){	// if even
//				printf("even edge found in (%d,%d)\n", j, i);
				mc2_rendtemplate2[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(mc2_endtemplate2[i][j]==mc2_EndOrdNum[0]-1){	//if max, then have to find the location of edge 1.
					connectingedge=1;
				}
				else{	//if not max, then find location the edgenum+1
					connectingedge=mc2_endtemplate2[i][j]+1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=mc2_EndOrdNum[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(mc2_endtemplate2[ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							mc2_rendtemplate2[ii][jj] = secondentry;
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




void mc2_fillrendtemplate4(){ //fills mc2_rendtemplate4 appropriately based on mc2_endtemplate4
	int i, j, ii, jj;
	int firstentry, secondentry;
	int connectingedge;

//	printf("mc2_fillrendtemplate4 called\n");
	for(i=0; i<=M; i++){	//initialize to zero
		for(j=0; j<=L; j++){
			mc2_rendtemplate4[i][j]=0;
			mc2_rendtemplate4[i][j]=0;
		}
	}

//	printf("initialized\n");
	firstentry=1;
	for(i=0; i<=M; i++){	
		for(j=0; j<=L; j++){
			if(mc2_endtemplate4[i][j] % 2 != 0 && mc2_rendtemplate4[i][j]==0){	//if odd
//				printf("odd edge found in (%d,%d)\n", j, i);
				mc2_rendtemplate4[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(mc2_endtemplate4[i][j]==1){	//if 1, then have to find the location of ledges or redges.
					connectingedge=mc2_EndOrdNum2[0]-1;
				}
				else{	//if not 1, then find location the edgenum-1
					connectingedge=mc2_endtemplate4[i][j]-1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=mc2_EndOrdNum2[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(mc2_endtemplate4[ii][jj]==connectingedge){
//							printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							mc2_rendtemplate4[ii][jj] = secondentry;
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
			else if(mc2_endtemplate4[i][j]>0 && mc2_rendtemplate4[i][j]==0){	// if even
//				printf("even edge found in (%d,%d)\n", j, i);
				mc2_rendtemplate4[i][j] = firstentry;
//				printf("recorded edge %d in rend\n", firstentry);
				if(mc2_endtemplate4[i][j]==mc2_EndOrdNum2[0]-1){	//if max, then have to find the location of edge 1.
					connectingedge=1;
				}
				else{	//if not max, then find location the edgenum+1
					connectingedge=mc2_endtemplate4[i][j]+1;
				}

				if(firstentry==1){	//if we putting in the first edge, the second number must be the max number.
					secondentry=mc2_EndOrdNum2[0]-1;
				}
				else{
					secondentry=firstentry+1;
				}
				for(ii=0; ii<=M; ii++){
					for(jj=0; jj<=L; jj++){
						if(mc2_endtemplate4[ii][jj]==connectingedge){
//								printf("found next edge (%d). located at (%d, %d)\n", connectingedge, jj, ii);
							mc2_rendtemplate4[ii][jj] = secondentry;
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

void mc2_printbuiltwalks(){
	int i, j;
	printf("There are %d walks in built_walks:\n", mc2_num_built_walks);
	for(i=0; i<=mc2_num_built_walks-1; i++){
		printf("Start: (%d, %d, %d). ", mc2_built_walks_start[i][0], mc2_built_walks_start[i][1], mc2_built_walks_start[i][2]);
		printf("Direcs:");
		j=0;
		while(mc2_built_walks_direcs[i][j]!=0){
			printf(" %d", mc2_built_walks_direcs[i][j]);
			j++;
		}
		printf(". End: (%d, %d, %d)\n", mc2_built_walks_end[i][0], mc2_built_walks_end[i][1], mc2_built_walks_end[i][2]);
	}
}

void mc2_printbuiltwalks_all(){
	int i, j;
	printf("There are %d walks in built_walks (printingall):\n", mc2_num_built_walks);
	for(i=0; i<=mc2_num_built_walks-1; i++){
		printf("Start: (%d, %d, %d). ", mc2_built_walks_start[i][0], mc2_built_walks_start[i][1], mc2_built_walks_start[i][2]);
		printf("Direcs:");
		for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
			printf(" %d", mc2_built_walks_direcs[i][j]);
		}
		printf(". End: (%d, %d, %d)\n", mc2_built_walks_end[i][0], mc2_built_walks_end[i][1], mc2_built_walks_end[i][2]);
	}
}

void mc2_printbuiltwalks2_all(){
	int i, j;
	printf("There are %d walks in built_walks2 (printingall):\n", mc2_num_built_walks2);
	for(i=0; i<=mc2_num_built_walks2-1; i++){
		printf("Start: (%d, %d, %d). ", mc2_built_walks_start2[i][0], mc2_built_walks_start2[i][1], mc2_built_walks_start2[i][2]);
		printf("Direcs:");
		for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
			printf(" %d", mc2_built_walks_direcs2[i][j]);
		}
		printf(". End: (%d, %d, %d)\n", mc2_built_walks_end2[i][0], mc2_built_walks_end2[i][1], mc2_built_walks_end2[i][2]);
	}
}

void mc2_add_to_built_walks(unsigned long int secnum, int nth_tspan){
	//adds the walks in the tspan sec1->nth_tspan to built_walks

	//3 main cases for a walk in a 2span.
	//case 1a. starts left, ends right --> just append walk onto appropriate walk. mc2_num_built_walks doesn't change.
	//case 1b. starts right, ends left --> just append walk onto appropriate walk. mc2_num_built_walks doesn't change.
	//case 2. starts left, ends left --> will connect two pre-existing walks. mc2_num_built_walks decreases by one.
	//case 3. starts right, ends right --> will add a new walk to built_walks. mc2_num_built_walks increases by one.

	int i,j;

	int walknum;
	int walktoadd[MAX_vMvL+1];
	int walktoadd_start[3];
	int walktoadd_end[3];

	int built_walknum;
	int length;

	//FIRST POLY
	for(walknum=0; walknum<=mc2_t_num_walks[secnum][nth_tspan]-1; walknum++){ //for each walk in the t-span
		//fill walktoadd
		for(i=0; i<=2; i++){
			walktoadd_start[i]=mc2_t_start[secnum][nth_tspan][i][walknum];
			walktoadd_end[i]=mc2_t_end[secnum][nth_tspan][i][walknum];
		}
		length=0;
		for(i=0; i<=vM*vL; i++){
			walktoadd[i] = mc2_t_walks[secnum][nth_tspan][i][walknum];
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
			for(built_walknum=0; built_walknum<=mc2_num_built_walks-1; built_walknum++){
				if(mc2_built_walks_start[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_start[built_walknum][2]==walktoadd_start[2]){
					//built_walk starts at same place walktoadd starts (going in opposite direction)
					//add reversed version of walktoadd to the start of built_walks

					//shift mc2_built_walks_direcs up "length-1" entries.
					i=0;
					while(mc2_built_walks_direcs[built_walknum][i+1]!=0){
						i++;
					}
					while(i>=0){
						mc2_built_walks_direcs[built_walknum][i+length-1] = mc2_built_walks_direcs[built_walknum][i];
						i--;
					}
					//entries 0->length-2 are bogus now. Need to fill them with reverse version of walktoadd
					for(i=0; i<=length-2; i++){
						mc2_built_walks_direcs[built_walknum][i] = mc2_reverse_direc(walktoadd[length-1-i]);
					}

					//change start point
					mc2_built_walks_start[built_walknum][1] = walktoadd_end[1];
					mc2_built_walks_start[built_walknum][2] = walktoadd_end[2];
				}
				else if(mc2_built_walks_end[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_end[built_walknum][2]==walktoadd_start[2]){
					//built_walk ends at same place walktoadd starts (going in same direction)

					i=0;
					while(mc2_built_walks_direcs[built_walknum][i]!=0){
						i++;
					}
					j=1;	//Note: skip overlapping first edge.
					while(walktoadd[j]!=0){
						mc2_built_walks_direcs[built_walknum][i] = walktoadd[j];
						i++;
						j++;
					}

					//change end point
					mc2_built_walks_end[built_walknum][1] = walktoadd_end[1];
					mc2_built_walks_end[built_walknum][2] = walktoadd_end[2];
				}
			}
		}
		else if(walktoadd_start[0]==2 && walktoadd_end[0]==0){	//case1b
//			printf("CASE 1b\n");
			//find which walk to append to
			for(built_walknum=0; built_walknum<=mc2_num_built_walks-1; built_walknum++){
				if(mc2_built_walks_start[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_start[built_walknum][2]==walktoadd_end[2]){
					//built_walk starts at same place walktoadd ends (going in same direction)
					//add walktoadd to the start of built_walks

					//shift mc2_built_walks_direcs up "length-1" entries.
					i=0;
					while(mc2_built_walks_direcs[built_walknum][i+1]!=0){
						i++;
					}
					while(i>=0){
						mc2_built_walks_direcs[built_walknum][i+length-1] = mc2_built_walks_direcs[built_walknum][i];
						i--;
					}
					//entries 0->length-2 are bogus now. Need to fill them with walktoadd
					for(i=0; i<=length-2; i++){
						mc2_built_walks_direcs[built_walknum][i] = walktoadd[i];
					}

					//change start point
					mc2_built_walks_start[built_walknum][1] = walktoadd_start[1];
					mc2_built_walks_start[built_walknum][2] = walktoadd_start[2];
				}
				else if(mc2_built_walks_end[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_end[built_walknum][2]==walktoadd_end[2]){
					//built_walk ends at same place walktoadd ends (going in opposite direction)

					i=0;
					while(mc2_built_walks_direcs[built_walknum][i]!=0){
						i++;
					}
					j=length-2; //Note: skip overlapping first edge.
					while(j>=0){
						mc2_built_walks_direcs[built_walknum][i] = mc2_reverse_direc(walktoadd[j]);
						i++;
						j--;
					}

					//change end point
					mc2_built_walks_end[built_walknum][1] = walktoadd_start[1];
					mc2_built_walks_end[built_walknum][2] = walktoadd_start[2];
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
			//mc2_num_built_walks--

			int a=-1;
			int b=-1;
			int newwalk_start[3];
			int newwalk_end[3];
			int newwalk[MAX_vMvL*(MAX_SPAN+1)];

			//find a and b.
			for(built_walknum=0; built_walknum<=mc2_num_built_walks-1; built_walknum++){
				if(mc2_built_walks_start[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_start[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}
				else if(mc2_built_walks_end[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_end[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}

				else if(mc2_built_walks_start[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_start[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
				else if(mc2_built_walks_end[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_end[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
			}


			if(a<0 || b<0){	//invalid--will close off a poly. return 0 (which means fail).
				printf("invalid connection. exitting\n");
				exit(1);
			}

			//fill newwalk with proper version of a. (start at opposite end of where walktoend starts)
			if(mc2_built_walks_end[a][1]==walktoadd_start[1] && mc2_built_walks_end[a][2]==walktoadd_start[2]){
				//correct direction. easy.
				for(i=0; i<=2; i++){
					newwalk_start[i]=mc2_built_walks_start[a][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					newwalk[i] = mc2_built_walks_direcs[a][i];
				}
			}
			else{
				//need to reverse a
				for(i=0; i<=2; i++){
					newwalk_start[i]=mc2_built_walks_end[a][i];
				}
				i=0;
				while(mc2_built_walks_direcs[a][i+1]!=0){
					i++;
				}
				j=0;
				while(i>=0){
					newwalk[j] = mc2_reverse_direc(mc2_built_walks_direcs[a][i]);
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
			if(mc2_built_walks_start[b][1]==walktoadd_end[1] && mc2_built_walks_start[b][2]==walktoadd_end[2]){
				//correct direction. easy.
				//append b to newwalk. //i is already in the correct position.
				j=1;	//Note: skip overlapping first edge.
				while(mc2_built_walks_direcs[b][j]!=0){
					newwalk[i] = mc2_built_walks_direcs[b][j];
					i++;
					j++;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = mc2_built_walks_end[b][i];
				}
			}
			else{
				//need to append reversed version of b
				j=0;
				while(mc2_built_walks_direcs[b][j+1]!=0){
					j++;
				}
				j--; //Note: skip overlapping first edge.
				while(j>=0){
					newwalk[i] = mc2_reverse_direc(mc2_built_walks_direcs[b][j]);
					i++;
					j--;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = mc2_built_walks_start[b][i];
				}
			}

			//replace a with newwalk.
			for(i=0; i<=2; i++){
				mc2_built_walks_start[a][i] = newwalk_start[i];
				mc2_built_walks_end[a][i] = newwalk_end[i];
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				mc2_built_walks_direcs[a][i] = newwalk[i];
			}

			//move walk b to end of array
			while(b<mc2_num_built_walks-1){
				for(i=0; i<=2; i++){
					mc2_built_walks_start[b][i] = mc2_built_walks_start[b+1][i];
					mc2_built_walks_end[b][i] = mc2_built_walks_end[b+1][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					mc2_built_walks_direcs[b][i] = mc2_built_walks_direcs[b+1][i];
				}
				b++;
			}

			//delete b
			for(i=0; i<=2; i++){
				mc2_built_walks_start[b][i] = -1;
				mc2_built_walks_end[b][i] = -1;
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				mc2_built_walks_direcs[b][i] = 0;
			}

			mc2_num_built_walks--;
		}
		else if(walktoadd_start[0]==2 && walktoadd_end[0]==2){	//case3
//			printf("CASE 3\n");
			//add a new walk built_walks (walktoadd)
			mc2_num_built_walks++;
			for(i=0; i<=2; i++){
				mc2_built_walks_start[mc2_num_built_walks-1][i] = walktoadd_start[i];
				mc2_built_walks_end[mc2_num_built_walks-1][i] = walktoadd_end[i];
			}
			for(i=0; i<=vM*vL; i++){
				mc2_built_walks_direcs[mc2_num_built_walks-1][i] = walktoadd[i];
			}
		}
		else{
			printf("Problem when adding a tspan. Exitting\n");
			exit(1);
		}
//		printf("one walk has been added, here is the result:\n");
//		mc2_printbuiltwalks_all();
	}


////////////////////////////////////////////

	//SECOND POLY
	for(walknum=0; walknum<=mc2_t_num_walks2[secnum][nth_tspan]-1; walknum++){ //for each walk in the t-span
		//fill walktoadd
		for(i=0; i<=2; i++){
			walktoadd_start[i]=mc2_t_start2[secnum][nth_tspan][i][walknum];
			walktoadd_end[i]=mc2_t_end2[secnum][nth_tspan][i][walknum];
		}
		length=0;
		for(i=0; i<=vM*vL; i++){
			walktoadd[i] = mc2_t_walks2[secnum][nth_tspan][i][walknum];
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
			for(built_walknum=0; built_walknum<=mc2_num_built_walks2-1; built_walknum++){
				if(mc2_built_walks_start2[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_start2[built_walknum][2]==walktoadd_start[2]){
					//built_walk2 starts at same place walktoadd starts (going in opposite direction)
					//add reversed version of walktoadd to the start of built_walks2

					//shift mc2_built_walks_direcs2 up "length-1" entries.
					i=0;
					while(mc2_built_walks_direcs2[built_walknum][i+1]!=0){
						i++;
					}
					while(i>=0){
						mc2_built_walks_direcs2[built_walknum][i+length-1] = mc2_built_walks_direcs2[built_walknum][i];
						i--;
					}
					//entries 0->length-2 are bogus now. Need to fill them with reverse version of walktoadd
					for(i=0; i<=length-2; i++){
						mc2_built_walks_direcs2[built_walknum][i] = mc2_reverse_direc(walktoadd[length-1-i]);
					}

					//change start point
					mc2_built_walks_start2[built_walknum][1] = walktoadd_end[1];
					mc2_built_walks_start2[built_walknum][2] = walktoadd_end[2];
				}
				else if(mc2_built_walks_end2[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_end2[built_walknum][2]==walktoadd_start[2]){
					//built_walk2 ends at same place walktoadd starts (going in same direction)

					i=0;
					while(mc2_built_walks_direcs2[built_walknum][i]!=0){
						i++;
					}
					j=1;	//Note: skip overlapping first edge.
					while(walktoadd[j]!=0){
						mc2_built_walks_direcs2[built_walknum][i] = walktoadd[j];
						i++;
						j++;
					}

					//change end point
					mc2_built_walks_end2[built_walknum][1] = walktoadd_end[1];
					mc2_built_walks_end2[built_walknum][2] = walktoadd_end[2];
				}
			}
		}
		else if(walktoadd_start[0]==2 && walktoadd_end[0]==0){	//case1b
//			printf("CASE 1b\n");
			//find which walk to append to
			for(built_walknum=0; built_walknum<=mc2_num_built_walks2-1; built_walknum++){
				if(mc2_built_walks_start2[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_start2[built_walknum][2]==walktoadd_end[2]){
					//built_walk2 starts at same place walktoadd ends (going in same direction)
					//add walktoadd to the start of built_walks2

					//shift mc2_built_walks_direcs2 up "length-1" entries.
					i=0;
					while(mc2_built_walks_direcs2[built_walknum][i+1]!=0){
						i++;
					}
					while(i>=0){
						mc2_built_walks_direcs2[built_walknum][i+length-1] = mc2_built_walks_direcs2[built_walknum][i];
						i--;
					}
					//entries 0->length-2 are bogus now. Need to fill them with walktoadd
					for(i=0; i<=length-2; i++){
						mc2_built_walks_direcs2[built_walknum][i] = walktoadd[i];
					}

					//change start point
					mc2_built_walks_start2[built_walknum][1] = walktoadd_start[1];
					mc2_built_walks_start2[built_walknum][2] = walktoadd_start[2];
				}
				else if(mc2_built_walks_end2[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_end2[built_walknum][2]==walktoadd_end[2]){
					//built_walk ends at same place walktoadd ends (going in opposite direction)

					i=0;
					while(mc2_built_walks_direcs2[built_walknum][i]!=0){
						i++;
					}
					j=length-2; //Note: skip overlapping first edge.
					while(j>=0){
						mc2_built_walks_direcs2[built_walknum][i] = mc2_reverse_direc(walktoadd[j]);
						i++;
						j--;
					}

					//change end point
					mc2_built_walks_end2[built_walknum][1] = walktoadd_start[1];
					mc2_built_walks_end2[built_walknum][2] = walktoadd_start[2];
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
			//mc2_num_built_walks2--

			int a=-1;
			int b=-1;
			int newwalk_start[3];
			int newwalk_end[3];
			int newwalk[MAX_vMvL*(MAX_SPAN+1)];

			//find a and b.
			for(built_walknum=0; built_walknum<=mc2_num_built_walks2-1; built_walknum++){
				if(mc2_built_walks_start2[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_start2[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}
				else if(mc2_built_walks_end2[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_end2[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}

				else if(mc2_built_walks_start2[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_start2[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
				else if(mc2_built_walks_end2[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_end2[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
			}

//			printf("a=%d,b=%d\n", a, b);

			if(a<0 || b<0){	//invalid--will close off a seperate poly that doesn't span tube. return 0 (which means fail).
				printf("invalid connection (2nd poly). exitting\n");
				exit(1);
			}

			//fill newwalk with proper version of a. (start at opposite end of where walktoend starts)
			if(mc2_built_walks_end2[a][1]==walktoadd_start[1] && mc2_built_walks_end2[a][2]==walktoadd_start[2]){
				//correct direction. easy.
				for(i=0; i<=2; i++){
					newwalk_start[i]=mc2_built_walks_start2[a][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					newwalk[i] = mc2_built_walks_direcs2[a][i];
				}
			}
			else{
				//need to reverse a
				for(i=0; i<=2; i++){
					newwalk_start[i]=mc2_built_walks_end2[a][i];
				}
				i=0;
				while(mc2_built_walks_direcs2[a][i+1]!=0){
					i++;
				}
				j=0;
				while(i>=0){
					newwalk[j] = mc2_reverse_direc(mc2_built_walks_direcs2[a][i]);
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
			if(mc2_built_walks_start2[b][1]==walktoadd_end[1] && mc2_built_walks_start2[b][2]==walktoadd_end[2]){
				//correct direction. easy.
				//append b to newwalk. //i is already in the correct position.
				j=1;	//Note: skip overlapping first edge.
				while(mc2_built_walks_direcs2[b][j]!=0){
					newwalk[i] = mc2_built_walks_direcs2[b][j];
					i++;
					j++;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = mc2_built_walks_end2[b][i];
				}
			}
			else{
				//need to append reversed version of b
				j=0;
				while(mc2_built_walks_direcs2[b][j+1]!=0){
					j++;
				}
				j--; //Note: skip overlapping first edge.
				while(j>=0){
					newwalk[i] = mc2_reverse_direc(mc2_built_walks_direcs2[b][j]);
					i++;
					j--;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = mc2_built_walks_start2[b][i];
				}
			}

			//replace a with newwalk.
			for(i=0; i<=2; i++){
				mc2_built_walks_start2[a][i] = newwalk_start[i];
				mc2_built_walks_end2[a][i] = newwalk_end[i];
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				mc2_built_walks_direcs2[a][i] = newwalk[i];
			}

			//move walk b to end of array
			while(b<mc2_num_built_walks2-1){
				for(i=0; i<=2; i++){
					mc2_built_walks_start2[b][i] = mc2_built_walks_start2[b+1][i];
					mc2_built_walks_end2[b][i] = mc2_built_walks_end2[b+1][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					mc2_built_walks_direcs2[b][i] = mc2_built_walks_direcs2[b+1][i];
				}
				b++;
			}

			//delete b
			for(i=0; i<=2; i++){
				mc2_built_walks_start2[b][i] = -1;
				mc2_built_walks_end2[b][i] = -1;
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				mc2_built_walks_direcs2[b][i] = 0;
			}

			mc2_num_built_walks2--;
		}
		else if(walktoadd_start[0]==2 && walktoadd_end[0]==2){	//case3
//			printf("CASE 3\n");
			//add a new walk built_walks2 (walktoadd)
			mc2_num_built_walks2++;
			for(i=0; i<=2; i++){
				mc2_built_walks_start2[mc2_num_built_walks2-1][i] = walktoadd_start[i];
				mc2_built_walks_end2[mc2_num_built_walks2-1][i] = walktoadd_end[i];
			}
			for(i=0; i<=vM*vL; i++){
				mc2_built_walks_direcs2[mc2_num_built_walks2-1][i] = walktoadd[i];
			}
		}
		else{
			printf("Problem when adding a tspan. Exitting\n");
			exit(1);
		}
//		printf("one walk has been added, here is the result:\n");
//		mc2_printbuiltwalks2_all();
	}


}

void mc2_add_right_endhinge(unsigned long int secnum, int nth_endhinge){
	//Each walk in the right endhinge will connect two existing walks. walk->walktoadd->walk.
	//Done until there is only one walk left.

	//ALGORITHM:
	//For mc2_num_built_walks>1:
		//Find walk numbers that will connect (call them a,b) a connects to walktoadd_start, b connects to walktoadd_end.
		//Fill array "newwalk" which will hold the final product. must cover all cases for ordering.
			//fill with a, then walktoadd, then b.
		//replace a with newwalk.
		//move walk b to end of array, then delete
		//mc2_num_built_walks--

	//When mc2_num_built_walks==1, finish up SAP.

	int i,j;

	int walknum;
	int walktoadd[MAX_vMvL+1];
	int walktoadd_start[3];
	int walktoadd_end[3];

	int built_walknum;
	int length;

	int a=-1;
	int b=-1;
	int newwalk_start[3];
	int newwalk_end[3];
	int newwalk[MAX_vMvL*(MAX_SPAN+1)];

	for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
		newwalk[i]=0;
	}
/*
	printf("mc2_Rend_num_walks=%d\n", mc2_Rend_num_walks[secnum][nth_endhinge]);
	printf("mc2_num_built_walks=%d\n", mc2_num_built_walks);
*/
	for(walknum=0; walknum<=mc2_Rend_num_walks[secnum][nth_endhinge]-1; walknum++){ //for each walk in the right endhinge
		//fill walktoadd
		for(i=0; i<=2; i++){
			walktoadd_start[i]=mc2_Rend_start[secnum][nth_endhinge][i][walknum];
			walktoadd_end[i]=mc2_Rend_end[secnum][nth_endhinge][i][walknum];
		}
		length=0;
		for(i=0; i<=vM*vL; i++){
			walktoadd[i] = mc2_Rend_walks[secnum][nth_endhinge][i][walknum];
			if(walktoadd[i]!=0){
				length++;
			}
		}

		if(mc2_num_built_walks>1){
			//find a and b.
			for(built_walknum=0; built_walknum<=mc2_num_built_walks-1; built_walknum++){
				if(mc2_built_walks_start[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_start[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}
				else if(mc2_built_walks_end[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_end[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}

				else if(mc2_built_walks_start[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_start[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
				else if(mc2_built_walks_end[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_end[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
			}

			if(a<0 || b<0){
				printf("didn't find a or b. Exitting\n");
				exit(1);
			}
			//fill newwalk with proper version of a. (start at opposite end of where walktoend starts)
			if(mc2_built_walks_end[a][1]==walktoadd_start[1] && mc2_built_walks_end[a][2]==walktoadd_start[2]){
				//correct direction. easy.
				for(i=0; i<=2; i++){
					newwalk_start[i]=mc2_built_walks_start[a][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					newwalk[i] = mc2_built_walks_direcs[a][i];
				}
			}
			else{
				//need to reverse a
				for(i=0; i<=2; i++){
					newwalk_start[i]=mc2_built_walks_end[a][i];
				}
				i=0;
				while(mc2_built_walks_direcs[a][i+1]!=0){
					i++;
				}
				j=0;
				while(i>=0){
					newwalk[j] = mc2_reverse_direc(mc2_built_walks_direcs[a][i]);
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
			if(mc2_built_walks_start[b][1]==walktoadd_end[1] && mc2_built_walks_start[b][2]==walktoadd_end[2]){
				//correct direction. easy.
				//append b to newwalk. //i is already in the correct position.
				j=1;	//Note: skip overlapping first edge.
				while(mc2_built_walks_direcs[b][j]!=0){
					newwalk[i] = mc2_built_walks_direcs[b][j];
					i++;
					j++;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = mc2_built_walks_end[b][i];
				}
			}
			else{
				//need to append reversed version of b
				j=0;
				while(mc2_built_walks_direcs[b][j+1]!=0){
					j++;
				}
				j--; //Note: skip overlapping first edge.
				while(j>=0){
					newwalk[i] = mc2_reverse_direc(mc2_built_walks_direcs[b][j]);
					i++;
					j--;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = mc2_built_walks_start[b][i];
				}
			}

			//replace a with newwalk.
			for(i=0; i<=2; i++){
				mc2_built_walks_start[a][i] = newwalk_start[i];
				mc2_built_walks_end[a][i] = newwalk_end[i];
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				mc2_built_walks_direcs[a][i] = newwalk[i];
			}

			//move walk b to end of array
			while(b<mc2_num_built_walks-1){
				for(i=0; i<=2; i++){
					mc2_built_walks_start[b][i] = mc2_built_walks_start[b+1][i];
					mc2_built_walks_end[b][i] = mc2_built_walks_end[b+1][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					mc2_built_walks_direcs[b][i] = mc2_built_walks_direcs[b+1][i];
				}
				b++;
			}

			//delete b
			for(i=0; i<=2; i++){
				mc2_built_walks_start[b][i] = -1;
				mc2_built_walks_end[b][i] = -1;
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				mc2_built_walks_direcs[b][i] = 0;
			}
			mc2_num_built_walks--;
		}
		else{	//close up SAP.
/*
printf("\ntime to close up the SAP:\n");
mc2_printbuiltwalks_all();
printf("walktoadd:");
i=0;
while(walktoadd[i]!=0){
	printf(" %d", walktoadd[i]);
	i++;
}
printf("\n");
*/			//either append walktoadd in same order or reverse
			i=0;
			while(mc2_built_walks_direcs[0][i]!=0){
				i++;
			}
			if(mc2_built_walks_end[0][1]==walktoadd_start[1] && mc2_built_walks_end[0][2]==walktoadd_start[2]){
				//just append
				j=1; //Note: skip overlapping first edge (which will be a 1).
				while(walktoadd[j]!=2){	//also skip overlapping last edge (which will be a 2)
					mc2_built_walks_direcs[0][i]=walktoadd[j];
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
					mc2_built_walks_direcs[0][i] = mc2_reverse_direc(walktoadd[j]);
					i++;
					j--;
				}
			}
			//set actual startpoint of SAP (and equivalent endpoint)
			mc2_built_walks_start[0][0]=totalspan;
			for(i=0; i<=2; i++){
				mc2_built_walks_end[0][i]=mc2_built_walks_start[0][i];
				mc2_built_walks_end[0][i]=mc2_built_walks_start[0][i];
			}
		}
/*
printf("%d-th part of endhinge added:\n", walknum+1);
mc2_printbuiltwalks_all();
*/
	}






	///////////2ND POLY

	for(walknum=0; walknum<=mc2_Rend_num_walks2[secnum][nth_endhinge]-1; walknum++){ //for each walk in the right endhinge
		//fill walktoadd
		for(i=0; i<=2; i++){
			walktoadd_start[i]=mc2_Rend_start2[secnum][nth_endhinge][i][walknum];

			walktoadd_end[i]=mc2_Rend_end2[secnum][nth_endhinge][i][walknum];
		}
		length=0;
		for(i=0; i<=vM*vL; i++){
			walktoadd[i] = mc2_Rend_walks2[secnum][nth_endhinge][i][walknum];
			if(walktoadd[i]!=0){
				length++;
			}
		}

		if(mc2_num_built_walks2>1){
			//find a and b.
			for(built_walknum=0; built_walknum<=mc2_num_built_walks2-1; built_walknum++){
				if(mc2_built_walks_start2[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_start2[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}
				else if(mc2_built_walks_end2[built_walknum][1]==walktoadd_start[1] && mc2_built_walks_end2[built_walknum][2]==walktoadd_start[2]){
					a=built_walknum;
				}

				else if(mc2_built_walks_start2[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_start2[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
				else if(mc2_built_walks_end2[built_walknum][1]==walktoadd_end[1] && mc2_built_walks_end2[built_walknum][2]==walktoadd_end[2]){
					b=built_walknum;
				}
			}

			if(a<0 || b<0){
				printf("didn't find a or b (2ndpoly). Exitting\n");
				exit(1);
			}
			//fill newwalk with proper version of a. (start at opposite end of where walktoend starts)
			if(mc2_built_walks_end2[a][1]==walktoadd_start[1] && mc2_built_walks_end2[a][2]==walktoadd_start[2]){
				//correct direction. easy.
				for(i=0; i<=2; i++){
					newwalk_start[i]=mc2_built_walks_start2[a][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					newwalk[i] = mc2_built_walks_direcs2[a][i];
				}
			}
			else{
				//need to reverse a
				for(i=0; i<=2; i++){
					newwalk_start[i]=mc2_built_walks_end2[a][i];
				}
				i=0;
				while(mc2_built_walks_direcs2[a][i+1]!=0){
					i++;
				}
				j=0;
				while(i>=0){
					newwalk[j] = mc2_reverse_direc(mc2_built_walks_direcs2[a][i]);
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
			if(mc2_built_walks_start2[b][1]==walktoadd_end[1] && mc2_built_walks_start2[b][2]==walktoadd_end[2]){
				//correct direction. easy.
				//append b to newwalk. //i is already in the correct position.
				j=1;	//Note: skip overlapping first edge.
				while(mc2_built_walks_direcs2[b][j]!=0){
					newwalk[i] = mc2_built_walks_direcs2[b][j];
					i++;
					j++;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = mc2_built_walks_end2[b][i];
				}
			}
			else{
				//need to append reversed version of b
				j=0;
				while(mc2_built_walks_direcs2[b][j+1]!=0){
					j++;
				}
				j--; //Note: skip overlapping first edge.
				while(j>=0){
					newwalk[i] = mc2_reverse_direc(mc2_built_walks_direcs2[b][j]);
					i++;
					j--;
				}
				for(i=0; i<=2; i++){
					newwalk_end[i] = mc2_built_walks_start2[b][i];
				}
			}

			//replace a with newwalk.
			for(i=0; i<=2; i++){
				mc2_built_walks_start2[a][i] = newwalk_start[i];
				mc2_built_walks_end2[a][i] = newwalk_end[i];
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				mc2_built_walks_direcs2[a][i] = newwalk[i];
			}

			//move walk b to end of array
			while(b<mc2_num_built_walks2-1){
				for(i=0; i<=2; i++){
					mc2_built_walks_start2[b][i] = mc2_built_walks_start2[b+1][i];
					mc2_built_walks_end2[b][i] = mc2_built_walks_end2[b+1][i];
				}
				for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
					mc2_built_walks_direcs2[b][i] = mc2_built_walks_direcs2[b+1][i];
				}
				b++;
			}

			//delete b
			for(i=0; i<=2; i++){
				mc2_built_walks_start2[b][i] = -1;
				mc2_built_walks_end2[b][i] = -1;
			}
			for(i=0; i<=vM*vL*(totalspan+1)-1; i++){
				mc2_built_walks_direcs2[b][i] = 0;
			}
			mc2_num_built_walks2--;
		}
		else{	//close up SAP.
/*
printf("\ntime to close up the SAP:\n");
mc2_printbuiltwalks_all();
printf("walktoadd:");
i=0;
while(walktoadd[i]!=0){
	printf(" %d", walktoadd[i]);
	i++;
}
printf("\n");
*/			//either append walktoadd in same order or reverse
			i=0;
			while(mc2_built_walks_direcs2[0][i]!=0){
				i++;
			}
			if(mc2_built_walks_end2[0][1]==walktoadd_start[1] && mc2_built_walks_end2[0][2]==walktoadd_start[2]){
				//just append
				j=1; //Note: skip overlapping first edge (which will be a 1).
				while(walktoadd[j]!=2){	//also skip overlapping last edge (which will be a 2)
					mc2_built_walks_direcs2[0][i]=walktoadd[j];
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
					mc2_built_walks_direcs2[0][i] = mc2_reverse_direc(walktoadd[j]);
					i++;
					j--;
				}
			}
			//set actual startpoint of SAP (and equivalent endpoint)
			mc2_built_walks_start2[0][0]=totalspan;
			for(i=0; i<=2; i++){
				mc2_built_walks_end2[0][i]=mc2_built_walks_start2[0][i];
				mc2_built_walks_end2[0][i]=mc2_built_walks_start2[0][i];
			}
		}
/*
printf("%d-th part of endhinge added:\n", walknum+1);
mc2_printbuiltwalks_all();
*/
	}
}


int mc2_reverse_direc(int direc){
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
		printf("error in mc2_reverse_direc(). Exitting\n");
		exit(1);
	}
}


void mc2_printtofile(){
	if(mc2_filetotal>=maxpolys){
		//need to start new file
		printf("Finished filling file %lu. It contains %lu 2SAPs. Creating a new one.\n", mc2_filenum, mc2_filetotal);
		fprintf(mc2_fp, "-999\n");
		fclose(mc2_fp);
		mc2_filenum++;
		mkdir("data", 0775);
	mkdir("data/MonteCarlo", 0775);
	mkdir("data/MonteCarlo/2SAPs", 0775);
	mc_checked_snprintf(mc2_filename, sizeof(mc2_filename), "data/MonteCarlo/2SAPs/MC2SAPsL%dM%dspan%drun%dnum%lu.txt", L, M, totalspan, runnum, mc2_filenum);
		mc2_fp = fopen(mc2_filename, "w");	//create or overwrite the file "mc2_filename

		if(mc2_fp != NULL){
			fprintf(mc2_fp, "UofS\n");	//first line in file is always "UofS"
			run_metadata_write(mc2_filename, "mc_master", "samples_uofs", mode, L, M, totalspan, seednum, runnum, dom_evalue);
			printf("printed UofS in file '%s'\n", mc2_filename);
		}
		else{
			printf("file pointer is pointing to NULL\n");
			exit(1);
		}
		mc2_filetotal=0;
		
	}
	//record
	fprintf(mc2_fp, "%d %d %d\n", mc2_built_walks_start[0][0], mc2_built_walks_start[0][1], mc2_built_walks_start[0][2]);	//starting point
	int i=0;
	while(mc2_built_walks_direcs[0][i] != 0){
		fprintf(mc2_fp, "%d\n", mc2_built_walks_direcs[0][i]);
		i++;
	}
	fprintf(mc2_fp, "-111\n");

	fprintf(mc2_fp, "%d %d %d\n", mc2_built_walks_start2[0][0], mc2_built_walks_start2[0][1], mc2_built_walks_start2[0][2]);	//starting point
	i=0;
	while(mc2_built_walks_direcs2[0][i] != 0){
		fprintf(mc2_fp, "%d\n", mc2_built_walks_direcs2[0][i]);
		i++;
	}
	fprintf(mc2_fp, "-111\n");
	mc2_filetotal++;
}



#define MC2_NAME(name) mc2_##name
#include "mc_deps_2sap.c"		/*mat_ent **mc2_matalloc(int rowlow, int rowhigh, int collow, int colhigh)*/
#undef MC2_NAME
/* used for dynamically allocating memory for matrices with enteries of type mat_ent */



#undef num_outsections



	










static double mc2_max_eval_LRvec(double fugacity)
{
	McTransitionSpectralInput input;

	memset(&input, 0, sizeof(input));
	input.max_keynum = max_keynum;
	input.max_tspans = (unsigned long int)max_tspans;
	input.num_outsections = mc2_num_outsections;
	input.tspans_outsection = mc2_t_outsection;
	input.tspans_edges = mc2_MC_tspans_edges;
	input.tspans_nrr = mc2_t_nrr;
	input.left[0] = mc2_MC_L_Evector[0];
	input.left[1] = mc2_MC_L_Evector[1];
	input.right[0] = mc2_MC_R_Evector[0];
	input.right[1] = mc2_MC_R_Evector[1];
	input.force = fval;
	input.L = L;
	input.M = M;
	input.hamiltonian = 0;

	return mc_spectral_max_eval_2sap(&input, fugacity);
}

typedef struct Mc2CreatorSeenNode {
	char *key;
	struct Mc2CreatorSeenNode *next;
} Mc2CreatorSeenNode;

#define MC2_CREATOR_HASH_SIZE 262144
static Mc2CreatorSeenNode *mc2_creator_seen[MC2_CREATOR_HASH_SIZE];

static unsigned long int mc2_creator_hash_string(const char *text)
{
	unsigned long int h = 14695981039346656037ULL;
	while (*text) {
		h ^= (unsigned char)*text++;
		h *= 1099511628211ULL;
	}
	return h;
}

static char *mc2_creator_strdup(const char *text)
{
	size_t len = strlen(text) + 1;
	char *copy = (char *)mc_xcalloc(len, sizeof(*copy), "creator key copy");
	memcpy(copy, text, len);
	return copy;
}

static void mc2_creator_seen_clear(void)
{
	unsigned long int i;
	for (i = 0; i < MC2_CREATOR_HASH_SIZE; i++) {
		Mc2CreatorSeenNode *node = mc2_creator_seen[i];
		while (node) {
			Mc2CreatorSeenNode *next = node->next;
			free(node->key);
			free(node);
			node = next;
		}
		mc2_creator_seen[i] = NULL;
	}
}

static int mc2_creator_seen_insert(const char *key)
{
	unsigned long int bucket = mc2_creator_hash_string(key) % MC2_CREATOR_HASH_SIZE;
	Mc2CreatorSeenNode *node = mc2_creator_seen[bucket];
	while (node) {
		if (strcmp(node->key, key) == 0) return 0;
		node = node->next;
	}
	node = (Mc2CreatorSeenNode *)mc_xcalloc(1, sizeof(*node), "creator seen node");
	node->key = mc2_creator_strdup(key);
	node->next = mc2_creator_seen[bucket];
	mc2_creator_seen[bucket] = node;
	return 1;
}

static char *mc2_creator_poly_string(int which)
{
	size_t cap = (size_t)(vM * vL * (totalspan + 1) + 8) * 24;
	char *buf = (char *)mc_xcalloc(cap, sizeof(*buf), "creator polygon string");
	size_t used = 0;
	int i = 0;

	used += (size_t)snprintf(buf + used, cap - used, "%d,%d,%d:",
		which == 0 ? mc2_built_walks_start[0][0] : mc2_built_walks_start2[0][0],
		which == 0 ? mc2_built_walks_start[0][1] : mc2_built_walks_start2[0][1],
		which == 0 ? mc2_built_walks_start[0][2] : mc2_built_walks_start2[0][2]);
	while ((which == 0 ? mc2_built_walks_direcs[0][i] : mc2_built_walks_direcs2[0][i]) != 0) {
		used += (size_t)snprintf(buf + used, cap - used, "%d,",
			which == 0 ? mc2_built_walks_direcs[0][i] : mc2_built_walks_direcs2[0][i]);
		if (used >= cap) {
			fprintf(stderr, "Fatal: creator polygon string exceeded buffer\n");
			exit(EXIT_FAILURE);
		}
		i++;
	}
	return buf;
}

static char *mc2_creator_pair_key(void)
{
	char *a = mc2_creator_poly_string(0);
	char *b = mc2_creator_poly_string(1);
	const char *first = a;
	const char *second = b;
	size_t len;
	char *key;

	if (strcmp(a, b) > 0) {
		first = b;
		second = a;
	}
	len = strlen(first) + strlen(second) + 2;
	key = (char *)mc_xcalloc(len, sizeof(*key), "creator pair key");
	snprintf(key, len, "%s|%s", first, second);
	free(a);
	free(b);
	return key;
}

static void mc2_creator_reset_built_walks(void)
{
	int i, j;
	for (i = 0; i <= vM * vL / 2 - 1; i++) {
		for (j = 0; j <= 2; j++) {
			mc2_built_walks_start[i][j] = -1;
			mc2_built_walks_end[i][j] = -1;
			mc2_built_walks_start2[i][j] = -1;
			mc2_built_walks_end2[i][j] = -1;
		}
		for (j = 0; j <= vM * vL * (totalspan + 1) - 1; j++) {
			mc2_built_walks_direcs[i][j] = 0;
			mc2_built_walks_direcs2[i][j] = 0;
		}
	}
}

static void mc2_creator_load_left_endhinge(unsigned long int section, int nth_endhinge)
{
	int nth_walk, i;

	mc2_num_built_walks = mc2_Lend_num_walks[section][nth_endhinge];
	mc2_num_built_walks2 = mc2_Lend_num_walks2[section][nth_endhinge];
	for (nth_walk = 0; nth_walk <= mc2_num_built_walks - 1; nth_walk++) {
		for (i = 0; i <= 2; i++) {
			mc2_built_walks_start[nth_walk][i] = mc2_Lend_start[section][nth_endhinge][i][nth_walk];
			mc2_built_walks_end[nth_walk][i] = mc2_Lend_end[section][nth_endhinge][i][nth_walk];
		}
		for (i = 0; i <= vM * vL; i++) mc2_built_walks_direcs[nth_walk][i] = mc2_Lend_walks[section][nth_endhinge][i][nth_walk];
	}
	for (nth_walk = 0; nth_walk <= mc2_num_built_walks2 - 1; nth_walk++) {
		for (i = 0; i <= 2; i++) {
			mc2_built_walks_start2[nth_walk][i] = mc2_Lend_start2[section][nth_endhinge][i][nth_walk];
			mc2_built_walks_end2[nth_walk][i] = mc2_Lend_end2[section][nth_endhinge][i][nth_walk];
		}
		for (i = 0; i <= vM * vL; i++) mc2_built_walks_direcs2[nth_walk][i] = mc2_Lend_walks2[section][nth_endhinge][i][nth_walk];
	}
}

static void mc2_creator_print_one_poly(FILE *out, int which)
{
	int i = 0;

	fprintf(out, "%d %d %d\n",
		which == 0 ? mc2_built_walks_start[0][0] : mc2_built_walks_start2[0][0],
		which == 0 ? mc2_built_walks_start[0][1] : mc2_built_walks_start2[0][1],
		which == 0 ? mc2_built_walks_start[0][2] : mc2_built_walks_start2[0][2]);
	while ((which == 0 ? mc2_built_walks_direcs[0][i] : mc2_built_walks_direcs2[0][i]) != 0) {
		fprintf(out, "%d\n", which == 0 ? mc2_built_walks_direcs[0][i] : mc2_built_walks_direcs2[0][i]);
		i++;
	}
	fprintf(out, "-111\n");
}

static void mc2_creator_print_pair(FILE *out)
{
	char *a = mc2_creator_poly_string(0);
	char *b = mc2_creator_poly_string(1);
	int first = strcmp(a, b) <= 0 ? 0 : 1;
	free(a);
	free(b);
	mc2_creator_print_one_poly(out, first);
	mc2_creator_print_one_poly(out, first == 0 ? 1 : 0);
}

typedef struct Mc2CreatorContext {
	FILE *fp;
	unsigned long int total;
	unsigned long int file_total;
	unsigned long int file_num;
	char outdir[160];
	char prefix[64];
	char last_path[512];
	unsigned long int *sequence_sections;
	int *sequence_nths;
} Mc2CreatorContext;

static void mc2_creator_open_next_file(Mc2CreatorContext *ctx)
{
	if (ctx->fp) {
		fprintf(ctx->fp, "-999\n");
		fclose(ctx->fp);
		ctx->fp = NULL;
	}
	snprintf(ctx->last_path, sizeof(ctx->last_path), "%s/%sL%dM%dspan%dnum%lu.txt",
		ctx->outdir, ctx->prefix, L, M, totalspan, ctx->file_num);
	ctx->fp = fopen(ctx->last_path, "w");
	if (!ctx->fp) {
		fprintf(stderr, "Fatal: could not open '%s': %s\n", ctx->last_path, strerror(errno));
		exit(EXIT_FAILURE);
	}
	fprintf(ctx->fp, "UofS\n");
	ctx->file_total = 0;
}

static void mc2_creator_prepare_write(Mc2CreatorContext *ctx, const char *outdir, const char *prefix)
{
	snprintf(ctx->outdir, sizeof(ctx->outdir), "%s", outdir);
	snprintf(ctx->prefix, sizeof(ctx->prefix), "%s", prefix);
	ctx->file_num = 1;
	mc2_creator_open_next_file(ctx);
}

static void mc2_creator_finish_write(Mc2CreatorContext *ctx)
{
	if (ctx->fp) {
		fprintf(ctx->fp, "-999\n");
		fclose(ctx->fp);
		ctx->fp = NULL;
	}
}

static void mc2_creator_write_pair(Mc2CreatorContext *ctx)
{
	if (!ctx->fp) return;
	if (ctx->file_total >= (unsigned long int)maxpolys) {
		ctx->file_num++;
		mc2_creator_open_next_file(ctx);
	}
	mc2_creator_print_pair(ctx->fp);
	ctx->file_total++;
}

static void mc2_creator_build_leaf(
	Mc2CreatorContext *ctx,
	unsigned long int left_section,
	int left_endhinge,
	int right_endhinge,
	int sequence_len)
{
	int idx;
	unsigned long int terminal_section;
	char *key;

	mc2_creator_reset_built_walks();
	mc2_creator_load_left_endhinge(left_section, left_endhinge);
	for (idx = 0; idx < sequence_len; idx++) {
		mc2_add_to_built_walks(ctx->sequence_sections[idx], ctx->sequence_nths[idx]);
	}
	terminal_section = mc2_t_outsection[ctx->sequence_sections[sequence_len - 1]][ctx->sequence_nths[sequence_len - 1]];
	mc2_add_right_endhinge(terminal_section, right_endhinge);

	key = mc2_creator_pair_key();
	if (mc2_creator_seen_insert(key)) {
		mc2_creator_write_pair(ctx);
		ctx->total++;
	}
	free(key);
}

static void mc2_creator_enumerate_paths(
	Mc2CreatorContext *ctx,
	unsigned long int left_section,
	int left_endhinge,
	unsigned long int current_section,
	int depth,
	int target_depth)
{
	unsigned long int nth;
	if (depth == target_depth) {
		for (nth = 1; nth <= mc2_num_right_endhinges[current_section]; nth++) {
			mc2_creator_build_leaf(ctx, left_section, left_endhinge, (int)nth, target_depth);
		}
		return;
	}
	for (nth = 1; nth <= mc2_num_outsections[current_section]; nth++) {
		ctx->sequence_sections[depth] = current_section;
		ctx->sequence_nths[depth] = (int)nth;
		mc2_creator_enumerate_paths(ctx, left_section, left_endhinge, mc2_t_outsection[current_section][nth], depth + 1, target_depth);
	}
}

static unsigned long int mc2_creator_enumerate(Mc2CreatorContext *write_ctx)
{
	Mc2CreatorContext ctx;
	unsigned long int section;
	int left;

	memset(&ctx, 0, sizeof(ctx));
	if (write_ctx) ctx = *write_ctx;
	ctx.sequence_sections = (unsigned long int *)mc_xcalloc((size_t)totalspan, sizeof(*ctx.sequence_sections), "creator section sequence");
	ctx.sequence_nths = (int *)mc_xcalloc((size_t)totalspan, sizeof(*ctx.sequence_nths), "creator transition sequence");
	mc2_creator_seen_clear();

	for (section = 1; section <= max_keynum; section++) {
		for (left = 1; left <= mc2_num_left_endhinges[section]; left++) {
			mc2_creator_enumerate_paths(&ctx, section, left, section, 0, totalspan - 1);
		}
	}
	free(ctx.sequence_sections);
	free(ctx.sequence_nths);
	if (write_ctx) {
		write_ctx->fp = ctx.fp;
		write_ctx->total = ctx.total;
		write_ctx->file_total = ctx.file_total;
		write_ctx->file_num = ctx.file_num;
		memcpy(write_ctx->last_path, ctx.last_path, sizeof(write_ctx->last_path));
	}
	return ctx.total;
}

static int mc2_run_creator_all_after_build(void)
{
	char outdir[160];
	Mc2CreatorContext write_ctx;
	unsigned long int count;

	if (totalspan < 2) {
		fprintf(stderr, "Error: Span (-s) must be at least 2 (received %d).\n", totalspan);
		return EXIT_FAILURE;
	}

	count = mc2_creator_enumerate(NULL);
	if (!mc2_creator_force && count > mc2_creator_limit) {
		fprintf(stderr,
			"Refusing to write %lu unordered 2SAPs; limit is %lu. Re-run with -C <limit> or -f to force.\n",
			count, mc2_creator_limit);
		mc2_creator_seen_clear();
		return EXIT_FAILURE;
	}
	mc2_creator_seen_clear();

	mkdir("data", 0775);
	mkdir("data/CreatorAll", 0775);
	snprintf(outdir, sizeof(outdir), "data/CreatorAll/All_2SAPs");
	mkdir(outdir, 0775);
	memset(&write_ctx, 0, sizeof(write_ctx));
	mc2_creator_prepare_write(&write_ctx, outdir, "All2SAPs");
	(void)mc2_creator_enumerate(&write_ctx);
	mc2_creator_finish_write(&write_ctx);
	mc2_creator_seen_clear();
	{
		char metadata_target[256];
		mc_checked_snprintf(metadata_target, sizeof(metadata_target), "%s/All2SAPsL%dM%dspan%d.summary", outdir, L, M, totalspan);
		run_metadata_write(metadata_target, "creator_all", "creator_all_summary", mode, L, M, totalspan, seednum, runnum, 0.0);
	}

	printf("Generated %lu unordered 2SAPs of exact span %d in %dx%d tube.\n", count, totalspan, L, M);
	printf("Wrote %lu file(s) under %s, with at most %d entries per file.\n", write_ctx.file_num, outdir, maxpolys);
	return EXIT_SUCCESS;
}

int run_integrated_2sap_creator_all(int argc, char *argv[], unsigned long int limit, int force)
{
	mc2_creator_all_mode = 1;
	mc2_creator_limit = limit;
	mc2_creator_force = force;
	return run_integrated_2sap_sampler(argc, argv);
}

#undef MC2_CREATOR_HASH_SIZE
