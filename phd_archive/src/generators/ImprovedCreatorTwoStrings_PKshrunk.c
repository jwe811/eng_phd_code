/*	This program will create all polygons in a LxM tube, with chosen span, that contain a "2-string pattern". This means their
start and end sections only have two edges, and any middle sections have >2 edges. Polygon must go the full span. ALSO-polys must have a min number of edges. 	*/
//24 since we're going to shrink and might go outside the tube ALSO, this program will try to shrink the polygon before printing it using BFACF moves (bfacf_shrink.c) 

//Changes:	
//			Create file called "tracking_v3.txt". Will output stuff as the program runs.

//changes:	Only interested in patterns, so only generate each pattern with one end-hinge, instead of all of them.
//			filenames tracking_v3.txt and _v3 on poly file names



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
#define	L 3				/* number of horizontal edges							*/
#define	M 1				/* number of vertical edges							*/
#define	totalspan 3		/* span of the polygon that will be generated (must be >=3)	*/

#define 	maxpolys 1000000	/* maximum number of polygons in each .txt file	*/





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





/***************************************************************************/
/****************** Structures used in this program ************************/
/***************************************************************************/


struct hinge_span {	/* This data structure is used to store information about two-spans */

	/*unsigned long int 	inorder */		/* no longer needed as it is implicitely defined */
	unsigned long int 	outorder	:26;	/* the section number of the righthand section   */	
	unsigned long int 	edgecount	:6;	/* the number of edges in the hinge		 */
	struct hinge_span 	*nexthinge;		/* pointer to another hinge_span data structure  */
	int				walk[vM*vL*8];	/* this will hold the random walks that the 2span creates (walks are seperated by -222 or -333)	*/
	int				startpoint[2];	//contains the (y,z) coordinates of where the walk is starting (we know x=0)


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


struct endhinge {			/* an endhinge is a linked list of edges			*/
	int endwalk[vM*vL*8];	/* this will hold the random walks in the endhinge (walks are seperated by -222 or -333)	*/
	int side;			//0 if a left endhinge, 1 if a right endhinge
	int start[2];		//contains the (y,z) coordinates that the walk starts
	struct endhinge* nextendhinge;	/* pointer to the next endhinge structure in the linked list	*/
};

struct endhinge* firstendhinge[vec_length];	/* Array of pointers to endhinge structures, where each	*/
										/* pointer in the array points to the first element in	*/
										/* the linked list of endhinge structures		*/
struct endhinge* currentendhinge[vec_length];	/* Array of dummy pointers that are used to traverse	*/
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

int				curwalk[vM*vL*4];	/* this will hold the current walks, to be recorded when the 2span is complete. (-222 or -333 separates walks)*/
int				startingy=0, startingz=0;
int				total_end_hinges=0;		/*the total number of possible endhinges	*/
unsigned long int	curindex=0;			//index of which endhinge is being created for ordertemplate (initially zero)
unsigned long int	edgesin_temp = 0;	//will count the number of edges in temp[]
unsigned short int	occupied[vM][vL];	//this will indicate which verticies are occupied (i,j)
int				leftside;				//leftside=1 if testing endhinges on the left; otherwise (testing endhinges on right) leftside=0.

int				pwalksize = vM*vL*(totalspan+1)*4;
int				filetotal=0;
char 			filename[100];			//this will hold the filename of the text file that we want to write to
FILE*			fp;	//file pointer. will point to the file that is to be written to.
int				filenum=1;
int				passedM=M;
int				passedL=L;
int curious=0;
int grand=0;
int grandtotal=0;

int endtemplate[vM][vL];		//contains the section associated with the end hinge
int endtemplate2[vM][vL];
int endstatus[vM][vL];		//contains info on which vertices in the end hinge are filled. 1 if filled, 0 if empty.
int EndOrdNum[2];		//contains info on # of edges in left section and endhinge
int num_right_endhinges=0;

int minedges=24;
unsigned long int num_millions=0;
unsigned long int num_uptomill=0;
int smallcounter=0;

FILE* fptwo;
char trackingname[50];


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
unsigned short int	LFlag_endhinge( int (*)[2]);
unsigned short int	RFlag( int (*)[3]);
unsigned short int	RFlag_endhinge( int (*)[2]);
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


unsigned short int testmatch(int i, int j);

void connect(unsigned long int small, unsigned long int big, unsigned long int cindex);
void goup(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int uindex);
void godown(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int dindex);
void goleft(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int lindex);
void goright(unsigned long int prevx, unsigned long int prevy, unsigned long int endx, unsigned long int endy, unsigned long int firstconnum, unsigned long int secondconnum, unsigned long int rindex);


void printcurwalk(void);
void addtowalk(int printablewalk[], int nextwalk[]);
void recursivegenerate(int curspan, int printablewalk[], struct hinge_span* tspantoadd, int firstsec, int x, int y);
int dupcheck(int x, int y, int dup, int firstsection);


void enterendhinge(int i, int j);
void leaveendhinge(int i, int j);
void endrowedges(int i, int j);
void endcoledges(int i, int j);
unsigned long int num_section_endhinge(void);
unsigned long int num_section_endhinge2(void);
void recordendtemplate(void);

int test2secs(int sectionnum);
int bfacf_shrink(int x, int y, int z, int polygon[], int length);

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
		currentendhinge[i] = newendhinge();
		firstendhinge[i] = currentendhinge[i];

	}
	printf("space created for hinge structures and endhinges\n\n");



	//initialize curwalk[] to be all zeros
	for(i=0; i<=vM*vL*4-1; i++){
		curwalk[i]=0;
	}

	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			printf("entering hinge at (y,z) = (%d,%d)\n", j, i);
			startingy=j;
			startingz=i;
			enterhinge(makeheaderij, makeheaderij, i, j, side, &ordNum);
//			printf("finished at i=%d, j=%d\n", i, j);
		}
	}


	printf("There were %d valid 2-spans\n\n", valid_2_spans);


	//generate the endhinges
	//we are going to create "right endhinges" (entering from the left, can't leave to the right)
	//When we record these right endhinges, then we will make appropriate changes to also record it's "mirror" left endhinge

	//initialize curwalk[] again to be all zeros
	for(i=0; i<=vM*vL*4-1; i++){
		curwalk[i]=0;
	}

	EndOrdNum[0]=1;	//first edge in left section will get the label 1
	EndOrdNum[1]=0;	//no edges currently in the endhinge

	//initialize endtemplate[][] to zeros
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			endtemplate[i][j]=0;
		}
	}

	
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("creating endhinges that enter from (y,z)=(%d,%d)\n", j, i);
			startingy = j;
			startingz = i;
			enterendhinge(i,j);
		}
	}

	printf("There were %d right endhinges ( and the same number of left endhinges)\n\n", num_right_endhinges);
/*
	j=0;
	for(i=1; i<=max_sections; i++){
		currentendhinge[i] = firstendhinge[i];
		while((*currentendhinge[i]).nextendhinge != NULL){
			currentendhinge[i] = (*currentendhinge[i]).nextendhinge;
			printf("We have an endhinge on side %d that connects to section %d\n", currentendhinge[i]->side, i);
			j++;
			
			
		}
	}
	printf("total of %d endhinges\n", j);
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//CREATE POLYGONS HERE
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	sprintf(filename, "ImpCreatorTwoStrings_PKshrunkL%dM%dspan%dnum%d.txt", L, M, totalspan, filenum);
	fp = fopen(filename, "w");	//create or overwrite the file "filename"

	if(fp != NULL){
		fprintf(fp, "UofS\n");	//first line in file is always "UofS"
		printf("printed UofS in file '%s'\n", filename);
	}
	else{
		printf("file pointer is pointing to NULL\n");
		exit(1);
	}

	sprintf(trackingname, "tracking_ImpCr2stringsPKshrL%dM%dspan%d.txt", L, M, totalspan);
	fptwo = fopen(trackingname, "w");
	if(fptwo != NULL){
		fprintf(fptwo, "Starting to create polys\n");
		printf("%s is created\n", trackingname);
	}
	else{
		printf("file pointer fptwo is pointing to NULL\n");
		exit(1);
	}
	fclose(fptwo);



	int printablewalk[pwalksize];	//this will store the walk that is to be printed
	for(i=0; i<=pwalksize-1; i++){	//initialize printable walk to be all zeros
		printablewalk[i]=0;
	}


	int x, y, z=0;							//this will be the coordinates of the starting point
	int cury, curz, curendhinge;				//going to use these variables to aid in adding the endhinge to printablewalk.
	int newdirec, prevdirec, direcsadded;
	int temp, temp2;
	int counter;
	int dummy;
	int firstsection;
	int nextsection;
	int memorywalk[pwalksize];	//this will remember what printablewalk is before recursive generate is ran
	int dup;
	int index;
	int temp1;
	int foundlefthinge;

	//loop through all 2-spans(next 2 lines), then check if it has a starting (left) endhinge
	for(firstsection=1; firstsection<=max_sections; firstsection++){
		if(test2secs(firstsection)==1){	//ONLY WANT POLYGONS WITH 2 SECTIONS IN THE FIRST SECTION
			fptwo = fopen(trackingname, "a");
			fprintf(fptwo, "Looking at polys with first section=%d\n", firstsection);
			fclose(fptwo);
	//		printf("firstsection=%d\n", firstsection);
			current_hinge_span[firstsection] = first_hinge_span[firstsection];	//get current pointing to first
			while((*current_hinge_span[firstsection]).nexthinge != NULL){
				current_hinge_span[firstsection] = (*current_hinge_span[firstsection]).nexthinge;	//current is now pointing at a possible starting 2-span
				//		printf("secondsection=%d\n", current_hinge_span[firstsection]->outorder);
				//Find a left endhinge
				foundlefthinge=0;
				currentendhinge[firstsection] = firstendhinge[firstsection];	//get current pointing to first
				while((*currentendhinge[firstsection]).nextendhinge != NULL && foundlefthinge==0){	//run this only if we haven't had a left hinge yet.
					currentendhinge[firstsection] = (*currentendhinge[firstsection]).nextendhinge;
					if((*currentendhinge[firstsection]).side==0){
						foundlefthinge=1;
	//					printf("We have a valid start 2-span and start endhinge\n");
	//					printf("firstsection=%d\n", firstsection);
						y = (*current_hinge_span[firstsection]).startpoint[0];	//record y-coord of startpoint
						z = (*current_hinge_span[firstsection]).startpoint[1];	//record z-coord of startpoint
						//if starting point isn't (0,0), need to check if it's a duplicate.
						//THIS IS THE DUPLICATE TEST
						dup=0;
						if(y>0 || z>0){
	//						printf("starting dup test. startpoint=(%d,%d). \n", y, z);
							dup = 0; //assume initially it is NOT a duplicate
							dup = dupcheck(y, z, dup, firstsection);
							if(dup==1){
	//							printf("found a duplicate. not printing it.\n");
							}
						}

						if(dup==0){

							//put walk[] into printablewalk[]
		//					printf("putting first 2-spans walk into printable walk: ");
							i=0;
							while((*current_hinge_span[firstsection]).walk[i] != 0){
								printablewalk[i] = (*current_hinge_span[firstsection]).walk[i];
		//						printf("%d ", printablewalk[i]);
								i++;
							}
							while(i<=pwalksize-1){	//fill rest with zeros
								printablewalk[i] =0;
								i++;
							}
		//					printf("\n");


		///////////////////////////////////////////////////////////////////////////////////////////////////////////
							//Now add in left endhinge walk
		/*					printf("adding left endhinge walk: \n");
							i=0;
							while((*currentendhinge[firstsection]).endwalk[i] != 0){;
								printf("%d ", (*currentendhinge[firstsection]).endwalk[i]);
								i++;
							}
							printf("\n");
							printf("printable walk is this: ");
							i=0;
							while(printablewalk[i]!=0){
								printf("%d ", printablewalk[i]);
								i++;
							}
							printf("\n");
		*/
							i=0;
							j=0;
							while(printablewalk[i]!=0){
								while(printablewalk[i]!=-333){
									i++;
								}
								//insert part of left endhinge walk here
		//						printf("inserting into spot %d\n", i);
								while((*currentendhinge[firstsection]).endwalk[j]!=-333){
									if((*currentendhinge[firstsection]).endwalk[j] == 3 ||(*currentendhinge[firstsection]).endwalk[j] == 4 ||(*currentendhinge[firstsection]).endwalk[j] == 5 ||(*currentendhinge[firstsection]).endwalk[j] == 6){
										//insert (*currentendhinge).endwalk[j] into printable walk[i]
										index=i;
		//								printf("RUNNING. i=%d, index=%d, j=%d\n", i, index, j);
										//shift printablewalk one over
										temp1=printablewalk[index];
										while(printablewalk[index]!=0){
											temp2=printablewalk[index+1];
											printablewalk[index+1]=temp1;
											temp1=temp2;
											index++;
										}
										//add nextwalk component
										printablewalk[i]=(*currentendhinge[firstsection]).endwalk[j];
										j++;
										i++;
									}
									else{
										j++;
									}
								}
		//						printf("Done inserting part of the endhinge. Can now remove the negative parts of printablewalk. i=%d, j=%d\n", i, j);
								while(printablewalk[i]<0){
		//							printf("deleting printablewalk[%d]=%d\n", i, printablewalk[i]);
									index=i;
									while(printablewalk[index]!=0){
										printablewalk[index] = printablewalk[index+1];
										index++;
									}
								}
								j++;
								i++;

							}
		/*					printf("printable walk is now: ");
							i=0;
							while(printablewalk[i]!=0){
								printf("%d ", printablewalk[i]);
								i++;
							}
							printf("\n");
		*/
		//////////////////////////////////////////////////////////////////////////////
							//finished adding in left endhinge
							//now add intermediate 2-spans

							//copy printablewalk into memorywalk
							for(i=0; i<=pwalksize-1; i++){
								memorywalk[i] = printablewalk[i];
							}

							nextsection=current_hinge_span[firstsection]->outorder;
							if( test2secs(nextsection)==0){	//second section can't be a 2string
								struct hinge_span* nexttspan;
								nexttspan = first_hinge_span[nextsection];	//point nextsection's current at first
								while(nexttspan->nexthinge != NULL){		//loop through all of nextsection's outorders.
									nexttspan = nexttspan->nexthinge;
			//						printf("\n\nfirst section=%d, second section=%d, third section =%d\n", firstsection, nextsection, nexttspan->outorder);
			//						printf("thirdsection=%d\n", nexttspan->outorder);
									recursivegenerate(2, printablewalk, nexttspan, nextsection, y, z);	//all polygons now generated with this starting 2-span.
									//reset printablewalk to memorywalk
									for(i=0; i<=pwalksize-1; i++){
										printablewalk[i]=memorywalk[i];
									}
								}
							}
						}
					}
				}
			}
		}
	}







//	exit(1);
	fprintf(fp, "-999\n");

	printf("L=%d, M=%d, span=%d, filetotal==%d, filenum=%d\n", L, M, totalspan, filetotal, filenum);
	printf("Only recording polygons with >=%d edges AFTER BEING SHRUNK\n", minedges);
	if(filenum==1){
		printf("Total Number of two-string polys = %d = %f = %e\n", filetotal + smallcounter, (float)filetotal + smallcounter, (float)filetotal + smallcounter);
		printf("Number of two-string Polys with >=%d edges AFTER BEING SHRUNK = %d = %f = %e!\n", minedges, filetotal, (float)filetotal, (float)filetotal);
	
	}
	else{
		printf("Total Number of two-string polys = %d = %f = %e\n", filetotal + ((filenum-1)*(maxpolys)) + smallcounter, (float)filetotal + ((filenum-1)*(maxpolys)) + smallcounter, (float)filetotal + ((filenum-1)*(maxpolys)) + smallcounter);
		printf("Number of two-string Polys with >=%d edges AFTER BEING SHRUNK = %d = %f = %e!\n", minedges, filetotal + ((filenum-1)*(maxpolys)), (float)filetotal + ((filenum-1)*(maxpolys)), (float)filetotal + ((filenum-1)*(maxpolys)));
	}

	printf("grand=%d\n", grand);

	printf("\nDone generating all 2-string polygons with >= %d edges AFTER BEING SHRUNK in %f seconds.\n\n", minedges, (double)clock()/(double)CLOCKS_PER_SEC);
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

	//add direction 1 or 2 to curwalk
	int index=0;
	while(curwalk[index]!=0){
		index++;
	}
	//if side==0, then we are entering from the left, so we are in direction posx (1)
	//else, we are on side==1, and entering from the right (negx, 2)
	if(side==0){
		curwalk[index]=1;
	}
	else{
		curwalk[index]=2;
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

	//add direction 1 or 2 to curwalk
	int index=0;
	while(curwalk[index]!=0){
		index++;
	}
	//if side==0, then we are exiting to the left, so we are in direction negz (6) (also add -333 for exiting left)
	//else, we are exiting to side==1 (right) (posz, 5) (also add -222 for exiting right)
	if(side==0){
		curwalk[index]=2;
		curwalk[index+1]=-333;
	}
	else{
		curwalk[index]=1;
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

				//add horizontal steps (-3's or -4's)
				if(jj>j){		//if reenterj>originalj, add -3's
					for(loopj=j+1; loopj<=jj; loopj++){
						curwalk[index]=-3;
						index++;
						totadded++;
					}

				}
				else if(jj<j){	//add -4's
					for(loopj=jj+1; loopj<=j; loopj++){
						curwalk[index]=-4;
						index++;
						totadded++;
					}
				}

				//add vertical steps (-5's or -6's)
				if(ii>i){		//if reenteri>originali, add -5's
					for(loopi=i+1; loopi<=ii; loopi++){
						curwalk[index]=-5;
						index++;
						totadded++;
					}
				}
				else if(ii<i){	//add -6's
					for(loopi=ii+1; loopi<=i; loopi++){
						curwalk[index]=-6;
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
	int             hmove;	/* -1 corresponds to moving "negy"  */
				/* +1 corresponds to moving "posy" */
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

			//add direction 3 or 4 to curwalk
			int index=0;
			while(curwalk[index]!=0){
				index++;
			}
			//if hmove==-1, add direction negy (4)
			//elseif hmove=1, add direction posy (3)
			if(hmove==-1){
				curwalk[index]=4;
			}
			else if(hmove==1){
				curwalk[index]=3;
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
	int             vmove;	/* -1 corresponds to moving "posz"   */
				/* +1 corresponds to moving "negz" */
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

			//add direction 5 or 6 to curwalk
			int index=0;
			while(curwalk[index]!=0){
				index++;
			}
			//if vmove==-1, add direction negy (6)
			//elseif vmove=1, add direction posy (5)
			if(vmove==-1){
				curwalk[index]=6;
			}
			else if(vmove==1){
				curwalk[index]=5;
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
	(*current_hinge_span[inNum]).startpoint[0]=startingy;
	(*current_hinge_span[inNum]).startpoint[1]=startingz;
//	printf("(%d, %d, 0), ", (*current_hinge_span[inNum]).startpoint[0], (*current_hinge_span[inNum]).startpoint[1]);

	//copy curwalk into the hinge_span structure (and print the walk)
	i=0;
	while(curwalk[i]!=0){
		(*current_hinge_span[inNum]).walk[i] = curwalk[i];
//		printf("%d, ", curwalk[i]);
		i++;
	}

//	printcurwalk();
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
#include "../topology/LFlag_endhinge.c"

#include "../topology/RFlag.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the right and returns 0 otherwise */
/* it also uses the global variable ordertemplate */
#include "../topology/RFlag_endhinge.c"

#if vM*vL<7
#include "../sections/Num_section_6V.c"  /*unsigned long int num_section(int side)*/
/* This function takes the arguement (int side) and assigns a unique number */
/* to the section on side = side */
#include "../sections/Num_section_6V_endhinge.c"
#elif vM*vL<9
#include "../sections/Num_section_8V.c"
#include "../sections/Num_section_8V_endhinge.c"
#elif vM*vL<11
#include "../sections/Num_section_10V.c"
#include "../sections/Num_section_10V_endhinge.c"
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


struct endhinge *
newendhinge(void)		/* returns a pointer to a endhinge structure */
/* this is a memory allocation function */
{
	struct endhinge *nextnewendhinge;

	nextnewendhinge = ((struct endhinge *) calloc(1, sizeof(struct endhinge)));
	if (nextnewendhinge == NULL) {
		fprintf(stderr, "unable to allocate memory");
		exit(1);
	}
	int i;
	for(i=0; i<=vM*vL*8-1; i++){
		(*nextnewendhinge).endwalk[i]=0;		//initialize to zero
	}
	return nextnewendhinge;
}

/***************************************************************************/

void printcurwalk(void){
	printf("curwalk: ");
	int i=0;
	while(curwalk[i]!=0){
		printf("%d, ", curwalk[i]);
		i++;
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





/*******************************************************************************************************/
//endhinge functions

void enterendhinge(int i, int j) {
/*	int blahi, blahj;
	printf("entering hinge at (y,z) (%d,%d). endtemplate looks like this:\n", j, i);
	for(blahi=0; blahi<=M; blahi++){
		for(blahj=0; blahj<=L; blahj++){
			printf("endtemplate[%d][%d] = %d\n", blahi, blahj, endtemplate[blahi][blahj]);
		}
	}

	printf("curwalk is this: ");
	blahi=0;
	while(curwalk[blahi] != 0){
		printf("%d ", curwalk[blahi]);
		blahi++;
	}
	printf("\n");
*/

	//entering hinge at (y,z) = (j,i). Already know it's empty.
	//add direction 1 curwalk
	int index=0;
	while(curwalk[index]!=0){
		index++;
	}
	curwalk[index]=1;

	endtemplate[i][j] = EndOrdNum[0];
		/* section edge through which the SAW enters the hinge is given a number */
		/* this number represents the order in which it occurs in the section    */
	EndOrdNum[0]++;
		/* the edge number through which the hinge is entered */
		/* is incremented so that the next edge gets a higher number */
	endstatus[i][j] = 1;	//vertex is now filled.

	/* the SAW can now do one of three things */
	/* it may explore edges in the horizontal direction => rowedges (...) is called */
	/* it may explore edges in the vertical direction => coleges (...) is called    */
	/* in each case the currrent state of the hinge is passed via the parameters    */

	endrowedges(i, j);
		/* horizontal edges are explored */
	endcoledges(i, j);
		/* vertical edges are explored */

////////////////////

	EndOrdNum[0]--;
		/* the edge number is decremented since, control is being passed back to  */
		/* the function that called enterhinge. That is to say, the program will  */
		/* now explore alternate ways of entering the hinge (i.e. from different  */
		/* locations). In essense the intial incremention above is nullified so	  */
		/* when enterhinge is called again the entering edge is numbered properly */
		
	endtemplate[i][j] = 0;
		/* the edge that was recorded must now be erased as this edge may not exist */
		/* as part of the SAW when alternate possibilities are explored. In the case */
		/* that it does not exist it needs to be numbered as "0" (zero).             */
	endstatus[i][j] = 0;

	//delete last entry in curwalk
	index=0;
	while(curwalk[index+1] != 0){
		index++;
	}
	curwalk[index]=0;

	return;

}


void endrowedges(int i, int j) {
	//at the point (y,z) = (j,i), want to explore moves in the y-direction
	if(j>0){ 	//if there's space to move in -y direction
		if(endstatus[i][j-1]==0){	//and if the space is empty
			//add direction 4
			int index=0;
			while(curwalk[index]!=0){
				index++;
			}
			curwalk[index]=4;

			EndOrdNum[1]++;	/* count the edge in the endhinge */
			endstatus[i][j-1] = 1;	//(i,j-1) is now occupied

			//now explore options from there:
			leaveendhinge(i, j-1);
			endrowedges(i, j-1);
			endcoledges(i, j-1);

			EndOrdNum[1]--;
				/* decrement the number of edges in the endhinge */
			endstatus[i][j-1] = 0;	//now empty

			//delete last entry in curwalk
			index=0;
			while(curwalk[index+1] != 0){
				index++;
			}
			curwalk[index]=0;
		}
	}

	if(j<L){ 	//if there's space to move in +y direction
		if(endstatus[i][j+1]==0){	//and if the space is empty
			//add direction 3
			int index=0;
			while(curwalk[index]!=0){
				index++;
			}
			curwalk[index]=3;

			EndOrdNum[1]++;	/* count the edge in the endhinge */
			endstatus[i][j+1] = 1;	//(i,j-1) is now occupied

			//now explore options from there:
			leaveendhinge(i, j+1);
			endrowedges(i, j+1);
			endcoledges(i, j+1);

			EndOrdNum[1]--;
				/* decrement the number of edges in the endhinge */
			endstatus[i][j+1] = 0;	//now empty

			//delete last entry in curwalk
			index=0;
			while(curwalk[index+1] != 0){
				index++;
			}
			curwalk[index]=0;
		}
	}
	return;
}


void endcoledges(int i, int j) {
	//at the point (y,z) = (j,i), want to explore moves in the z-direction
	if(i>0){ 	//if there's space to move in -z direction
		if(endstatus[i-1][j]==0){	//and if the space is empty
			//add direction 6
			int index=0;
			while(curwalk[index]!=0){
				index++;
			}
			curwalk[index]=6;

			EndOrdNum[1]++;	/* count the edge in the endhinge */
			endstatus[i-1][j] = 1;	//(i-1,j) is now occupied

			//now explore options from there:
			leaveendhinge(i-1, j);
			endrowedges(i-1, j);
			endcoledges(i-1, j);

			EndOrdNum[1]--;
				/* decrement the number of edges in the endhinge */
			endstatus[i-1][j] = 0;	//now empty

			//delete last entry in curwalk
			index=0;
			while(curwalk[index+1] != 0){
				index++;
			}
			curwalk[index]=0;
		}
	}

	if(i<M){ 	//if there's space to move in +z direction
		if(endstatus[i+1][j]==0){	//and if the space is empty
			//add direction 5
			int index=0;
			while(curwalk[index]!=0){
				index++;
			}
			curwalk[index]=5;

			EndOrdNum[1]++;	/* count the edge in the endhinge */
			endstatus[i+1][j] = 1;	//(i,j-1) is now occupied

			//now explore options from there:
			leaveendhinge(i+1, j);
			endrowedges(i+1, j);
			endcoledges(i+1, j);

			EndOrdNum[1]--;
				/* decrement the number of edges in the endhinge */
			endstatus[i+1][j] = 0;	//now empty

			//delete last entry in curwalk
			index=0;
			while(curwalk[index+1] != 0){
				index++;
			}
			curwalk[index]=0;
		}
	}
	return;
}



void leaveendhinge(int i, int j)
/* called by: rowedges,coledges,enterhinge */
{

	//add direction 2 to curwalk
	int index=0;
	while(curwalk[index]!=0){
		index++;
	}
	curwalk[index]=2;
	curwalk[index+1]=-333;

	int       ii;			/*vertex number in vertical direction*/
	int		jj;			/*vertex number in horizontal direction */
	int		shift = 0;	/*keeps track of how far the headernode is shifted */
						/*from the 0,0 position going left to right first */
						/*then down as each row of vertices is exhausted  */
	int loopi, loopj, looper;
	int fullhinge=1;
	int totadded;

	endtemplate[i][j] = EndOrdNum[0];	/* record the section
								 * edge (i.e. the
								 * leaving edge) */
	EndOrdNum[0]++;
		/* the edge number is incremented so that the next edge gets a higher number */

	if(RFlag_endhinge(&EndOrdNum)==1){
		//record endtemplate
		recordendtemplate();
		num_right_endhinges++;
	}
/*	else{
		printf("Rflag_endhinge failed on this endhinge:\n");
		for(loopi=0; loopi<=M; loopi++){
			for(loopj=0; loopj<=L; loopj++){
				printf("endtemplate[%d][%d] = %d\n", loopi, loopj, endtemplate[loopi][loopj]);
			}
		}

	}
*/
	//try to re-enter at an empty space
	for(ii=0; ii<=M; ii++){
		for(jj=0; jj<=L; jj++){
			if(endstatus[ii][jj] == 0){
				//need to add neg numbers outside of hinge.
				//left endhinge at j,i, and reentering at jj,ii

				index=0;
				while(curwalk[index]!=0){
					index++;
				}

				totadded=0;	//will keep track of how many steps outside of the hinge are added, so we can delete them later

				//add y-steps (-3's or -4's)
				if(jj>j){		//if reenterj>originalj, add -3's
					for(loopj=j+1; loopj<=jj; loopj++){
						curwalk[index]=-3;
						index++;
						totadded++;
					}

				}
				else if(jj<j){	//add -4's
					for(loopj=jj+1; loopj<=j; loopj++){
						curwalk[index]=-4;
						index++;
						totadded++;
					}
				}

				//add z-steps (-5's or -6's)
				if(ii>i){		//if reenteri>originali, add -5's
					for(loopi=i+1; loopi<=ii; loopi++){
						curwalk[index]=-5;
						index++;
						totadded++;
					}
				}
				else if(ii<i){	//add -6's
					for(loopi=ii+1; loopi<=i; loopi++){
						curwalk[index]=-6;
						index++;
						totadded++;
					}
				}

				enterendhinge(ii,jj);

				//delete the last "totadded" from curwalk
				index=0;
				while(curwalk[index]!=0){
					index++;
				}
				for(looper=1; looper<=totadded; looper++){
					curwalk[index-looper]=0;
				}


			}
		}
	}

	EndOrdNum[0]--;
	endtemplate[i][j]=0;
	//delete last 2 entries in curwalk
	index=0;
	while(curwalk[index+2] != 0){
		index++;
	}
	curwalk[index]=0;
	curwalk[index+1]=0;


	return;
}





void recordendtemplate(void)
/* This function records the pertinent information about an endtemplate */
{
	int i, j;
//	printf("recording an endhinge\n");
	//get section number
	int secnum = num_section_endhinge();
//	printf("This endhinge connects to section %d\n", secnum);
/*	printf("endtemplate[][] looks like this:\n");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("endtemplate[%d][%d] = %d\n", i, j, endtemplate[i][j]);
		}
	}
*/
	//create newendhinge
	(*currentendhinge[secnum]).nextendhinge = newendhinge(); 
	currentendhinge[secnum] = (*currentendhinge[secnum]).nextendhinge;
	(*currentendhinge[secnum]).side = 1;
	(*currentendhinge[secnum]).start[0] = startingy;
	(*currentendhinge[secnum]).start[1] = startingz;
	//record walk
	i=0;
	while(curwalk[i]!=0){
		(*currentendhinge[secnum]).endwalk[i] = curwalk[i];
		i++;
	}

/*	printf("This is an endhinge with side=%d\n", (*currentendhinge[secnum]).side);
	printf("This endhinge starts at (y,z) = (%d, %d)\n", (*currentendhinge[secnum]).start[0], (*currentendhinge[secnum]).start[1]);
	printf("Here is the walk: ");
	i=0;
	while((*currentendhinge[secnum]).endwalk[i] != 0){
		printf("%d ", (*currentendhinge[secnum]).endwalk[i]);
		i++;
	}
	printf("\n");
*/

	//now record corresponding left endhinge.
	//get new section number
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

	int secnum2 = num_section_endhinge2();
	//create newendhinge
	(*currentendhinge[secnum2]).nextendhinge = newendhinge(); 
	currentendhinge[secnum2] = (*currentendhinge[secnum2]).nextendhinge;
	(*currentendhinge[secnum2]).side = 0;
	(*currentendhinge[secnum2]).start[0] = startingy;
	(*currentendhinge[secnum2]).start[1] = startingz;
	//record walk
	i=0;
	while(curwalk[i]!=0){
		if(curwalk[i]==1){
			(*currentendhinge[secnum2]).endwalk[i] = 2;
		}
		else if (curwalk[i]==2){
			(*currentendhinge[secnum2]).endwalk[i] = 1;
		}
		else{
			(*currentendhinge[secnum2]).endwalk[i] = curwalk[i];
		}
		i++;
	}

/*
	printf("This endhinge connects to section %d\n", secnum2);
	printf("endtemplate[][] looks like this:\n");
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			printf("endtemplate2[%d][%d] = %d\n", i, j, endtemplate2[i][j]);
		}
	}


	printf("This is an endhinge with side=%d\n", (*currentendhinge[secnum2]).side);
	printf("This endhinge starts at (y,z) = (%d, %d)\n", (*currentendhinge[secnum2]).start[0], (*currentendhinge[secnum2]).start[1]);
	printf("Here is the walk: ");
	i=0;
	while((*currentendhinge[secnum2]).endwalk[i] != 0){
		printf("%d ", (*currentendhinge[secnum2]).endwalk[i]);
		i++;
	}
	printf("\n\n");
*/




//	printf("finished recording an endhinge\n\n");


}


int dupcheck(int y, int z, int dup, int firstsection){

	//will recursively check if (*current_hinge_span[firstsection]).walk[i] will be a duplicate
	//check previous "starting" point prior to (x,y). (move there, then check if first step is a 1 or last step is a 2).
	if(dup==1){
		return dup;
	}

	//recursively check through all previous starting points before (y,z)   ordering of starting points is (0,0) (1,0) .... (L,0) (0,1) (1,1) ... (L,1) ... (0,2)... (L,2) ... (L,M)
	//move to new point, then check if first step is 1, or last step is 2. if it is, then we have a duplicate. (it must've already been generated)
	if(y>0){
		y--;
	}
	else if(z>0){
		y=L;
		z--;
	}
	else{//y==0 and z==0;
		return dup;
	}
	//now found proper y,z we need to test.


	int cury=(*current_hinge_span[firstsection]).startpoint[0];
	int curz=(*current_hinge_span[firstsection]).startpoint[1];
	int curx=0;

//	printf("checking for duplicate at point (%d, %d), startpoint is (%d,%d)\n", y, z, cury, curz);

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

	//move along dupcheckwalk[] until we get to (0,y,z)
//	printf("going to start moving through dupcheckwalk[]\n");
	i=0;
	while(((cury!=y || curz!=z) || curx!=0) || (dupcheckwalk[i-1]<0 && dupcheckwalk[i]<0)){	//while not at (0,y,z),
//		printf("curx=%d, cury=%d, curz=%d. Not at (0, %d, %d) yet!!!\n", curx, cury, curz, y, z);
		if(!((cury!=y || curz!=z) || curx!=0) && (dupcheckwalk[i-1]<0 && dupcheckwalk[i]<0)){
			curious++;
//			printf("At check point (0,%d,%d), but pretty sure it's a fake. curious=%d.\n", y, z, curious);
//			if(curious==4){
//				exit(1);
//			}
		}
		if(dupcheckwalk[i]==1){
			curx++;
		}
		else if(dupcheckwalk[i]==2){
			curx--;
		}
		else if(dupcheckwalk[i]==3 || dupcheckwalk[i]==-3){
			cury++;
		}
		else if(dupcheckwalk[i]==4 || dupcheckwalk[i]==-4){
			cury--;
		}
		else if(dupcheckwalk[i]==5 || dupcheckwalk[i]==-5){
			curz++;
		}
		else if(dupcheckwalk[i]==6 || dupcheckwalk[i]==-6){
			curz--;
		}
		else if(dupcheckwalk[i]==0){	//reached end of the walk[] (never gets to (0, y, z))
			i=i+1;	//just to be safe
			break;
		}
		i++;
//		printf("curx=%d, cury=%d, curz=%d\n", curx, cury, curz);
	}


	//now at (0,y,z) (before executing dupcheckwalk[i])
	if(dupcheckwalk[i-1]==2){	//just moved negx, so we have a duplicate
//		printf("duplicate found with last move was negx (y=%d, z=%d)\n", y, z);
		return 1;	//we have a duplicate
	}
	else if(dupcheckwalk[i]==1){	//next move was to move posx (add the +1 because of the neg.)
//		printf("duplicate found with next move of posx(y=%d, z=%d)\n", y, z);
		return 1;	//we have a duplicate
	}
	else{
		dup=dupcheck(y, z, dup, firstsection);
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









/***************************************************************************/







#include "../utils/unsgnlong_vecalloc.c" /*unsigned long int unsgnlong_vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with unsigned long integer enteries */

#include "../utils/vecalloc.c"		/*vec_ent * vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with enteries of type vec_ent */

#include "../utils/matalloc.c"		/*mat_ent **matalloc(int rowlow, int rowhigh, int collow, int colhigh)*/
/* used for dynamically allocating memory for matrices with enteries of type mat_ent */



#include "Improvedaddtowalk.c"
#include "ImprovedRecurGen_TwoString.c"
#include "bfacf_shrink_v2.c"

	









