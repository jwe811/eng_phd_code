/* This program will count the total knottypes of a series of knot files	*/


#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>


/** START OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDYIED **/

#define	L 4				/* number of horizontal edges  (x-direction)	*/
#define	M 1				/* number of vertical edges (y-direction)	*/
#define	totalspan 50			/* span of generated polygon (z-direction)	*/
#define	tot_num_files 100
#define runnum 1

#define	vM (M+1)				/* number of vertices in vertical direction   = M + 1      */
#define	vL (L+1)				/* number of vertices in horizontal direction = L + 1      */


/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/
int main(void)
{
	system("date");  /* prints the date and time */
	clock();	/* Start clock to time program */

	//	int i, j, k;

	int unknot=0;

	int postref=0;
	int negtref=0;

	int threeoneG_threeoneG=0;
	int threeoneG_threeoneL=0;
	int threeoneL_threeoneL=0;

	int threeone_sG_threeone_sG=0;
	int threeone_sG_threeone_sL=0;
	int threeone_sL_threeone_sL=0;

	int threeoneG_threeone_sG=0;
	int threeoneG_threeone_sL=0;
	int threeoneL_threeone_sG=0;
	int threeoneL_threeone_sL=0;

	int figeight=0;
	int fiveone=0;
	int fiveone_s=0;
	int fivetwo=0;
	int fivetwo_s=0;
	int sixone=0;
	int sixone_s=0;
	int sixtwo=0;
	int sixtwo_s=0;
	int sixthree=0;

	int sixcone=0;
	int sixcone_s=0;
	int sixctwo=0;

	int sevenone=0;
	int sevenone_s=0;
	int seventwo=0;
	int seventwo_s=0;
	int seventhree=0;
	int seventhree_s=0;
	int sevenfour=0;
	int sevenfour_s=0;
	int sevenfive=0;
	int sevenfive_s=0;
	int sevensix=0;
	int sevensix_s=0;
	int sevenseven=0;
	int sevenseven_s=0;
	int sevencone=0;
	int sevencone_s=0;
	int eightone=0;
	int eightone_s=0;
	int eighttwo=0;
	int eighttwo_s=0;
	int eightthree=0;
	int eightthree_s=0;
	int eightfour=0;
	int eightfour_s=0;
	int eightfive=0;
	int eightfive_s=0;
	int eightsix=0;
	int eightsix_s=0;
	int eightseven=0;
	int eightseven_s=0;
	int eighteight=0;
	int eighteight_s=0;
	int eightnine=0;
	int eightnine_s=0;
	int eightten=0;
	int eightten_s=0;
	int eighteleven=0;
	int eighteleven_s=0;
	int eighttwelve=0;
	int eighttwelve_s=0;
	int eightthirteen=0;
	int eightthirteen_s=0;
	int eightfourteen=0;
	int eightfourteen_s=0;
	int eightfifteen=0;
	int eightfifteen_s=0;
	int eightsixteen=0;
	int eightsixteen_s=0;
	int eightseventeen=0;
	int eightseventeen_s=0;
	int eighteightteen=0;
	int eighteightteen_s=0;
	int eightnineteen=0;
	int eightnineteen_s=0;
	int eighttwenty=0;
	int eighttwenty_s=0;
	int eighttwentyone=0;
	int eighttwentyone_s=0;
	int eightcone=0;
	int eightcone_s=0;
	int eightctwo=0;
	int eightctwo_s=0;
	int eightcthree=0;
	int eightcthree_s=0;
	int eightcfour=0;
	int eightcfour_s=0;
	int eightcfive=0;
	int eightcfive_s=0;

	int postrefL=0;
	int negtrefL=0;
	int figeightL=0;
	int fiveoneL=0;
	int fiveone_sL=0;
	int fivetwoL=0;
	int fivetwo_sL=0;
	int sixoneL=0;
	int sixone_sL=0;
	int sixtwoL=0;
	int sixtwo_sL=0;
	int sixthreeL=0;

	int sixconeL=0;
	int sixcone_sL=0;
	int sixctwoL=0;

	int sevenoneL=0;
	int sevenone_sL=0;
	int seventwoL=0;
	int seventwo_sL=0;
	int seventhreeL=0;
	int seventhree_sL=0;
	int sevenfourL=0;
	int sevenfour_sL=0;
	int sevenfiveL=0;
	int sevenfive_sL=0;
	int sevensixL=0;
	int sevensix_sL=0;
	int sevensevenL=0;
	int sevenseven_sL=0;
	int sevenconeL=0;
	int sevencone_sL=0;
	int eightoneL=0;
	int eightone_sL=0;
	int eighttwoL=0;
	int eighttwo_sL=0;
	int eightthreeL=0;
	int eightthree_sL=0;
	int eightfourL=0;
	int eightfour_sL=0;
	int eightfiveL=0;
	int eightfive_sL=0;
	int eightsixL=0;
	int eightsix_sL=0;
	int eightsevenL=0;
	int eightseven_sL=0;
	int eighteightL=0;
	int eighteight_sL=0;
	int eightnineL=0;
	int eightnine_sL=0;
	int eighttenL=0;
	int eightten_sL=0;
	int eightelevenL=0;
	int eighteleven_sL=0;
	int eighttwelveL=0;
	int eighttwelve_sL=0;
	int eightthirteenL=0;
	int eightthirteen_sL=0;
	int eightfourteenL=0;
	int eightfourteen_sL=0;
	int eightfifteenL=0;
	int eightfifteen_sL=0;
	int eightsixteenL=0;
	int eightsixteen_sL=0;
	int eightseventeenL=0;
	int eightseventeen_sL=0;
	int eighteightteenL=0;
	int eighteightteen_sL=0;
	int eightnineteenL=0;
	int eightnineteen_sL=0;
	int eighttwentyL=0;
	int eighttwenty_sL=0;
	int eighttwentyoneL=0;
	int eighttwentyone_sL=0;
	int eightconeL=0;
	int eightcone_sL=0;
	int eightctwoL=0;
	int eightctwo_sL=0;
	int eightcthreeL=0;
	int eightcthree_sL=0;
	int eightcfourL=0;
	int eightcfour_sL=0;
	int eightcfiveL=0;
	int eightcfive_sL=0;


	int unreg=0;
	int unknown=0;
	int unknownL=0;
	char line[20];
	//	char dummy[20];

	int totalfiles=tot_num_files;
	int filenum=1;
	char filename[100];	//this will hold the filename of the current text file of knottypes that will be counted
	int counter=0;
	long unsigned int grand=0;


	for(filenum=1; filenum<=totalfiles; filenum++){	//for each knot file
	  counter=0;
	  sprintf(filename, "MCpolysHamL%dM%dspan%drun%dnum%d_LP.txt-id", L, M, totalspan, runnum, filenum); //name of knot file
		FILE* fp = fopen(filename, "r");	//open "filename" (a knot file)

		if(fp != NULL){
			printf("file %s opened successfully\n", filename);
		}
		else{
		  printf("file pointer is pointing to NULL. filename was %s\n", filename);
			break;
		}

		while( fgets(line, 20, fp)!=NULL ){	//while there's a line to read		fgets(line, 7, fp);
			if( strcmp(line, "\n")==0 ){	//we have an unknot
				unknot++;
			}
////////////////////////////
			else if ( strcmp(line, "3_1G\n")==0 ){
				postref++;
			}
			else if ( strcmp(line, "3_1sG\n")==0 ){
				negtref++;
			}
////////////
			else if ( strcmp(line, "3_1G#3_1G\n")==0 ){
				threeoneG_threeoneG++;
			}
			else if ( strcmp(line, "3_1G#3_1L\n")==0 ){
				threeoneG_threeoneL++;
			}
			else if ( strcmp(line, "3_1L#3_1G\n")==0 ){
				threeoneG_threeoneL++;
			}
			else if ( strcmp(line, "3_1L#3_1L\n")==0 ){
				threeoneL_threeoneL++;
			}

			else if ( strcmp(line, "3_1sG#3_1sG\n")==0 ){
				threeone_sG_threeone_sG++;
			}
			else if ( strcmp(line, "3_1sG#3_1sL\n")==0 ){
				threeone_sG_threeone_sL++;
			}
			else if ( strcmp(line, "3_1sL#3_1sG\n")==0 ){
				threeone_sG_threeone_sL++;
			}
			else if ( strcmp(line, "3_1sL#3_1sL\n")==0 ){
				threeone_sL_threeone_sL++;
			}

			else if ( strcmp(line, "3_1G#3_1sG\n")==0 ){
				threeoneG_threeone_sG++;
			}
			else if ( strcmp(line, "3_1G#3_1sL\n")==0 ){
				threeoneG_threeone_sL++;
			}
			else if ( strcmp(line, "3_1L#3_1sG\n")==0 ){
				threeoneL_threeone_sG++;
			}
			else if ( strcmp(line, "3_1L#3_1sL\n")==0 ){
				threeoneL_threeone_sL++;
			}

			else if ( strcmp(line, "3_1sG#3_1G\n")==0 ){
				threeoneG_threeone_sG++;
			}
			else if ( strcmp(line, "3_1sL#3_1G\n")==0 ){
				threeoneG_threeone_sL++;
			}
			else if ( strcmp(line, "3_1sG#3_1L\n")==0 ){
				threeoneL_threeone_sG++;
			}
			else if ( strcmp(line, "3_1sL#3_1G\n")==0 ){
				threeoneL_threeone_sL++;
			}
/////////////

			else if ( strcmp(line, "4_1G\n")==0 ){	//we have a figure eight
				figeight++;
			}
			else if ( strcmp(line, "5_1G\n")==0 ){	//we have a 5_1
				fiveone++;
			}
			else if ( strcmp(line, "5_1sG\n")==0 ){	//we have a 5_1s
				fiveone_s++;
			}
			else if ( strcmp(line, "5_2G\n")==0 ){	//we have a 5_2
				fivetwo++;
			}
			else if ( strcmp(line, "5_2sG\n")==0 ){	//we have a 5_2s
				fivetwo_s++;
			}
			else if ( strcmp(line, "6_1G\n")==0 ){	//we have a 6_1
				sixone++;
			}
			else if ( strcmp(line, "6_1sG\n")==0 ){	//we have a 6_1s
				sixone_s++;
			}
			else if ( strcmp(line, "6_2G\n")==0 ){
			  sixtwo++;
			}
			else if ( strcmp(line, "6_2sG\n")==0 ){
			  sixtwo_s++;
			}
			else if ( strcmp(line, "6_3G\n")==0 ){
			  sixthree++;
			}

			else if ( strcmp(line, "6c_1G\n")==0 ){
			  sixcone++;
			}
			else if ( strcmp(line, "6c_1sG\n")==0 ){
			  sixcone_s++;
			}
			else if ( strcmp(line, "6c_2G\n")==0 ){
			  sixctwo++;
			}

			else if ( strcmp(line, "7_1G\n")==0 ){
			  sevenone++;
			}
			else if ( strcmp(line, "7_1sG\n")==0 ){
			  sevenone_s++;
			}
			else if ( strcmp(line, "7_2G\n")==0 ){
			  seventwo++;
			}
			else if ( strcmp(line, "7_2sG\n")==0 ){
			  seventwo_s++;
			}
			else if ( strcmp(line, "7_3G\n")==0 ){
			  seventhree++;
			}
			else if ( strcmp(line, "7_3sG\n")==0 ){
			  seventhree_s++;
			}
			else if ( strcmp(line, "7_4G\n")==0 ){
			  sevenfour++;
			}
			else if ( strcmp(line, "7_4sG\n")==0 ){
			  sevenfour_s++;
			}
			else if ( strcmp(line, "7_5G\n")==0 ){
			  sevenfive++;
			}
			else if ( strcmp(line, "7_5sG\n")==0 ){
			  sevenfive_s++;
			}
			else if ( strcmp(line, "7_6G\n")==0 ){
			  sevensix++;
			}
			else if ( strcmp(line, "7_6sG\n")==0 ){
			  sevensix_s++;
			}
			else if ( strcmp(line, "7_7G\n")==0 ){
			  sevenseven++;
			}
			else if ( strcmp(line, "7_7sG\n")==0 ){
			  sevenseven_s++;
			}
			else if ( strcmp(line, "7c_1G\n")==0 ){
			  sevencone++;
			}
			else if ( strcmp(line, "7c_1sG\n")==0 ){
			  sevencone_s++;
			}
			else if ( strcmp(line, "8_1G\n")==0 ){
			  eightone++;
			}
			else if ( strcmp(line, "8_1sG\n")==0 ){
			  eightone_s++;
			}
			else if ( strcmp(line, "8_2G\n")==0 ){
			  eighttwo++;
			}
			else if ( strcmp(line, "8_2sG\n")==0 ){
			  eighttwo_s++;
			}
			else if ( strcmp(line, "8_3G\n")==0 ){
			  eightthree++;
			}
			else if ( strcmp(line, "8_3sG\n")==0 ){
			  eightthree_s++;
			}
			else if ( strcmp(line, "8_4G\n")==0 ){
			  eightfour++;
			}
			else if ( strcmp(line, "8_4sG\n")==0 ){
			  eightfour_s++;
			}
			else if ( strcmp(line, "8_5G\n")==0 ){
			  eightfive++;
			}
			else if ( strcmp(line, "8_5sG\n")==0 ){
			  eightfive_s++;
			}
			else if ( strcmp(line, "8_6G\n")==0 ){
			  eightsix++;
			}
			else if ( strcmp(line, "8_6sG\n")==0 ){
			  eightsix_s++;
			}
			else if ( strcmp(line, "8_7G\n")==0 ){
			  eightseven++;
			}
			else if ( strcmp(line, "8_7sG\n")==0 ){
			  eightseven_s++;
			}
			else if ( strcmp(line, "8_8G\n")==0 ){
			  eighteight++;
			}
			else if ( strcmp(line, "8_8sG\n")==0 ){
			  eighteight_s++;
			}
			else if ( strcmp(line, "8_9G\n")==0 ){
			  eightnine++;
			}
			else if ( strcmp(line, "8_9sG\n")==0 ){
			  eightnine_s++;
			}
			else if ( strcmp(line, "8_10G\n")==0 ){
			  eightten++;
			}
			else if ( strcmp(line, "8_10sG\n")==0 ){
			  eightten_s++;
			}
			else if ( strcmp(line, "8_11G\n")==0 ){
			  eighteleven++;
			}
			else if ( strcmp(line, "8_11sG\n")==0 ){
			  eighteleven_s++;
			}
			else if ( strcmp(line, "8_12G\n")==0 ){
			  eighttwelve++;
			}
			else if ( strcmp(line, "8_12sG\n")==0 ){
			  eighttwelve_s++;
			}
			else if ( strcmp(line, "8_13G\n")==0 ){
			  eightthirteen++;
			}
			else if ( strcmp(line, "8_13sG\n")==0 ){
			  eightthirteen_s++;
			}
			else if ( strcmp(line, "8_14G\n")==0 ){
			  eightfourteen++;
			}
			else if ( strcmp(line, "8_14sG\n")==0 ){
			  eightfourteen_s++;
			}
			else if ( strcmp(line, "8_15G\n")==0 ){
			  eightfifteen++;
			}
			else if ( strcmp(line, "8_15sG\n")==0 ){
			  eightfifteen_s++;
			}
			else if ( strcmp(line, "8_16G\n")==0 ){
			  eightsixteen++;
			}
			else if ( strcmp(line, "8_16sG\n")==0 ){
			  eightsixteen_s++;
			}
			else if ( strcmp(line, "8_17G\n")==0 ){
			  eightseventeen++;
			}
			else if ( strcmp(line, "8_17sG\n")==0 ){
			  eightseventeen_s++;
			}
			else if ( strcmp(line, "8_18G\n")==0 ){
			  eighteightteen++;
			}
			else if ( strcmp(line, "8_18sG\n")==0 ){
			  eighteightteen_s++;
			}
			else if ( strcmp(line, "8_19G\n")==0 ){
				eightnineteen++;
			}
			else if ( strcmp(line, "8_19sG\n")==0 ){
				eightnineteen_s++;
			}
			else if ( strcmp(line, "8_20G\n")==0 ){
			  eighttwenty++;
			}
			else if ( strcmp(line, "8_20sG\n")==0 ){
			  eighttwenty_s++;
			}
			else if ( strcmp(line, "8_21G\n")==0 ){
			  eighttwentyone++;
			}
			else if ( strcmp(line, "8_21sG\n")==0 ){
			  eighttwentyone_s++;
			}
			else if ( strcmp(line, "8c_1G\n")==0 ){
			  eightcone++;
			}
			else if ( strcmp(line, "8c_1sG\n")==0 ){
			  eightcone_s++;
			}
			else if ( strcmp(line, "8c_2G\n")==0 ){
			  eightctwo++;
			}
			else if ( strcmp(line, "8c_2sG\n")==0 ){
			  eightctwo_s++;
			}
			else if ( strcmp(line, "8c_3G\n")==0 ){
			  eightcthree++;
			}
			else if ( strcmp(line, "8c_3sG\n")==0 ){
			  eightcthree_s++;
			}
			else if ( strcmp(line, "8c_4G\n")==0 ){
			  eightcfour++;
			}
			else if ( strcmp(line, "8c_4sG\n")==0 ){
			  eightcfour_s++;
			}
			else if ( strcmp(line, "8c_5G\n")==0 ){
			  eightcfive++;
			}
			else if ( strcmp(line, "8c_5sG\n")==0 ){
			  eightcfive_s++;
			}
///////////////////////////////////////
			else if ( strcmp(line, "3_1L\n")==0 ){	//we have a positive trefoil
				postrefL++;
			}
			else if ( strcmp(line, "3_1sL\n")==0 ){	//we have a negative trefoil
				negtrefL++;
			}
			else if ( strcmp(line, "4_1L\n")==0 ){	//we have a figure eight
				figeightL++;
			}
			else if ( strcmp(line, "5_1L\n")==0 ){	//we have a 5_1
				fiveoneL++;
			}
			else if ( strcmp(line, "5_1sL\n")==0 ){	//we have a 5_1s
				fiveone_sL++;
			}
			else if ( strcmp(line, "5_2L\n")==0 ){	//we have a 5_2
				fivetwoL++;
			}
			else if ( strcmp(line, "5_2sL\n")==0 ){	//we have a 5_2s
				fivetwo_sL++;
			}
			else if ( strcmp(line, "6_1L\n")==0 ){	//we have a 6_1
				sixoneL++;
			}
			else if ( strcmp(line, "6_1sL\n")==0 ){	//we have a 6_1s
				sixone_sL++;
			}
			else if ( strcmp(line, "6_2L\n")==0 ){
			  sixtwoL++;
			}
			else if ( strcmp(line, "6_2sL\n")==0 ){
			  sixtwo_sL++;
			}
			else if ( strcmp(line, "6_3L\n")==0 ){
			  sixthreeL++;
			}
			else if ( strcmp(line, "6c_1L\n")==0 ){
			  sixconeL++;
			}
			else if ( strcmp(line, "6c_1sL\n")==0 ){
			  sixcone_sL++;
			}
			else if ( strcmp(line, "6c_2L\n")==0 ){
			  sixctwoL++;
			}
			else if ( strcmp(line, "7_1L\n")==0 ){
			  sevenoneL++;
			}
			else if ( strcmp(line, "7_1sL\n")==0 ){
			  sevenone_sL++;
			}
			else if ( strcmp(line, "7_2L\n")==0 ){
			  seventwoL++;
			}
			else if ( strcmp(line, "7_2sL\n")==0 ){
			  seventwo_sL++;
			}
			else if ( strcmp(line, "7_3L\n")==0 ){
			  seventhreeL++;
			}
			else if ( strcmp(line, "7_3sL\n")==0 ){
			  seventhree_sL++;
			}
			else if ( strcmp(line, "7_4L\n")==0 ){
			  sevenfourL++;
			}
			else if ( strcmp(line, "7_4sL\n")==0 ){
			  sevenfour_sL++;
			}
			else if ( strcmp(line, "7_5L\n")==0 ){
			  sevenfiveL++;
			}
			else if ( strcmp(line, "7_5sL\n")==0 ){
			  sevenfive_sL++;
			}
			else if ( strcmp(line, "7_6L\n")==0 ){
			  sevensixL++;
			}
			else if ( strcmp(line, "7_6sL\n")==0 ){
			  sevensix_sL++;
			}
			else if ( strcmp(line, "7_7L\n")==0 ){
			  sevensevenL++;
			}
			else if ( strcmp(line, "7_7sL\n")==0 ){
			  sevenseven_sL++;
			}
			else if ( strcmp(line, "7c_1L\n")==0 ){
			  sevenconeL++;
			}
			else if ( strcmp(line, "7c_1sL\n")==0 ){
			  sevencone_sL++;
			}
			else if ( strcmp(line, "8_1L\n")==0 ){
			  eightoneL++;
			}
			else if ( strcmp(line, "8_1sL\n")==0 ){
			  eightone_sL++;
			}
			else if ( strcmp(line, "8_2L\n")==0 ){
			  eighttwoL++;
			}
			else if ( strcmp(line, "8_2sL\n")==0 ){
			  eighttwo_sL++;
			}
			else if ( strcmp(line, "8_3L\n")==0 ){
			  eightthreeL++;
			}
			else if ( strcmp(line, "8_3sL\n")==0 ){
			  eightthree_sL++;
			}
			else if ( strcmp(line, "8_4L\n")==0 ){
			  eightfourL++;
			}
			else if ( strcmp(line, "8_4sL\n")==0 ){
			  eightfour_sL++;
			}
			else if ( strcmp(line, "8_5L\n")==0 ){
			  eightfiveL++;
			}
			else if ( strcmp(line, "8_5sL\n")==0 ){
			  eightfive_sL++;
			}
			else if ( strcmp(line, "8_6L\n")==0 ){
			  eightsixL++;
			}
			else if ( strcmp(line, "8_6sL\n")==0 ){
			  eightsix_sL++;
			}
			else if ( strcmp(line, "8_7L\n")==0 ){
			  eightsevenL++;
			}
			else if ( strcmp(line, "8_7sL\n")==0 ){
			  eightseven_sL++;
			}
			else if ( strcmp(line, "8_8L\n")==0 ){
			  eighteightL++;
			}
			else if ( strcmp(line, "8_8sL\n")==0 ){
			  eighteight_sL++;
			}
			else if ( strcmp(line, "8_9L\n")==0 ){
			  eightnineL++;
			}
			else if ( strcmp(line, "8_9sL\n")==0 ){
			  eightnine_sL++;
			}
			else if ( strcmp(line, "8_10L\n")==0 ){
			  eighttenL++;
			}
			else if ( strcmp(line, "8_10sL\n")==0 ){
			  eightten_sL++;
			}
			else if ( strcmp(line, "8_11L\n")==0 ){
			  eightelevenL++;
			}
			else if ( strcmp(line, "8_11sL\n")==0 ){
			  eighteleven_sL++;
			}
			else if ( strcmp(line, "8_12L\n")==0 ){
			  eighttwelveL++;
			}
			else if ( strcmp(line, "8_12sL\n")==0 ){
			  eighttwelve_sL++;
			}
			else if ( strcmp(line, "8_13L\n")==0 ){
			  eightthirteenL++;
			}
			else if ( strcmp(line, "8_13sL\n")==0 ){
			  eightthirteen_sL++;
			}
			else if ( strcmp(line, "8_14L\n")==0 ){
			  eightfourteenL++;
			}
			else if ( strcmp(line, "8_14sL\n")==0 ){
			  eightfourteen_sL++;
			}
			else if ( strcmp(line, "8_15L\n")==0 ){
			  eightfifteenL++;
			}
			else if ( strcmp(line, "8_15sL\n")==0 ){
			  eightfifteen_sL++;
			}
			else if ( strcmp(line, "8_16L\n")==0 ){
			  eightsixteenL++;
			}
			else if ( strcmp(line, "8_16sL\n")==0 ){
			  eightsixteen_sL++;
			}
			else if ( strcmp(line, "8_17L\n")==0 ){
			  eightseventeenL++;
			}
			else if ( strcmp(line, "8_17sL\n")==0 ){
			  eightseventeen_sL++;
			}
			else if ( strcmp(line, "8_18L\n")==0 ){
			  eighteightteenL++;
			}
			else if ( strcmp(line, "8_18sL\n")==0 ){
			  eighteightteen_sL++;
			}
			else if ( strcmp(line, "8_19L\n")==0 ){
				eightnineteenL++;
			}
			else if ( strcmp(line, "8_19sL\n")==0 ){
				eightnineteen_sL++;
			}
			else if ( strcmp(line, "8_20L\n")==0 ){
			  eighttwentyL++;
			}
			else if ( strcmp(line, "8_20sL\n")==0 ){
			  eighttwenty_sL++;
			}
			else if ( strcmp(line, "8_21L\n")==0 ){
			  eighttwentyoneL++;
			}
			else if ( strcmp(line, "8_21sL\n")==0 ){
			  eighttwentyone_sL++;
			}
			else if ( strcmp(line, "8c_1L\n")==0 ){
			  eightconeL++;
			}
			else if ( strcmp(line, "8c_1sL\n")==0 ){
			  eightcone_sL++;
			}
			else if ( strcmp(line, "8c_2L\n")==0 ){
			  eightctwoL++;
			}
			else if ( strcmp(line, "8c_2sL\n")==0 ){
			  eightctwo_sL++;
			}
			else if ( strcmp(line, "8c_3L\n")==0 ){
			  eightcthreeL++;
			}
			else if ( strcmp(line, "8c_3sL\n")==0 ){
			  eightcthree_sL++;
			}
			else if ( strcmp(line, "8c_4L\n")==0 ){
			  eightcfourL++;
			}
			else if ( strcmp(line, "8c_4sL\n")==0 ){
			  eightcfour_sL++;
			}
			else if ( strcmp(line, "8c_5L\n")==0 ){
			  eightcfiveL++;
			}
			else if ( strcmp(line, "8c_5sL\n")==0 ){
			  eightcfive_sL++;
			}
////////////////////////
			else if ( strcmp(line, "unknownG\n")==0){
			  unknown++;
			}
			else if ( strcmp(line, "unknownL\n")==0){
			  unknownL++;
			}
			else{	//unregistered knot type
				printf("unregistered found, it is: %s", line);
//				scanf("%s", dummy);
				unreg++;
			}
			counter++;
			grand++;
		}
		printf("file %d contains %d polys\n", filenum, counter);



	}

	printf("0_1: %d\n", unknot);
	printf("3_1G: %d\n", postref);
	printf("3_1L: %d\n", postrefL);
	printf("3_1sG: %d\n", negtref);
	printf("3_1sL: %d\n", negtrefL);

	printf("3_1G#3_1G: %d\n", threeoneG_threeoneG);
	printf("3_1G#3_1L: %d\n", threeoneG_threeoneL);
	printf("3_1L#3_1L: %d\n", threeoneL_threeoneL);

	printf("3_1sG#3_1sG: %d\n", threeone_sG_threeone_sG);
	printf("3_1sG#3_1sL: %d\n", threeone_sG_threeone_sL);
	printf("3_1sL#3_1sL: %d\n", threeone_sL_threeone_sL);

	printf("3_1G#3_1sG: %d\n", threeoneG_threeone_sG);
	printf("3_1G#3_1sL: %d\n", threeoneG_threeone_sL);
	printf("3_1L#3_1sG: %d\n", threeoneL_threeone_sG);
	printf("3_1L#3_1sL: %d\n", threeoneL_threeone_sL);

	printf("4_1G: %d\n", figeight);
	printf("4_1L: %d\n", figeightL);

	printf("5_1G: %d\n", fiveone);
	printf("5_1L: %d\n", fiveoneL);
	printf("5_1sG: %d\n", fiveone_s);
	printf("5_1sL: %d\n", fiveone_sL);
	printf("5_2G: %d\n", fivetwo);
	printf("5_2L: %d\n", fivetwoL);
	printf("5_2sG: %d\n", fivetwo_s);
	printf("5_2sL: %d\n", fivetwo_sL);

	printf("6_1G: %d\n", sixone);
	printf("6_1L: %d\n", sixoneL);
	printf("6_1sG: %d\n", sixone_s);
	printf("6_1sL: %d\n", sixone_sL);
	printf("6_2G: %d\n", sixtwo);
	printf("6_2L: %d\n", sixtwoL);
	printf("6_2sG: %d\n", sixtwo_s);
	printf("6_2sL: %d\n", sixtwo_sL);
	printf("6_3G: %d\n", sixthree);
	printf("6_3L: %d\n", sixthreeL);

	printf("6c_1G: %d\n", sixcone);
	printf("6c_1L: %d\n", sixconeL);
	printf("6c_1sG: %d\n", sixcone_s);
	printf("6c_1sL: %d\n", sixcone_sL);
	printf("6c_2G: %d\n", sixctwo);
	printf("6c_2L: %d\n", sixctwoL);

	printf("7_1G: %d\n", sevenone);
	printf("7_1L: %d\n", sevenoneL);
	printf("7_1sG: %d\n", sevenone_s);
	printf("7_1sL: %d\n", sevenone_sL);
	printf("7_2G: %d\n", seventwo);
	printf("7_2L: %d\n", seventwoL);
	printf("7_2sG: %d\n", seventwo_s);
	printf("7_2sL: %d\n", seventwo_sL);
	printf("7_3G: %d\n", seventhree);
	printf("7_3L: %d\n", seventhreeL);
	printf("7_3sG: %d\n", seventhree_s);
	printf("7_3sL: %d\n", seventhree_sL);
	printf("7_4G: %d\n", sevenfour);
	printf("7_4L: %d\n", sevenfourL);
	printf("7_4sG: %d\n", sevenfour_s);
	printf("7_4sL: %d\n", sevenfour_sL);
	printf("7_5G: %d\n", sevenfive);
	printf("7_5L: %d\n", sevenfiveL);
	printf("7_5sG: %d\n", sevenfive_s);
	printf("7_5sL: %d\n", sevenfive_sL);
	printf("7_6G: %d\n", sevensix);
	printf("7_6L: %d\n", sevensixL);
	printf("7_6sG: %d\n", sevensix_s);
	printf("7_6sL: %d\n", sevensix_sL);
	printf("7_7G: %d\n", sevenseven);
	printf("7_7L: %d\n", sevensevenL);
	printf("7_7sG: %d\n", sevenseven_s);
	printf("7_7sL: %d\n", sevenseven_sL);
	printf("7c_1G: %d\n", sevencone);
	printf("7c_1L: %d\n", sevenconeL);
	printf("7c_1sG: %d\n", sevencone_s);
	printf("7c_1sL: %d\n", sevencone_sL);

	printf("8_1: %d\n", eightone);
	printf("8_1s: %d\n", eightone_s);
	printf("8_2: %d\n", eighttwo);
	printf("8_2s: %d\n", eighttwo_s);
	printf("8_3: %d\n", eightthree);
	printf("8_3s: %d\n", eightthree_s);
	printf("8_4: %d\n", eightfour);
	printf("8_4s: %d\n", eightfour_s);
	printf("8_5: %d\n", eightfive);
	printf("8_5s: %d\n", eightfive_s);
	printf("8_6: %d\n", eightsix);
	printf("8_6s: %d\n", eightsix_s);
	printf("8_7: %d\n", eightseven);
	printf("8_7s: %d\n", eightseven_s);
	printf("8_8: %d\n", eighteight);
	printf("8_8s: %d\n", eighteight_s);
	printf("8_9: %d\n", eightnine);
	printf("8_9s: %d\n", eightnine_s);
	printf("8_10: %d\n", eightten);
	printf("8_10s: %d\n", eightten_s);
	printf("8_11: %d\n", eighteleven);
	printf("8_11s: %d\n", eighteleven_s);
	printf("8_12: %d\n", eighttwelve);
	printf("8_12s: %d\n", eighttwelve_s);
	printf("8_13: %d\n", eightthirteen);
	printf("8_13s: %d\n", eightthirteen_s);
	printf("8_14: %d\n", eightfourteen);
	printf("8_14s: %d\n", eightfourteen_s);
	printf("8_15: %d\n", eightfifteen);
	printf("8_15s: %d\n", eightfifteen_s);
	printf("8_16: %d\n", eightsixteen);
	printf("8_16s: %d\n", eightsixteen_s);
	printf("8_17: %d\n", eightseventeen);
	printf("8_17s: %d\n", eightseventeen_s);
	printf("8_18: %d\n", eighteightteen);
	printf("8_18s: %d\n", eighteightteen_s);
	printf("8_19: %d\n", eightnineteen);
	printf("8_19s: %d\n", eightnineteen_s);
	printf("8_20: %d\n", eighttwenty);
	printf("8_20s: %d\n", eighttwenty_s);
	printf("8_21: %d\n", eighttwentyone);
	printf("8_21s: %d\n", eighttwentyone_s);
	printf("8c_1: %d\n", eightcone);
	printf("8c_1s: %d\n", eightcone_s);
	printf("8c_2: %d\n", eightctwo);
	printf("8c_2s: %d\n", eightctwo_s);
	printf("8c_3: %d\n", eightcthree);
	printf("8c_3s: %d\n", eightcthree_s);
	printf("8c_4: %d\n", eightcfour);
	printf("8c_4s: %d\n", eightcfour_s);
	printf("8c_5: %d\n", eightcfive);
	printf("8c_5s: %d\n", eightcfive_s);
	printf("unknownG: %d\n", unknown);
	printf("unknownL: %d\n", unknownL);
	printf("unregistered: %d\n", unreg);
	printf("total polys read: %lu\n", grand);


	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	system("date");  /* prints the date and time */

	return 0;
}				/* end of main */
/***************************************************************************/


