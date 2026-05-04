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
#define	totalspan 150			/* span of generated polygon (z-direction)	*/
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
	int sixcone=0;
	int sixcone_s=0;
	int unreg=0;
	int unknown=0;
	char line[20];
	//	char dummy[20];

	int totalfiles=tot_num_files;
	int filenum=1;
	char filename[100];	//this will hold the filename of the current text file of knottypes that will be counted
	int counter=0;
	long unsigned int grand=0;


	for(filenum=1; filenum<=totalfiles; filenum++){	//for each knot file
	  counter=0;
	  sprintf(filename, "MCpolysHamL%dM%dspan%drun%dnum%d.txt_knot", L, M, totalspan, runnum, filenum); //name of knot file
		FILE* fp = fopen(filename, "r");	//open "filename" (a knot file)

		if(fp != NULL){
			printf("file %s opened successfully\n", filename);
		}
		else{
		  printf("file pointer is pointing to NULL. filename was %s\n", filename);
			break;
		}

		while( fgets(line, 20, fp)!=NULL ){	//while there's a line to read		fgets(line, 7, fp);
			if( strcmp(line, "0_1\n")==0 ){	//we have an unknot
				unknot++;
			}
			else if ( strcmp(line, "3_1\n")==0 ){	//we have a positive trefoil
				postref++;
			}
			else if ( strcmp(line, "3_1s\n")==0 ){	//we have a negative trefoil
				negtref++;
			}
			else if ( strcmp(line, "4_1\n")==0 ){	//we have a figure eight
				figeight++;
			}
			else if ( strcmp(line, "5_1\n")==0 ){	//we have a 5_1
				fiveone++;
			}
			else if ( strcmp(line, "5_1s\n")==0 ){	//we have a 5_1s
				fiveone_s++;
			}
			else if ( strcmp(line, "5_2\n")==0 ){	//we have a 5_2
				fivetwo++;
			}
			else if ( strcmp(line, "5_2s\n")==0 ){	//we have a 5_2s
				fivetwo_s++;
			}
			else if ( strcmp(line, "6_1\n")==0 ){	//we have a 6_1
				sixone++;
			}
			else if ( strcmp(line, "6_1s\n")==0 ){	//we have a 6_1s
				sixone_s++;
			}
			else if ( strcmp(line, "6_2\n")==0 ){
			  sixtwo++;
			}
			else if ( strcmp(line, "6_2s\n")==0 ){
			  sixtwo_s++;
			}
			else if ( strcmp(line, "6c_2\n")==0 ){	//we have a 6c_2
				sixctwo++;
			}
			else if ( strcmp(line, "6_3\n")==0 ){
			  sixthree++;
			}
			else if ( strcmp(line, "7_1\n")==0 ){
			  sevenone++;
			}
			else if ( strcmp(line, "7_1s\n")==0 ){
			  sevenone_s++;
			}
			else if ( strcmp(line, "7_2\n")==0 ){
			  seventwo++;
			}
			else if ( strcmp(line, "7_2s\n")==0 ){
			  seventwo_s++;
			}
			else if ( strcmp(line, "7_3\n")==0 ){
			  seventhree++;
			}
			else if ( strcmp(line, "7_3s\n")==0 ){
			  seventhree_s++;
			}
			else if ( strcmp(line, "7_4\n")==0 ){
			  sevenfour++;
			}
			else if ( strcmp(line, "7_4s\n")==0 ){
			  sevenfour_s++;
			}
			else if ( strcmp(line, "7_5\n")==0 ){
			  sevenfive++;
			}
			else if ( strcmp(line, "7_5s\n")==0 ){
			  sevenfive_s++;
			}
			else if ( strcmp(line, "7_6\n")==0 ){
			  sevensix++;
			}
			else if ( strcmp(line, "7_6s\n")==0 ){
			  sevensix_s++;
			}
			else if ( strcmp(line, "7_7\n")==0 ){
			  sevenseven++;
			}
			else if ( strcmp(line, "7_7s\n")==0 ){
			  sevenseven_s++;
			}
			else if ( strcmp(line, "7c_1\n")==0 ){
			  sevencone++;
			}
			else if ( strcmp(line, "7c_1s\n")==0 ){
			  sevencone_s++;
			}
			else if ( strcmp(line, "8_1\n")==0 ){
			  eightone++;
			}
			else if ( strcmp(line, "8_1s\n")==0 ){
			  eightone_s++;
			}
			else if ( strcmp(line, "8_2\n")==0 ){
			  eighttwo++;
			}
			else if ( strcmp(line, "8_2s\n")==0 ){
			  eighttwo_s++;
			}
			else if ( strcmp(line, "8_3\n")==0 ){
			  eightthree++;
			}
			else if ( strcmp(line, "8_3s\n")==0 ){
			  eightthree_s++;
			}
			else if ( strcmp(line, "8_4\n")==0 ){
			  eightfour++;
			}
			else if ( strcmp(line, "8_4s\n")==0 ){
			  eightfour_s++;
			}
			else if ( strcmp(line, "8_5\n")==0 ){
			  eightfive++;
			}
			else if ( strcmp(line, "8_5s\n")==0 ){
			  eightfive_s++;
			}
			else if ( strcmp(line, "8_6\n")==0 ){
			  eightsix++;
			}
			else if ( strcmp(line, "8_6s\n")==0 ){
			  eightsix_s++;
			}
			else if ( strcmp(line, "8_7\n")==0 ){
			  eightseven++;
			}
			else if ( strcmp(line, "8_7s\n")==0 ){
			  eightseven_s++;
			}
			else if ( strcmp(line, "8_8\n")==0 ){
			  eighteight++;
			}
			else if ( strcmp(line, "8_8s\n")==0 ){
			  eighteight_s++;
			}
			else if ( strcmp(line, "8_9\n")==0 ){
			  eightnine++;
			}
			else if ( strcmp(line, "8_9s\n")==0 ){
			  eightnine_s++;
			}
			else if ( strcmp(line, "8_10\n")==0 ){
			  eightten++;
			}
			else if ( strcmp(line, "8_10s\n")==0 ){
			  eightten_s++;
			}
			else if ( strcmp(line, "8_11\n")==0 ){
			  eighteleven++;
			}
			else if ( strcmp(line, "8_11s\n")==0 ){
			  eighteleven_s++;
			}
			else if ( strcmp(line, "8_12\n")==0 ){
			  eighttwelve++;
			}
			else if ( strcmp(line, "8_12s\n")==0 ){
			  eighttwelve_s++;
			}
			else if ( strcmp(line, "8_13\n")==0 ){
			  eightthirteen++;
			}
			else if ( strcmp(line, "8_13s\n")==0 ){
			  eightthirteen_s++;
			}
			else if ( strcmp(line, "8_14\n")==0 ){
			  eightfourteen++;
			}
			else if ( strcmp(line, "8_14s\n")==0 ){
			  eightfourteen_s++;
			}
			else if ( strcmp(line, "8_15\n")==0 ){
			  eightfifteen++;
			}
			else if ( strcmp(line, "8_15s\n")==0 ){
			  eightfifteen_s++;
			}
			else if ( strcmp(line, "8_16\n")==0 ){
			  eightsixteen++;
			}
			else if ( strcmp(line, "8_16s\n")==0 ){
			  eightsixteen_s++;
			}
			else if ( strcmp(line, "8_17\n")==0 ){
			  eightseventeen++;
			}
			else if ( strcmp(line, "8_17s\n")==0 ){
			  eightseventeen_s++;
			}
			else if ( strcmp(line, "8_18\n")==0 ){
			  eighteightteen++;
			}
			else if ( strcmp(line, "8_18s\n")==0 ){
			  eighteightteen_s++;
			}
			else if ( strcmp(line, "8_19\n")==0 ){
				eightnineteen++;
			}
			else if ( strcmp(line, "8_19s\n")==0 ){
				eightnineteen_s++;
			}
			else if ( strcmp(line, "8_20\n")==0 ){
			  eighttwenty++;
			}
			else if ( strcmp(line, "8_20s\n")==0 ){
			  eighttwenty_s++;
			}
			else if ( strcmp(line, "8_21\n")==0 ){
			  eighttwentyone++;
			}
			else if ( strcmp(line, "8_21s\n")==0 ){
			  eighttwentyone_s++;
			}
			else if ( strcmp(line, "8c_1\n")==0 ){
			  eightcone++;
			}
			else if ( strcmp(line, "8c_1s\n")==0 ){
			  eightcone_s++;
			}
			else if ( strcmp(line, "8c_2\n")==0 ){
			  eightctwo++;
			}
			else if ( strcmp(line, "8c_2s\n")==0 ){
			  eightctwo_s++;
			}
			else if ( strcmp(line, "8c_3\n")==0 ){
			  eightcthree++;
			}
			else if ( strcmp(line, "8c_3s\n")==0 ){
			  eightcthree_s++;
			}
			else if ( strcmp(line, "8c_4\n")==0 ){
			  eightcfour++;
			}
			else if ( strcmp(line, "8c_4s\n")==0 ){
			  eightcfour_s++;
			}
			else if ( strcmp(line, "8c_5\n")==0 ){
			  eightcfive++;
			}
			else if ( strcmp(line, "8c_5s\n")==0 ){
			  eightcfive_s++;
			}
			else if ( strcmp(line, "6c_1\n")==0 ){
			  sixcone++;
			}
			else if ( strcmp(line, "6c_1s\n")==0){
			  sixcone_s++;
			}
			else if ( strcmp(line, "unknown\n")==0){
			  unknown++;
			}
			else{	//unregistered knot type
				printf("unregistered found, it is: %s\n", line);
//				scanf("%s", dummy);
				unreg++;
			}
			counter++;
			grand++;
		}
		if(counter<10000){
		  printf("file %d contains <1000000 polys, it only contains %d\n", filenum, counter);
		}



	}

	printf("0_1: %d\n", unknot);
	printf("3_1: %d\n", postref);
	printf("3_1s: %d\n", negtref);
	printf("4_1: %d\n", figeight);
	printf("5_1: %d\n", fiveone);
	printf("5_1s: %d\n", fiveone_s);
	printf("5_2: %d\n", fivetwo);
	printf("5_2s: %d\n", fivetwo_s);
	printf("6_1: %d\n", sixone);
	printf("6_1s: %d\n", sixone_s);
	printf("6_2: %d\n", sixtwo);
	printf("6_2s: %d\n", sixtwo_s);
	printf("6_3: %d\n", sixthree);
	printf("6c_1: %d\n", sixcone);
	printf("6c_1s: %d\n", sixcone_s);
	printf("6c_2: %d\n", sixctwo);
	printf("7_1: %d\n", sevenone);
	printf("7_1s: %d\n", sevenone_s);
	printf("7_2: %d\n", seventwo);
	printf("7_2s: %d\n", seventwo_s);
	printf("7_3: %d\n", seventhree);
	printf("7_3s: %d\n", seventhree_s);
	printf("7_4: %d\n", sevenfour);
	printf("7_4s: %d\n", sevenfour_s);
	printf("7_5: %d\n", sevenfive);
	printf("7_5s: %d\n", sevenfive_s);
	printf("7_6: %d\n", sevensix);
	printf("7_6s: %d\n", sevensix_s);
	printf("7_7: %d\n", sevenseven);
	printf("7_7s: %d\n", sevenseven_s);
	printf("7c_1: %d\n", sevencone);
	printf("7c_1s: %d\n", sevencone_s);
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
	printf("unknown: %d\n", unknown);
	printf("unregistered: %d\n", unreg);
	printf("total polys read: %lu\n", grand);


	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	system("date");  /* prints the date and time */

	return 0;
}				/* end of main */
/***************************************************************************/


