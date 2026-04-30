//Takes goes through folders "span10 -> span20 -> ... -> span100", goes into their knot_list files, and tallies up the 9 crossing link-types


#include <stdio.h>		/* Standard input/output		*/
#include <stdlib.h>		/* Need Ascii to integer function	*/
#include <time.h>		/* Used to time the program		*/
#include <math.h>		/* Standard math functions (i.e. pow)	*/
#include <sys/time.h>
#include <string.h>



#define L 4
#define M 1
#define startingspan 125	//only span
#define runnum 1
#define totalfiles 100



int main(void) {
	printf("\nmain starting\n\n");

	char filename[200];		//files with the links

	char line[16];		//holds ID
	char line2[16];	//holds the next ID

	int i, j;

	long int Tally[61];	//index is link type (9_?) (0->60) (61 types)
	for(j=0; j<=60; j++){
		Tally[j]=0;
	}

	int span;
	int filenum;

	FILE* fp;

	unsigned long int totalread=0;
	int filecounter;
	int linecounter;

	int crossnumber;
	int shouldbetwo;
	int index;

	for(span=startingspan; span<=startingspan; span=span+10){
		for(filenum=1; filenum<=totalfiles; filenum++){	//for each link file
			filecounter=0;
			sprintf(filename, "./span%d/knot_lists/MC2SAPsHamL%dM%dspan%drun%dnum%d_LP.txt_knot_new", span, L, M, span, runnum, filenum); //name of knot file
			fp = fopen(filename, "r");	//open "filename" (a knot file)

			if(fp != NULL){
				printf("file %s opened successfully\n", filename);
			}
			else{
				printf("file pointer is pointing to NULL. filename is %s\n", filename);
				break;
			}

			filecounter++;
			linecounter=0;

			fgets(line2, 15, fp);
			while( strcpy(line, line2)!=NULL ){
				if( fgets(line2, 15, fp) != NULL){
					linecounter++;
					totalread++;
//					printf("line %d read\n", linecounter);
					if( atoi(&line[0])==9 ){
						sscanf(line, "%d.%d.%d", &crossnumber, &shouldbetwo, &index);
						if(crossnumber==9){
//							printf("found a crossing number 9, index is %d\n", index);
//							printf("line=%s", line);
//							printf("line2=%s", line2);
							if(linecounter<=10000){
								Tally[index-1]++;
							}
						}
						else{
						  printf("problem. Exitting\n");
						  exit(1);
						}
					}
				}
				else{
					break;
				}
			}
			//			printf("linecounter=%d\n", linecounter);
			fclose(fp);
		}
	}

	int ninecount=0;
	printf("Final Tally:\n");
	for(j=0; j<=60; j++){
		if(Tally[j]>0){
			printf("L=%d, M=%d, Span=%d, 9.2.%d = %li\n", L, M, i*10+10, j+1, Tally[j]);
			ninecount = ninecount + Tally[j];
		}
	}
	printf("ninecount=%d\n", ninecount);
	printf("\n");


	printf("totalread=%lu\n", totalread);
	printf("\nmain ending\n\n");

	return 0;

}






