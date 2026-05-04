/*
This program will count the total counts of knot types in a series of knotID files.
*/


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
#define	runnum 1
#define	tot_num_files 100

#define	vM (M+1)				/* number of vertices in vertical direction   = M + 1      */
#define	vL (L+1)				/* number of vertices in horizontal direction = L + 1      */


/** END OF DEFINITIONS THAT CHANGE DEPENDING ON SYSTEM BEING STUDIED **/
int main(void) {
	system("date");  /* prints the date and time */
	clock();	/* Start clock to time program */

	int totalfiles=tot_num_files;
	int filenum=1;
	char filename[60];	//this will hold the filename of the current text file of knottypes that will be counted
	long long int counter=0;

	int filecounter=0;
	int linecounter=0;
	int filesread=0;

	for(filenum=1; filenum<=totalfiles; filenum++){	//for each knot file
		filecounter=0;
		sprintf(filename, "MCpolysHamL%dM%dspan%drun%dnum%d_LP.txt-id", L, M, totalspan, runnum, filenum); //name of knot file
		FILE* fp = fopen(filename, "r");	//open "filename" (a knot file)
		if(fp != NULL){
			printf("file %s opened successfully\n", filename);
		}
		else{
			printf("file pointer is pointing to NULL. num=%d\n", filenum);
		  	break;
		}

		filesread++;
		linecounter=0;
		while( fgets(line, 15, fp)!=NULL ){	//while there's a line to read
			linecounter++;
			if( strcmp(line, "0.2.1\n")==0 || strcmp(line, "0.2.1a\n")==0){	//we have an unlink
				unlink++;
				prevtype=1;
			}
			else if ( strcmp(line, "2.2.1a\n")==0 ){	//we have a positive hopf
				poshopf++;
				prevtype=2;
			}
			else if ( strcmp(line, "2.2.1b\n")==0 ){	//we have a negative hopf
			  neghopf++;
			  prevtype=3;
			}
			else if ( strcmp(line, "4.2.1d\n")==0 ){
			  fouroned++;
			  prevtype=4;
			}
		        else if( strcmp(line, "4.2.1b\n")==0 ){
			  fouroneb++;
			  prevtype=5;
			}
		        else if( strcmp(line, "4.2.1c\n")==0 ){
			  fouronec++;
			  prevtype=6;
			}
		        else if ( strcmp(line, "4.2.1a\n")==0 ){
			  fouronea++;
			  prevtype=7;
			}
		        else if ( strcmp(line, "5.2.1a\n")==0 ){
			  fiveonea++;
			  prevtype=8;
			}
		        else if ( strcmp(line, "5.2.1b\n")==0 ){
			  fiveoneb++;
			  prevtype=9;
			}
		        else if ( strcmp(line, "5.2.1c\n")==0 ){
			  fiveonec++;
			  prevtype=10;
			}
		        else if ( strcmp(line, "5.2.1d\n")==0 ){
			  fiveoned++;
			  prevtype=11;
			}
		        else if ( strcmp(line, "6.2.1a\n")==0 ){
			  sixonea++;
			  prevtype=12;
			}
		        else if ( strcmp(line, "6.2.1b\n")==0 ){
			  sixoneb++;
			  prevtype=13;
			}
		        else if ( strcmp(line, "6.2.1c\n")==0 ){
			  sixonec++;
			  prevtype=14;
			}
		        else if ( strcmp(line, "6.2.1d\n")==0 ){
			  sixoned++;
			  prevtype=15;
			}
		        else if ( strcmp(line, "6.2.3a\n")==0 ){
			  sixthreea++;
			  prevtype=16;
			}
		        else if ( strcmp(line, "6.2.3b\n")==0 ){
			  sixthreeb++;
			  prevtype=17;
			}
		        else if ( strcmp(line, "6.2.3c\n")==0 ){
			  sixthreec++;
			  prevtype=18;
			}
		        else if ( strcmp(line, "6.2.3d\n")==0 ){
			  sixthreed++;
			  prevtype=19;
			}		  
		        else if ( strcmp(line, "l3.1a\n")==0 ){
			  l3_1a++;
			  prevtype=20;
			}
		        else if ( strcmp(line, "l3.1b\n")==0 ){
			  l3_1b++;
		          prevtype=21;
			}
		        else if ( strcmp(line, "l3.1c\n")==0 ){
			  l3_1c++;
			  prevtype=22;
			}
		        else if ( strcmp(line, "l3.1d\n")==0 ){
			  l3_1d++;
			  prevtype=23;
			}
		        else if ( strcmp(line, "u3.1a\n")==0 ){
			  u3_1a++;
			  prevtype=24;
			}
		        else if ( strcmp(line, "u3.1b\n")==0 ){
			  u3_1b++;
			  prevtype=25;
			}
		        else if ( strcmp(line, "u3.1c\n")==0 ){
			  u3_1c++;
			  prevtype=26;
			}
		        else if ( strcmp(line, "u3.1d\n")==0 ){
			  u3_1d++;
			  prevtype=27;
			}
		        else if ( strcmp(line, "unknown\n")==0 ){
			  unknown++;
			  prevtype=28;
			}
		        else if ( strcmp(line, "6.2.2a\n")==0 ){
			  sixtwoa++;
			  prevtype=29;
			}
		        else if ( strcmp(line, "6.2.2b\n")==0 ){
			  sixtwob++;
			  prevtype=30;
			}
		        else if ( strcmp(line, "7.2.1a\n")==0 ){
			  sevenonea++;
			  prevtype=31;
			}
		        else if ( strcmp(line, "7.2.1b\n")==0 ){
			  sevenoneb++;
			  prevtype=32;
			}
		        else if ( strcmp(line, "7.2.1c\n")==0 ){
			  sevenonec++;
			  prevtype=33;
			}
		        else if ( strcmp(line, "7.2.1d\n")==0 ){
			  sevenoned++;
			  prevtype=34;
			}
		        else if ( strcmp(line, "7.2.2a\n")==0 ){
			  seventwoa++;
			  prevtype=35;
			}
		        else if ( strcmp(line, "7.2.2b\n")==0 ){
			  seventwob++;
			  prevtype=36;
			}
		        else if ( strcmp(line, "7.2.2c\n")==0 ){
			  seventwoc++;
			  prevtype=37;
			}
		        else if ( strcmp(line, "7.2.2d\n")==0 ){
			  seventwod++;
			  prevtype=38;
			}
		        else if ( strcmp(line, "7.2.3a\n")==0 ){
			  seventhreea++;
			  prevtype=39;
			}
		        else if ( strcmp(line, "7.2.3b\n")==0 ){
			  seventhreeb++;
			  prevtype=40;
			}
		        else if ( strcmp(line, "7.2.3c\n")==0 ){
			  seventhreec++;
			  prevtype=41;
			}
		        else if ( strcmp(line, "7.2.3d\n")==0 ){
			  seventhreed++;
			  prevtype=42;
			}
		        else if ( strcmp(line, "7.2.4a\n")==0 ){
			  sevenfoura++;
			  prevtype=43;
			}
		        else if ( strcmp(line, "7.2.4b\n")==0 ){
			  sevenfourb++;
			  prevtype=44;
			}
		        else if ( strcmp(line, "7.2.4c\n")==0 ){
			  sevenfourc++;
			  prevtype=45;
			}
		        else if ( strcmp(line, "7.2.4d\n")==0 ){
			  sevenfourd++;
			  prevtype=46;
			}
		        else if ( strcmp(line, "7.2.5a\n")==0 ){
			  sevenfivea++;
			  prevtype=47;
			}
		        else if ( strcmp(line, "7.2.5b\n")==0 ){
			  sevenfiveb++;
			  prevtype=48;
			}
		        else if ( strcmp(line, "7.2.5c\n")==0 ){
			  sevenfivec++;
			  prevtype=49;
			}
		        else if ( strcmp(line, "7.2.5d\n")==0 ){
			  sevenfived++;
			  prevtype=50;
			}
		        else if ( strcmp(line, "7.2.6a\n")==0 ){
			  sevensixa++;
			  prevtype=51;
			}
		        else if ( strcmp(line, "7.2.6b\n")==0 ){
			  sevensixb++;
			  prevtype=52;
			}
		        else if ( strcmp(line, "7.2.6c\n")==0 ){
			  sevensixc++;
			  prevtype=53;
			}
		        else if ( strcmp(line, "7.2.6d\n")==0 ){
			  sevensixd++;
			  prevtype=54;
			}
		        else if ( strcmp(line, "7.2.7a\n")==0 ){
			  sevensevena++;
			  prevtype=55;
			}
		        else if ( strcmp(line, "7.2.7b\n")==0 ){
			  sevensevenb++;
			  prevtype=56;
			}
		        else if ( strcmp(line, "7.2.7c\n")==0 ){
			  sevensevenc++;
			  prevtype=57;
			}
		        else if ( strcmp(line, "7.2.7d\n")==0 ){
			  sevensevend++;
			  prevtype=58;
			}
		        else if ( strcmp(line, "7.2.8a\n")==0 ){
			  seveneighta++;
			  prevtype=59;
			}
		        else if ( strcmp(line, "7.2.8b\n")==0 ){
			  seveneightb++;
			  prevtype=60;
			}
		        else if ( strcmp(line, "7.2.8c\n")==0 ){
			  seveneightc++;
			  prevtype=61;
			}
		        else if ( strcmp(line, "7.2.8d\n")==0 ){
			  seveneightd++;
			  prevtype=62;
			}
		        else if ( strcmp(line, "8.2.1a\n")==0 ){
			  eightone_a++;
			  prevtype=63;
			}
		        else if ( strcmp(line, "8.2.1b\n")==0 ){
			  eightone_b++;
			  prevtype=64;
			}
		        else if ( strcmp(line, "8.2.1c\n")==0 ){
			  eightone_c++;
			  prevtype=65;
			}
		        else if ( strcmp(line, "8.2.1d\n")==0 ){
			  eightone_d++;
			  prevtype=66;
			}
		        else if ( strcmp(line, "8.2.2a\n")==0 ){
			  eighttwo_a++;
			  prevtype=67;
			}
		        else if ( strcmp(line, "8.2.2b\n")==0 ){
			  eighttwo_b++;
			  prevtype=68;
			}
		        else if ( strcmp(line, "8.2.2c\n")==0 ){
			  eighttwo_c++;
			  prevtype=69;
			}
		        else if ( strcmp(line, "8.2.2d\n")==0 ){
			  eighttwo_d++;
			  prevtype=70;
			}
		        else if ( strcmp(line, "8.2.3a\n")==0 ){
			  eightthree_a++;
			  prevtype=71;
			}
		        else if ( strcmp(line, "8.2.3b\n")==0 ){
			  eightthree_b++;
			  prevtype=72;
			}
		        else if ( strcmp(line, "8.2.3c\n")==0 ){
			  eightthree_c++;
			  prevtype=73;
			}
		        else if ( strcmp(line, "8.2.3d\n")==0 ){
			  eightthree_d++;
			  prevtype=74;
			}
		        else if ( strcmp(line, "8.2.4a\n")==0 ){
			  eightfour_a++;
			  prevtype=75;
			}
		        else if ( strcmp(line, "8.2.4b\n")==0 ){
			  eightfour_b++;
			  prevtype=76;
			}
		        else if ( strcmp(line, "8.2.4c\n")==0 ){
			  eightfour_c++;
			  prevtype=77;
			}
		        else if ( strcmp(line, "8.2.4d\n")==0 ){
			  eightfour_d++;
			  prevtype=78;
			}
		        else if ( strcmp(line, "8.2.5a\n")==0 ){
			  eightfive_a++;
			  prevtype=79;
			}
		        else if ( strcmp(line, "8.2.5b\n")==0 ){
			  eightfive_b++;
			  prevtype=80;
			}
		        else if ( strcmp(line, "8.2.5c\n")==0 ){
			  eightfive_c++;
			  prevtype=81;
			}
		        else if ( strcmp(line, "8.2.5d\n")==0 ){
			  eightfive_d++;
			  prevtype=82;
			}
		        else if ( strcmp(line, "8.2.6a\n")==0 ){
			  eightsix_a++;
			  prevtype=83;
			}
		        else if ( strcmp(line, "8.2.6b\n")==0 ){
			  eightsix_b++;
			  prevtype=84;
			}
		        else if ( strcmp(line, "8.2.6c\n")==0 ){
			  eightsix_c++;
			  prevtype=85;
			}
		        else if ( strcmp(line, "8.2.6d\n")==0 ){
			  eightsix_d++;
			  prevtype=86;
			}
		        else if ( strcmp(line, "8.2.7a\n")==0 ){
			  eightseven_a++;
			  prevtype=87;
			}
		        else if ( strcmp(line, "8.2.7b\n")==0 ){
			  eightseven_b++;
			  prevtype=88;
			}
		        else if ( strcmp(line, "8.2.7c\n")==0 ){
			  eightseven_c++;
			  prevtype=89;
			}
		        else if ( strcmp(line, "8.2.7d\n")==0 ){
			  eightseven_d++;
			  prevtype=90;
			}
		        else if ( strcmp(line, "8.2.8a\n")==0 ){
			  eighteight_a++;
			  prevtype=91;
			}
		        else if ( strcmp(line, "8.2.8b\n")==0 ){
			  eighteight_b++;
			  prevtype=92;
			}
		        else if ( strcmp(line, "8.2.8c\n")==0 ){
			  eighteight_c++;
			  prevtype=93;
			}
		        else if ( strcmp(line, "8.2.8d\n")==0 ){
			  eighteight_d++;
			  prevtype=94;
			}
		        else if ( strcmp(line, "8.2.9a\n")==0 ){
			  eightnine_a++;
			  prevtype=95;
			}
		        else if ( strcmp(line, "8.2.9b\n")==0 ){
			  eightnine_b++;
			  prevtype=96;
			}
		        else if ( strcmp(line, "8.2.9c\n")==0 ){
			  eightnine_c++;
			  prevtype=97;
			}
		        else if ( strcmp(line, "8.2.9d\n")==0 ){
			  eightnine_d++;
			  prevtype=98;
			}
		        else if ( strcmp(line, "8.2.10a\n")==0 ){
			  eightten_a++;
			  prevtype=99;
			}
		        else if ( strcmp(line, "8.2.10b\n")==0 ){
			  eightten_b++;
			  prevtype=100;
			}
		        else if ( strcmp(line, "8.2.10c\n")==0 ){
			  eightten_c++;
			  prevtype=101;
			}
		        else if ( strcmp(line, "8.2.10d\n")==0 ){
			  eightten_d++;
			  prevtype=102;
			}
		        else if ( strcmp(line, "8.2.11a\n")==0 ){
			  eighteleven_a++;
			  prevtype=103;
			}
		        else if ( strcmp(line, "8.2.11b\n")==0 ){
			  eighteleven_b++;
			  prevtype=104;
			}
		        else if ( strcmp(line, "8.2.11c\n")==0 ){
			  eighteleven_c++;
			  prevtype=105;
			}
		        else if ( strcmp(line, "8.2.11d\n")==0 ){
			  eighteleven_d++;
			  prevtype=106;
			}



			else{	//newtype link type
				printf("newtype: %s", line);
//				printf("found on line %d\n", linecounter);
//				scanf("%s", dummy);
				newtype++;
				prevtype=-1;
			}
			counter++;
			filecounter++;
		}
		if(prevtype==1){
		  unlink--;
		}
		else if(prevtype==2){
		  poshopf--;
		}
		else if(prevtype==3){
		  neghopf--;
		}
		else if(prevtype==4){
		  fouroned--;
		}
		else if(prevtype==5){
		  fouroneb--;
		}
		else if(prevtype==6){
		  fouronec--;
		}
		else if (prevtype==7){
		  fouronea--;
		}
		else if (prevtype==8){
		  fiveonea--;
		}
		else if (prevtype==9){
		  fiveoneb--;
		}
		else if (prevtype==10){
		  fiveonec--;
		}
		else if (prevtype==11){
		  fiveoned--;
		}
		else if (prevtype==12){
		  sixonea--;
		}
		else if (prevtype==13){
		  sixoneb--;
		}
		else if (prevtype==14){
		  sixonec--;
		}
		else if (prevtype==15){
		  sixoned--;
		}
		else if (prevtype==16){
		  sixthreea--;
		}
		else if (prevtype==17){
		  sixthreeb--;
		}
		else if (prevtype==18){
		  sixthreec--;
		}
		else if (prevtype==19){
		  sixthreed--;
		}
		else if (prevtype==20){
		  l3_1a--;
		}
		else if (prevtype==21){
		  l3_1b--;
		}
		else if (prevtype==22){
		  l3_1c--;
		}
		else if (prevtype==23){
		  l3_1d--;
		}
		else if (prevtype==24){
		  u3_1a--;
		}
		else if (prevtype==25){
		  u3_1b--;
		}
		else if (prevtype==26){
		  u3_1c--;
		}
		else if (prevtype==27){
		  u3_1d--;
		}
		else if (prevtype==28){
		  unknown--;
		}
		else if (prevtype==29){
		  sixtwoa--;
		}
		else if (prevtype==30){
		  sixtwob--;
		}
		else if (prevtype==31){
		  sevenonea--;
		}
		else if (prevtype==32){
		  sevenoneb--;
		}
		else if (prevtype==33){
		  sevenonec--;
		}
		else if (prevtype==34){
		  sevenoned--;
		}
		else if (prevtype==35){
		  seventwoa--;
		}
		else if (prevtype==36){
		  seventwob--;
		}
		else if (prevtype==37){
		  seventwoc--;
		}
		else if (prevtype==38){
		  seventwod--;
		}
		else if (prevtype==39){
		  seventhreea--;
		}
		else if (prevtype==40){
		  seventhreeb--;
		}
		else if (prevtype==41){
		  seventhreec--;
		}
		else if (prevtype==42){
		  seventhreed--;
		}
		else if (prevtype==43){
		  sevenfoura--;
		}
		else if (prevtype==44){
		  sevenfourb--;
		}
		else if (prevtype==45){
		  sevenfourc--;
		}
		else if (prevtype==46){
		  sevenfourd--;
		}
		else if (prevtype==47){
		  sevenfivea--;
		}
		else if (prevtype==48){
		  sevenfiveb--;
		}
		else if (prevtype==49){
		  sevenfivec--;
		}
		else if (prevtype==50){
		  sevenfived--;
		}
		else if (prevtype==51){
		  sevensixa--;
		}
		else if (prevtype==52){
		  sevensixb--;
		}
		else if (prevtype==53){
		  sevensixc--;
		}
		else if (prevtype==54){
		  sevensixd--;
		}
		else if (prevtype==55){
		  sevensevena--;
		}
		else if (prevtype==56){
		  sevensevenb--;
		}
		else if (prevtype==57){
		  sevensevenc--;
		}
		else if (prevtype==58){
		  sevensevend--;
		}
		else if (prevtype==59){
		  seveneighta--;
		}
		else if (prevtype==60){
		  seveneightb--;
		}
		else if (prevtype==61){
		  seveneightc--;
		}
		else if (prevtype==62){
		  seveneightd--;
		}
		else if (prevtype==63){
		  eightone_a--;
		}
		else if (prevtype==64){
		  eightone_b--;
		}
		else if (prevtype==65){
		  eightone_c--;
		}
		else if (prevtype==66){
		  eightone_d--;
		}
		else if (prevtype==67){
		  eighttwo_a--;
		}
		else if (prevtype==68){
		  eighttwo_b--;
		}
		else if (prevtype==69){
		  eighttwo_c--;
		}
		else if (prevtype==70){
		  eighttwo_d--;
		}
		else if (prevtype==71){
		  eightthree_a--;
		}
		else if (prevtype==72){
		  eightthree_b--;
		}
		else if (prevtype==73){
		  eightthree_c--;
		}
		else if (prevtype==74){
		  eightthree_d--;
		}
		else if (prevtype==75){
		  eightfour_a--;
		}
		else if (prevtype==76){
		  eightfour_b--;
		}
		else if (prevtype==77){
		  eightfour_c--;
		}
		else if (prevtype==78){
		  eightfour_d--;
		}
		else if (prevtype==79){
		  eightfive_a--;
		}
		else if (prevtype==80){
		  eightfive_b--;
		}
		else if (prevtype==81){
		  eightfive_c--;
		}
		else if (prevtype==82){
		  eightfive_d--;
		}
		else if (prevtype==83){
		  eightsix_a--;
		}
		else if (prevtype==84){
		  eightsix_b--;
		}
		else if (prevtype==85){
		  eightsix_c--;
		}
		else if (prevtype==86){
		  eightsix_d--;
		}
		else if (prevtype==87){
		  eightseven_a--;
		}
		else if (prevtype==88){
		  eightseven_b--;
		}
		else if (prevtype==89){
		  eightseven_c--;
		}
		else if (prevtype==90){
		  eightseven_d--;
		}
		else if (prevtype==91){
		  eighteight_a--;
		}
		else if (prevtype==92){
		  eighteight_b--;
		}
		else if (prevtype==93){
		  eighteight_c--;
		}
		else if (prevtype==94){
		  eighteight_d--;
		}
		else if (prevtype==95){
		  eightnine_a--;
		}
		else if (prevtype==96){
		  eightnine_b--;
		}
		else if (prevtype==97){
		  eightnine_c--;
		}
		else if (prevtype==98){
		  eightnine_d--;
		}
		else if (prevtype==99){
		  eightten_a--;
		}
		else if (prevtype==100){
		  eightten_b--;
		}
		else if (prevtype==101){
		  eightten_c--;
		}
		else if (prevtype==102){
		  eightten_d--;
		}
		else if (prevtype==103){
		  eighteleven_a--;
		}
		else if (prevtype==104){
		  eighteleven_b--;
		}
		else if (prevtype==105){
		  eighteleven_c--;
		}
		else if (prevtype==106){
		  eighteleven_d--;
		}

		
		else{
		  printf("last link was newtype\n");
		}
		if(filecounter!=10001){
		  printf("file %d has < 10000 links. it has %d only\n", filenum, filecounter-1);
		}



	}

	printf("0_1: %lld\n", unlink);
	printf("2_1a: %d\n", poshopf);
	printf("2_1b: %d\n", neghopf);
	printf("4_1a: %d\n", fouronea);
	printf("4_1b: %d\n", fouroneb);
	printf("4_1c: %d\n", fouronec);
	printf("4_1d: %d\n", fouroned);
	printf("5_1a: %d\n", fiveonea);
	printf("5_1b: %d\n", fiveoneb);
	printf("5_1c: %d\n", fiveonec);
	printf("5_1d: %d\n", fiveoned);
	printf("6_1a: %d\n", sixonea);
	printf("6_1b: %d\n", sixoneb);
	printf("6_1c: %d\n", sixonec);
	printf("6_1d: %d\n", sixoned);
	printf("6_2a: %d\n", sixtwoa);
	printf("6_2b: %d\n", sixtwob);
	printf("6_3a: %d\n", sixthreea);
	printf("6_3b: %d\n", sixthreeb);
	printf("6_3c: %d\n", sixthreec);
	printf("6_3d: %d\n", sixthreed);

	printf("7_1a: %d\n", sevenonea);
	printf("7_1b: %d\n", sevenoneb);
	printf("7_1c: %d\n", sevenonec);
	printf("7_1d: %d\n", sevenoned);
	printf("7_2a: %d\n", seventwoa);
	printf("7_2b: %d\n", seventwob);
	printf("7_2c: %d\n", seventwoc);
	printf("7_2d: %d\n", seventwod);
	printf("7_3a: %d\n", seventhreea);
	printf("7_3b: %d\n", seventhreeb);
	printf("7_3c: %d\n", seventhreec);
	printf("7_3d: %d\n", seventhreed);
	printf("7_4a: %d\n", sevenfoura);
	printf("7_4b: %d\n", sevenfourb);
	printf("7_4c: %d\n", sevenfourc);
	printf("7_4d: %d\n", sevenfourd);
	printf("7_5a: %d\n", sevenfivea);
	printf("7_5b: %d\n", sevenfiveb);
	printf("7_5c: %d\n", sevenfivec);
	printf("7_5d: %d\n", sevenfived);
	printf("7_6a: %d\n", sevensixa);
	printf("7_6b: %d\n", sevensixb);
	printf("7_6c: %d\n", sevensixc);
	printf("7_6d: %d\n", sevensixd);
	printf("7_7a: %d\n", sevensevena);
	printf("7_7b: %d\n", sevensevenb);
	printf("7_7c: %d\n", sevensevenc);
	printf("7_7d: %d\n", sevensevend);
	printf("7_8a: %d\n", seveneighta);
	printf("7_8b: %d\n", seveneightb);
	printf("7_8c: %d\n", seveneightc);
	printf("7_8d: %d\n", seveneightd);

	printf("8_1a: %d\n", eightone_a);
	printf("8_1b: %d\n", eightone_b);
	printf("8_1c: %d\n", eightone_c);
	printf("8_1d: %d\n", eightone_d);
	printf("8_2a: %d\n", eighttwo_a);
	printf("8_2b: %d\n", eighttwo_b);
	printf("8_2c: %d\n", eighttwo_c);
	printf("8_2d: %d\n", eighttwo_d);
	printf("8_3a: %d\n", eightthree_a);
	printf("8_3b: %d\n", eightthree_b);
	printf("8_3c: %d\n", eightthree_c);
	printf("8_3d: %d\n", eightthree_d);
	printf("8_4a: %d\n", eightfour_a);
	printf("8_4b: %d\n", eightfour_b);
	printf("8_4c: %d\n", eightfour_c);
	printf("8_4d: %d\n", eightfour_d);
	printf("8_5a: %d\n", eightfive_a);
	printf("8_5b: %d\n", eightfive_b);
	printf("8_5c: %d\n", eightfive_c);
	printf("8_5d: %d\n", eightfive_d);
	printf("8_6a: %d\n", eightsix_a);
	printf("8_6b: %d\n", eightsix_b);
	printf("8_6c: %d\n", eightsix_c);
	printf("8_6d: %d\n", eightsix_d);
	printf("8_7a: %d\n", eightseven_a);
	printf("8_7b: %d\n", eightseven_b);
	printf("8_7c: %d\n", eightseven_c);
	printf("8_7d: %d\n", eightseven_d);
	printf("8_8a: %d\n", eighteight_a);
	printf("8_8b: %d\n", eighteight_b);
	printf("8_8c: %d\n", eighteight_c);
	printf("8_8d: %d\n", eighteight_d);
	printf("8_9a: %d\n", eightnine_a);
	printf("8_9b: %d\n", eightnine_b);
	printf("8_9c: %d\n", eightnine_c);
	printf("8_9d: %d\n", eightnine_d);
	printf("8_10a: %d\n", eightten_a);
	printf("8_10b: %d\n", eightten_b);
	printf("8_10c: %d\n", eightten_c);
	printf("8_10d: %d\n", eightten_d);
	printf("8_11a: %d\n", eighteleven_a);
	printf("8_11b: %d\n", eighteleven_b);
	printf("8_11c: %d\n", eighteleven_c);
	printf("8_11d: %d\n", eighteleven_d);

	printf("l3_1a: %d\n", l3_1a);
	printf("l3_1b: %d\n", l3_1b);
	printf("l3_1c: %d\n", l3_1c);
	printf("l3_1d: %d\n", l3_1d);
	printf("u3_1a: %d\n", u3_1a);
	printf("u3_1b: %d\n", u3_1b);
	printf("u3_1c: %d\n", u3_1c);
	printf("u3_1d: %d\n", u3_1d);

	printf("unknown: %d\n", unknown);
	printf("newtype: %d\n", newtype);
       
	//	printf("counter=%d, filecounter=%d\n", counter, filecounter);
	printf("total links read: %lld\n", counter-filesread);

	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);
	system("date");  /* prints the date and time */
	return(1);
}				/* end of main */
/***************************************************************************/


