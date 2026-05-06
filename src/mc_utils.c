#include "mc_globals.h"
#include "mc_runtime.h"
#include "run_metadata.h"

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
			exit(EXIT_FAILURE);
		}

		t_end[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_end[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(EXIT_FAILURE);
		}

		t_walks[section_num] = (int***)malloc((num_outsections[section_num]+1) * sizeof(int**));
		if (t_walks[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(EXIT_FAILURE);
		}

		for (i = 1; i <= num_outsections[section_num]; i++){
			t_start[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (t_start[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(EXIT_FAILURE);
			}
			t_end[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (t_end[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(EXIT_FAILURE);
			}
			t_walks[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (t_walks[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(EXIT_FAILURE);
			}

			for (j = 0; j < 3; j++){
				t_start[section_num][i][j] = (int*)malloc(t_num_walks[section_num][i] * sizeof(int));
		   		if (t_start[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(EXIT_FAILURE);
				}
				t_end[section_num][i][j] = (int*)malloc(t_num_walks[section_num][i] * sizeof(int));
		   		if (t_end[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(EXIT_FAILURE);
				}
			}
			for (j = 0; j <= vM*vL; j++){
				t_walks[section_num][i][j] = (int*)malloc(t_num_walks[section_num][i] * sizeof(int));
		   		if (t_walks[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(EXIT_FAILURE);
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
			exit(EXIT_FAILURE);
		}

		Lend_end[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_end[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(EXIT_FAILURE);
		}

		Lend_walks[section_num] = (int***)malloc((num_left_endhinges[section_num]+1) * sizeof(int**));
		if (Lend_walks[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(EXIT_FAILURE);
		}
/////////
		Rend_start[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_start[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(EXIT_FAILURE);
		}

		Rend_end[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_end[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(EXIT_FAILURE);
		}

		Rend_walks[section_num] = (int***)malloc((num_right_endhinges[section_num]+1) * sizeof(int**));
		if (Rend_walks[section_num] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(EXIT_FAILURE);
		}
/////////
		for (i = 1; i <= num_left_endhinges[section_num]; i++){
			Lend_start[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Lend_start[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(EXIT_FAILURE);
			}
			Lend_end[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Lend_end[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(EXIT_FAILURE);
			}
			Lend_walks[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (Lend_walks[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(EXIT_FAILURE);
			}

			for (j = 0; j < 3; j++){
				Lend_start[section_num][i][j] = (int*)malloc(Lend_num_walks[section_num][i] * sizeof(int));
		   		if (Lend_start[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(EXIT_FAILURE);
				}
				Lend_end[section_num][i][j] = (int*)malloc(Lend_num_walks[section_num][i] * sizeof(int));
		   		if (Lend_end[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(EXIT_FAILURE);
				}
			}
			for (j = 0; j <= vM*vL; j++){
				Lend_walks[section_num][i][j] = (int*)malloc(Lend_num_walks[section_num][i] * sizeof(int));
		   		if (Lend_walks[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(EXIT_FAILURE);
				}
			}
		}
/////////
		for (i = 1; i <= num_right_endhinges[section_num]; i++){
			Rend_start[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Rend_start[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(EXIT_FAILURE);
			}
			Rend_end[section_num][i] = (int**)malloc(3 * sizeof(int*));
			if (Rend_end[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(EXIT_FAILURE);
			}
			Rend_walks[section_num][i] = (int**)malloc((vM*vL+1) * sizeof(int*));
			if (Rend_walks[section_num][i] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(EXIT_FAILURE);
			}

			for (j = 0; j < 3; j++){
				Rend_start[section_num][i][j] = (int*)malloc(Rend_num_walks[section_num][i] * sizeof(int));
		   		if (Rend_start[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(EXIT_FAILURE);
				}
				Rend_end[section_num][i][j] = (int*)malloc(Rend_num_walks[section_num][i] * sizeof(int));
		   		if (Rend_end[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(EXIT_FAILURE);
				}
			}
			for (j = 0; j <= vM*vL; j++){
				Rend_walks[section_num][i][j] = (int*)malloc(Rend_num_walks[section_num][i] * sizeof(int));
		   		if (Rend_walks[section_num][i][j] == NULL) {
					fprintf(stderr, "Out of memory");
					exit(EXIT_FAILURE);
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

void printtofile(){
	if(filetotal>=maxpolys){
		//need to start new file
		printf("Finished filling file %lu. It contains %lu polys. Creating a new one.\n", filenum, filetotal);
		fprintf(fp, "-999\n");
		fclose(fp);
		filenum++;
		const char *file_prefix = (ham_check) ? "MCpolysHam" : "MCpolys";
		mc_checked_snprintf(filename, sizeof(filename), "%s/%sL%dM%dspan%drun%dnum%lu.txt", output_dir, file_prefix, L, M, totalspan, runnum, filenum);
		fp = fopen(filename, "w");	//create or overwrite the file "filename"

		if(fp != NULL){
			fprintf(fp, "UofS\n");	//first line in file is always "UofS"
			run_metadata_write(filename, "mc_master", "samples_uofs", mode, L, M, totalspan, seednum, runnum, dom_evalue);
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
