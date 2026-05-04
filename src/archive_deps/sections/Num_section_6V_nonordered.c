/***************************************************************************/
unsigned long int
num_section_norder(int side)
{
	unsigned long int i,j; 			/*used as dummy variables in for loops */
	unsigned short int edge_num;		/*the number of the edge being considered in for loop*/
						/*this is zero if there is no edge there */
	unsigned short int loc=1;		/*location of the edge being considered in for loop */
	unsigned short int test_edge_num;	/*tests for edges with lower number */
	unsigned short int slot_num[7];	/*of the positions available to the edge determined by*/
						/*the array index slot_num[index] records which available*/
						/*position (or slot) the edge took                       */
						/*							*/
						/*example:  consider the slit with edges: */

						/* 4 3 1 2 0 */
						/* at first all positions are available */
						/* # # # # #                            */
						/* 1 takes the third available position */
						/* # # 1 # #				 */

						/* 2 takes the third position available now that edge 1 is in */
						/* note: this is the location of the vertex minus 1 (since 1 < 2) */
						/* # # 1 2 # */

						/* 3 takes the second position available to it */
						/* # 3 1 2 # */

						/* 4 takes the first position available to it */
						/* 4 3 1 2 0 */

	unsigned short int prev_less;		/*the number of edges with numbers which are less than a given edge number*/
						/*in the vertices considered so far					  */
						/*In the above example: this function considers the enters in the following */
						/*manner */
						/*loop through the vertices (ie. in "reading order" left-right top-to-bottom) */
						/* 4 # # # # */
						/* slot_num[4] = loc - prev_less = 1 - 0 = 1 (ie. first available space) */
						/* 4 3 # # # */
						/* slot_num[3] = loc - prev_less = 2 - 0 = 2 (ie. second available space) */
						/* the fact that the first slot is filled is irrelavent since 4 > 3 */
						/* 4 3 1 # # */
						/* slot_num[1] = 3 */
						/* 4 3 1 2 # */
						/* slot_num[2] = loc - prev_less = 4 - 1 = 3  (ie. third available space */
						/* since 1 < 2 , that is there is a previous edge with a lesser number)  */

for (i=1;i <= 6;i++){slot_num[i]=0;} /* This is necessary to initialize the array to contain zero elements */

	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			edge_num=reordertemplate[side][i][j];/*assign number of edge determined by i and j to edge_num */
							   /*this will be zero if there is no edge */

			if (edge_num != 0){ /*if there is an edge at this position */

				prev_less=0; /*previous occurences of edges in section with a lower number than edge_num*/
				test_edge_num=1; /*start looking with test_edge_num = 1 */

				while (test_edge_num < edge_num){ /*while loop only considers edges with lower numbers */
								  /*and only the vertices considered so far are those  */
								  /*looped through by the i and j for loops	       */

					if (slot_num[test_edge_num] != 0){/*this means that of the vertices considered */
								    /*so far there is an occurence of an edge which */
								    /*has a lower number	                     */
						prev_less++;	/*record that there is an occurence of an edge with */
								/*a lower number than edge_num */
					}
					test_edge_num++;	/*increment test_edge_num to consider other edge numbers */
				}
				slot_num[edge_num] = (loc - prev_less); /*record which available slot the edge edge_num took*/
			}
			loc++;
		}
	}

/*The following is explained in the Maple program which was used to come up with the formulas */
/*basically these formulas assign a unique number to the section given the way in which the */
/*edges where filled */

if (slot_num[5] !=0){
	return(185.0*slot_num[1]-232.0+37.0*slot_num[2]+9.0*slot_num[3]+3.0*slot_num[4]+slot_num[5]);
}
if (slot_num[3] !=0){
	return(185.0*slot_num[1]-232.0+37.0*slot_num[2]+9.0*slot_num[3]+3.0*slot_num[4]);
}
	return(185.0*slot_num[1]-221.0+37.0*slot_num[2]);
}/*end of function num_section*/
/***************************************************************************/























 
