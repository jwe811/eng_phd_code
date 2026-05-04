unsigned short int		/* returns 1 if testhinge connects, 0 otherwise */
LFlag2(int (*pointordNum)[6])
{
	unsigned short int i, j, b; //looping variables
	if ((*pointordNum)[0]+(*pointordNum)[3] - 2 < vM*vL){
		//not all vertices are filled
		return 1;
	}
	else { //all vertices are filled.

		if(L==1 && M==1){ //then if they don't criss cross, return 1.
			if(ordertemplate[0][0][0] != 0 && ordertemplate[0][1][1] == 0){	//then they dont crisscross and we're ok.
				return(1);
			}
			else if(ordertemplate2[0][0][0] != 0 && ordertemplate2[0][1][1] == 0){	//then 2nd poly doesnt crisscross and we're ok
				return(1);
			}
		}
		else{

			//check 1st polygon to see if any pairs are right beside each other. There must be >2 edges in first polygon to do this, or closing it off will kill poly
			if((*pointordNum)[0]-1 >2){

				for (i = 0; i <= M; i++) {
					for (j = 0; j < L; j++) { /*the i and j loops serve to loop through*/
								/* all posible pairs of horizontally neighboring edges */
								/* that is i,j and i,j+1  hence the j loop ends at L-1 */

						switch ((ordertemplate[0][i][j]) - (b=ordertemplate[0][i][j + 1])) {
							/* look at the difference between the numbers of the two edges */

						case (1):{	if ( b!= 0 && (b & 1) == 0) /*b is not zero but is congruent to 0 mod 2*/
								return(1); break; }
								/*these conditions are satisfied if 3 neighbors 2 */
								/* or 5 neighbors 4, or 7 neighbors 6, etc... */
								/* note: it doesn't make sense for 1 to neighbor 0 */
								/* so this case is elliminated by  b != 0 in if statement*/

						case (-1):{	if ( b!=1 && (b & 1) == 1) /*b is not one but is congruent to 1 mod 2*/
								return(1); break; }
								/*these conditions are satisfied if 2 neighbors 3 */
								/* or 4 neighbors 5, or 6 neighbors 7, etc... */
								/* note: it doesn't make sense for 0 to neighbor 1 */
								/* so this case is elliminated by b != 1 in if statement*/
						} /*end switch*/

						if(ordertemplate[0][i][j]!=0 && ordertemplate[0][i][j + 1]!=0){
							if(ordertemplate[0][i][j] - ordertemplate[0][i][j + 1]==(*pointordNum)[0]-1-1){
								return(1);
							}
							else if(ordertemplate[0][i][j] - ordertemplate[0][i][j + 1]==1-((*pointordNum)[0]-1)){
								return(1);
							}
						}
					}
				}
				for (i = 0; i < M; i++) {
					for (j = 0; j <= L; j++) {/*the i and j loops serve to loop through*/
								/* all posible pairs of vertically neighboring edges */
								/* that is i,j and i+1,j  hence the i loop ends at M-1 */

						switch ((ordertemplate[0][i][j]) - (b=ordertemplate[0][i+1][j])) {

						case (1):{	if ( b!= 0 && (b & 1) == 0) /*b is not zero but is congruent to 0 mod 2*/
								return(1); break; }
								/*these conditions are satisfied if 3 neighbors 2 */
								/* or 5 neighbors 4, or 7 neighbors 6, etc... */
								/* note: it doesn't make sense for 1 to neighbor 0 */
								/* so this case is elliminated by b != 0 in if statement*/

						case (-1):{	if ( b!=1 && (b & 1) == 1) /*b is not one but is congruent to 1 mod 2*/
								return(1); break; }
								/*these conditions are satisfied if 2 neighbors 3 */
								/* or 4 neighbors 5, or 6 neighbors 7, etc... */
								/* note: it doesn't make sense for 0 to neighbor 1 */
								/* so this case is elliminated by b != 1 in if statement*/
						} /*end switch*/

						if(ordertemplate[0][i][j]!=0 && ordertemplate[0][i+1][j]!=0){
							if(ordertemplate[0][i][j] - ordertemplate[0][i+1][j]==(*pointordNum)[0]-1-1){
								return(1);
							}
							else if(ordertemplate[0][i][j] - ordertemplate[0][i+1][j]==1-((*pointordNum)[0]-1)){
								return(1);
							}
						}
					}
				}
			}
			//check 2nd polygon to see if any pairs are right beside each other. There must be >2 edges in first polygon to do this, or closing it off will kill poly
			if((*pointordNum)[3]-1 >2){

				for (i = 0; i <= M; i++) {
					for (j = 0; j < L; j++) { /*the i and j loops serve to loop through*/
								/* all posible pairs of horizontally neighboring edges */
								/* that is i,j and i,j+1  hence the j loop ends at L-1 */

						switch ((ordertemplate2[0][i][j]) - (b=ordertemplate2[0][i][j + 1])) {
							/* look at the difference between the numbers of the two edges */

						case (1):{	if ( b!= 0 && (b & 1) == 0) /*b is not zero but is congruent to 0 mod 2*/
								return(1); break; }
								/*these conditions are satisfied if 3 neighbors 2 */
								/* or 5 neighbors 4, or 7 neighbors 6, etc... */
								/* note: it doesn't make sense for 1 to neighbor 0 */
								/* so this case is elliminated by  b != 0 in if statement*/

						case (-1):{	if ( b!=1 && (b & 1) == 1) /*b is not one but is congruent to 1 mod 2*/
								return(1); break; }
								/*these conditions are satisfied if 2 neighbors 3 */
								/* or 4 neighbors 5, or 6 neighbors 7, etc... */
								/* note: it doesn't make sense for 0 to neighbor 1 */
								/* so this case is elliminated by b != 1 in if statement*/
						} /*end switch*/

						if(ordertemplate2[0][i][j]!=0 && ordertemplate2[0][i][j + 1]!=0){
							if(ordertemplate2[0][i][j] - ordertemplate2[0][i][j + 1]==(*pointordNum)[3]-1-1){
								return(1);
							}
							else if(ordertemplate2[0][i][j] - ordertemplate2[0][i][j + 1]==1-((*pointordNum)[3]-1)){
								return(1);
							}
						}
					}
				}

				for (i = 0; i < M; i++) {
					for (j = 0; j <= L; j++) {/*the i and j loops serve to loop through*/
								/* all posible pairs of vertically neighboring edges */
								/* that is i,j and i+1,j  hence the i loop ends at M-1 */

						switch ((ordertemplate2[0][i][j]) - (b=ordertemplate2[0][i+1][j])) {

						case (1):{	if ( b!= 0 && (b & 1) == 0) /*b is not zero but is congruent to 0 mod 2*/
								return(1); break; }
								/*these conditions are satisfied if 3 neighbors 2 */
								/* or 5 neighbors 4, or 7 neighbors 6, etc... */
								/* note: it doesn't make sense for 1 to neighbor 0 */
								/* so this case is elliminated by b != 0 in if statement*/

						case (-1):{	if ( b!=1 && (b & 1) == 1) /*b is not one but is congruent to 1 mod 2*/
								return(1); break; }
								/*these conditions are satisfied if 2 neighbors 3 */
								/* or 4 neighbors 5, or 6 neighbors 7, etc... */
								/* note: it doesn't make sense for 0 to neighbor 1 */
								/* so this case is elliminated by b != 1 in if statement*/

						} /*end switch*/
						if(ordertemplate2[0][i][j]!=0 && ordertemplate2[0][i+1][j]!=0){
							if(ordertemplate2[0][i][j] - ordertemplate2[0][i+1][j]==(*pointordNum)[3]-1-1){
								return(1);
							}
							else if(ordertemplate2[0][i][j] - ordertemplate2[0][i+1][j]==1-((*pointordNum)[3]-1)){
								return(1);
							}
						}
					}
				}
			}
		}	
	}
	return (0); /*if after checking all neighboring edges there has been no */
			/*occurance of a close pair such as 2 and 3 */
			/*then the section being considered does not connect to phi */
			/*recall that all edges are filled! see if statement above */
}
