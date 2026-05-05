#include "mc_globals.h"

void
enterhinge(PolygonState *poly, int i, int j, int side, int (*pointordNum)[3], int curlength)
{
	if( !(poly->alreadyentered[i][j]==1 && side==0) ){
	//	printf("entered hinge (i=%d, j=%d, side=%d\n", i, j, side);

		curlength++;
		(*poly->num_walks)++;
		if(side==0){
			poly->curstart[0][(*poly->num_walks)-1]=0;
			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 1;
		}
		else{ //side==1
			poly->curstart[0][(*poly->num_walks)-1]=2;
			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 2;
		}
		poly->curstart[1][(*poly->num_walks)-1]=j;
		poly->curstart[2][(*poly->num_walks)-1]=i;


		poly->hingestatus[i][j]=1;
		poly->ordertemplate[side][i][j] = (*pointordNum)[side];
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

		leavehinge(poly, i, j, side ^ 1, pointordNum, curlength);
			/* the SAW leaves the hinge out the opposite side, hence (side ^ 1) is passed */
			/* if side = 1 then side ^ 1 = 0 */
			/* if side = 0 then side ^ 1 = 1 */
		rowedges(poly, i, j, pointordNum, curlength);
			/* horizontal edges are explored */
		coledges(poly, i, j, pointordNum, curlength);
			/* vertical edges are explored */

		(*pointordNum)[side]--;
			/* the edge number is decremented since, control is being passed back to  */
			/* the function that called enterhinge. That is to say, the program will  */
			/* now explore alternate ways of entering the hinge (i.e. from different  */
			/* locations). In essense the intial incremention above is nullified so	  */
			/* when enterhinge is called again the entering edge is numbered properly */
		
		poly->ordertemplate[side][i][j] = 0;
			/* the edge that was recorded must now be erased as this edge may not exist */
			/* as part of the SAW when alternate possibilities are explored. In the case */
			/* that it does not exist it needs to be numbered as "0" (zero).             */
		poly->hingestatus[i][j]=0;

		poly->curstart[0][(*poly->num_walks)-1]=-1;
		poly->curstart[1][(*poly->num_walks)-1]=-1;
		poly->curstart[2][(*poly->num_walks)-1]=-1;
		poly->curwalks[curlength-1][(*poly->num_walks)-1] = 0;

		curlength--;
		(*poly->num_walks)--;
	//	printf("finished enterhinge\n");
	}
	return;

}				/* end of function enterhinge */

/***************************************************************************/

void
leavehinge(PolygonState *poly, int i, int j, int side, int (*pointordNum)[3], int curlength)
/* called by: rowedges,coledges,enterhinge */
{
	if( !(poly->alreadyentered[i][j]==1 && side==0) ){
	//	printf("left hinge (i=%d, j=%d, side=%d\n", i, j, side);

		curlength++;
		if(side==0){
			poly->curend[0][(*poly->num_walks)-1]=0;
			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 2;
		}
		else{ //side==1
			poly->curend[0][(*poly->num_walks)-1]=2;
			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 1;
		}
		poly->curend[1][(*poly->num_walks)-1]=j;
		poly->curend[2][(*poly->num_walks)-1]=i;


		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		poly->ordertemplate[side][i][j] = (*pointordNum)[side];	/* record the section
									 * edge (i.e. the
									 * leaving edge) */
		(*pointordNum)[side]++;
			/* the edge number is incremented so that the next edge gets a higher number */

		if (side == 0 && (*pointordNum)[1] > 1) {
	 	/* if the leaving edge is on the left and there is at least one edge on the right section */
		/* this may be a valid two span */
			if (LFlag(pointordNum) /*&& RFlag(pointordNum)*/){
				int isHam = 1;
				if (ham_check) {
					for (int ii = 0; ii <= M; ii++) {
						for (int jj = 0; jj <= L; jj++) {
							if (poly->hingestatus[ii][jj] == 0) { isHam = 0; break; }
						}
						if (!isHam) break;
					}
				}

				if (isHam) {
					fillreordertemplate((*pointordNum)[0]-1, (*pointordNum)[1]-1);
					valid_2_spans++; /*This is a valid 2 span so it is counted. Includes duplicates. */
					recordtemplate(pointordNum); /*the information contained in ordertemplate is recorded*/
				}
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
				if (poly->hingestatus[ii][jj] == 0) {
				/* if the vertex is not occupied then explore entering the hinge at this vertex */
				/* this maintains the self avoiding nature of the walk */
					enterhinge(poly, ii, jj, side, pointordNum,0);
				}
			}
		}
		(*pointordNum)[side]--;
			/* the edge number is decremented since, control is being passed back to  */
			/* the function that called leavehinge. That is to say, the program will  */
			/* now explore alternate ways of leaving the hinge (i.e. via different    */
			/* paths. this may be as simple as going out the other side). In essense  */
			/* the intial incremention above is nullified so when leavehinge is 	  */
			/* called again the leaving edge is numbered properly  			  */
		
		poly->ordertemplate[side][i][j] = 0;
			/* the edge that was recorded must now be erased as this edge may not exist  */
			/* as part of the SAW when alternate possibilities are explored. In the case */
			/* that it does not exist it needs to be numbered as "0" (zero).             */
	//		printf("finished leavehinge\n");


		poly->curend[0][(*poly->num_walks)-1]=-1;
		poly->curend[1][(*poly->num_walks)-1]=-1;
		poly->curend[2][(*poly->num_walks)-1]=-1;
		poly->curwalks[curlength-1][(*poly->num_walks)-1] = 0;
		curlength--;
	}
	return;
}
/***************************************************************************/

void rowedges(PolygonState *poly, int i, int j, int (*pointordNum)[3], int curlength)  //y-direction
{
//	printf("rowedges called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(poly->hingestatus[i][j-1]==0){
			poly->hingestatus[i][j-1]=1;

			curlength++;
			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 4;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			poly->rowhingeedges[i][j-1]=1;
			leavehinge(poly, i, j-1, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(poly, i, j-1, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			rowedges(poly, i, j-1, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			coledges(poly, i, j-1, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			poly->hingestatus[i][j-1]=0;
			poly->rowhingeedges[i][j-1]=0;

			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 0;
			curlength--;
		}
	}
	if(j<L){
		if(poly->hingestatus[i][j+1]==0){
			poly->hingestatus[i][j+1]=1;

			curlength++;
			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 3;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			poly->rowhingeedges[i][j]=1;
			leavehinge(poly, i, j+1, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(poly, i, j+1, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			rowedges(poly, i, j+1, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			coledges(poly, i, j+1, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			poly->hingestatus[i][j+1]=0;
			poly->rowhingeedges[i][j]=0;

			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 0;
			curlength--;
		}
	}
//	printf("finished rowedges\n");
	return;
}				/* end of function rowedges */

/***************************************************************************/

void coledges(PolygonState *poly, int i, int j, int (*pointordNum)[3], int curlength)  //z-direction
/* called by rowedges,coledges */
{
//	printf("coledges called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(poly->hingestatus[i-1][j]==0){
			poly->hingestatus[i-1][j]=1;

			curlength++;
			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 6;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			poly->colhingeedges[i-1][j]=1;
			leavehinge(poly, i-1, j, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(poly, i-1, j, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			rowedges(poly, i-1, j, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			coledges(poly, i-1, j, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			poly->hingestatus[i-1][j]=0;
			poly->colhingeedges[i-1][j]=0;

			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 0;
			curlength--;
		}
	}
	if(i<M){
		if(poly->hingestatus[i+1][j]==0){
			poly->hingestatus[i+1][j]=1;

			curlength++;
			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 5;

			(*pointordNum)[2]++;	/* count the edge in the hinge */
						/* if there is ever a need to record the actual configuration of */
						/* the hinge this would be an appropriate place to do so */
			poly->colhingeedges[i][j]=1;
			leavehinge(poly, i+1, j, 0, pointordNum, curlength);
				/* exit the hinge on side 0 with j moved*/
			leavehinge(poly, i+1, j, 1, pointordNum, curlength);
				/* exit the hinge on side 1 with j moved*/
			rowedges(poly, i+1, j, pointordNum, curlength);
				/* explore horizontal moves within the hinge */
			coledges(poly, i+1, j, pointordNum, curlength);
				/* explore vertical moves within the hinge */
			(*pointordNum)[2]--;
				/* decrement the number of edges in the hinge */
			poly->hingestatus[i+1][j]=0;
			poly->colhingeedges[i][j]=0;

			poly->curwalks[curlength-1][(*poly->num_walks)-1] = 0;
			curlength--;
		}
	}
//	printf("finished coledges\n");
	return;
}				/* end of function coledges */

void
enterendhinge(int i, int j, int curlength)
{
	if( alreadyentered[i][j]!=1 ){
	//	printf("entered endhinge (i=%d, j=%d\n", i, j);

		curlength++;
		num_walks++;
		curstart[0][num_walks-1]=0;
		curstart[1][num_walks-1]=j;
		curstart[2][num_walks-1]=i;

		curwalks[curlength-1][num_walks-1] = 1;

		hingestatus[i][j]=1;
		endtemplate[i][j] = EndOrdNum[0];
		EndOrdNum[0]++;


		endrowedges(i, j, curlength);
			/* horizontal edges are explored */
		endcoledges(i, j, curlength);
			/* vertical edges are explored */

		//reset stuff
		EndOrdNum[0]--;
		endtemplate[i][j] = 0;
		hingestatus[i][j]=0;

		curstart[0][num_walks-1]=-1;
		curstart[1][num_walks-1]=-1;
		curstart[2][num_walks-1]=-1;
		curwalks[curlength-1][num_walks-1] = 0;

		curlength--;
		num_walks--;
	//	printf("finished enterhinge\n");
	}
	return;

}				/* end of function enterhinge */

/***************************************************************************/

void
leaveendhinge(int i, int j, int curlength)
/* called by: rowedges,coledges,enterhinge */
{
	if( alreadyentered[i][j]!=1 ){
	//	printf("left hinge (i=%d, j=%d\n", i, j);

		curlength++;
		curend[0][num_walks-1]=0;
		curend[1][num_walks-1]=j;
		curend[2][num_walks-1]=i;
		curwalks[curlength-1][num_walks-1] = 2;

		int ii;	/*vertex number in vertical direction*/
		int jj;		/*vertex number in horizontal direction */

		endtemplate[i][j] = EndOrdNum[0];
		EndOrdNum[0]++;

		if (LFlag_endhinge(&EndOrdNum)){
			int isHam = 1;
			if (ham_check) {
				for (int ii = 0; ii <= M; ii++) {
					for (int jj = 0; jj <= L; jj++) {
						if (hingestatus[ii][jj] == 0) { isHam = 0; break; }
					}
					if (!isHam) break;
				}
			}

			if (isHam) {
				fillrendtemplate();
				recordendtemplate(); /*the information contained in endtemplate is recorded*/
			}
		}




		for (ii = 0; ii <= M; ii++) {
		/* these are the "vertical" locations of the vertices being considered */
			for (jj = 0; jj <= L; jj++) {
			/* these are the horizontal locations of the vertices being considered */
				if (hingestatus[ii][jj] == 0) {
				/* if the vertex is not occupied then explore entering the hinge at this vertex */
				/* this maintains the self avoiding nature of the walk */
					enterendhinge(ii, jj, 0);
				}
			}
		}

		EndOrdNum[0]--;
		endtemplate[i][j] = 0;

		curend[0][num_walks-1]=-1;
		curend[1][num_walks-1]=-1;
		curend[2][num_walks-1]=-1;
		curwalks[curlength-1][num_walks-1] = 0;
		curlength--;
	}
	return;
}
/***************************************************************************/

void endrowedges(int i, int j, int curlength)  //y-direction
{
//	printf("endrowedges called (i=%d, j=%d\n", i, j);
	if(j>0){
		if(hingestatus[i][j-1]==0){
			hingestatus[i][j-1]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 4;

			EndOrdNum[1]++;
			rowhingeedges[i][j-1]=1;

			leaveendhinge(i, j-1, curlength);
			endrowedges(i, j-1, curlength);
			endcoledges(i, j-1, curlength);

			rowhingeedges[i][j-1]=0;
			EndOrdNum[1]--;
			hingestatus[i][j-1]=0;
			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
	if(j<L){
		if(hingestatus[i][j+1]==0){
			hingestatus[i][j+1]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 3;

			EndOrdNum[1]++;
			rowhingeedges[i][j]=1;

			leaveendhinge(i, j+1, curlength);
			endrowedges(i, j+1, curlength);
			endcoledges(i, j+1, curlength);

			rowhingeedges[i][j]=0;
			EndOrdNum[1]--;
			hingestatus[i][j+1]=0;
			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
//	printf("finished endrowedges\n");
	return;
}				/* end of function rowedges */

/***************************************************************************/

void endcoledges(int i, int j, int curlength)  //z-direction
{
//	printf("endcoledges called (i=%d, j=%d\n", i, j);
	if(i>0){
		if(hingestatus[i-1][j]==0){
			hingestatus[i-1][j]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 6;

			EndOrdNum[1]++;
			colhingeedges[i-1][j]=1;

			leaveendhinge(i-1, j, curlength);
			endrowedges(i-1, j, curlength);
			endcoledges(i-1, j, curlength);

			colhingeedges[i-1][j]=0;
			EndOrdNum[1]--;
			hingestatus[i-1][j]=0;
			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
	if(i<M){
		if(hingestatus[i+1][j]==0){
			hingestatus[i+1][j]=1;

			curlength++;
			curwalks[curlength-1][num_walks-1] = 5;

			EndOrdNum[1]++;
			colhingeedges[i][j]=1;

			leaveendhinge(i+1, j, curlength);
			endrowedges(i+1, j, curlength);
			endcoledges(i+1, j, curlength);

			colhingeedges[i][j]=0;
			EndOrdNum[1]--;
			hingestatus[i+1][j]=0;
			curwalks[curlength-1][num_walks-1] = 0;
			curlength--;
		}
	}
//	printf("finished endcoledges\n");
	return;
}				/* end of function coledges */
