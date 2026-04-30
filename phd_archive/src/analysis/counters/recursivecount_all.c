void recursivegenerate(int curspan, struct hinge_span* tspantoadd, int firstsec, int x, int y){

/*this is the recursive function that will loop through all intermediate 2-spans and count all possible polygons.
Will also have to add the right endhinges at this point.
curspan is the current span of the polygon.
tspantoadd is the 2span that comes next in the polygon. the walk in this 2-span must be added and then all of the following
2-spans must be looped through (using recursion)
*/

	int i, j;
	int nextsection;
	int curx, cury;
	int prevdirec, direcsadded, newdirec;
	int temp, temp2;
	int counter;
	int mem1;
	int mem2;

//	printf("Adding %dth 2-span's edges(section %d). Previously=%d, now: ", curspan, firstsec, curlength);
	//add tspantoadd's edge count
	curlength = curlength + tspantoadd->edge_one;
	curlength = curlength + tspantoadd->edge_two;
	mem1 = tspantoadd->edge_one;
	mem2 = tspantoadd->edge_two;
//	printf("%d\n", curlength);

	//ALGORITHM: Add right-endhinges, if any match, record.
	//if length (before adding right-endhinges) is smaller, then call recursive generate again.

	if(tspantoadd->num_con_endhinges_right>0){	//then add appropriate edges and see if it has the correct amount.
//		printf("there are %d right endhinges that work with this polygon.\n", tspantoadd->num_con_endhinges_right);
		int curendhinge;
		for(curendhinge=1; curendhinge <= tspantoadd->num_con_endhinges_right; curendhinge++){	//for each endhinge that connects
//			printf("Adding last section (section %d) and right endhinge (endhinge %d) edges at span=%d. Previously=%d, now: ", tspantoadd->outorder, tspantoadd->con_endhinges_right[curendhinge-1], curspan+1, curlength);
			curlength = curlength + tspantoadd->edge_three; //last section's edges.
			curlength = curlength + edges_in_endhinge[tspantoadd->con_endhinges_right[curendhinge-1]];
//			printf("%d\n", curlength);
			if(curlength==totaledges){
//				printf("\nCounting a poly with %d edges\n", curlength);
				if(grand==999999){
					millions++;
					grand=0;
				}
				else{
					grand++;
				}
//				printf("grand=%d\n\n", grand);
			}
			curlength = curlength - tspantoadd->edge_three; //last section's edges.
			curlength = curlength - edges_in_endhinge[tspantoadd->con_endhinges_right[curendhinge-1]];
		}
	}
	if(curlength<=totaledges-3){	//try adding another span
		nextsection = tspantoadd->outorder;
		struct hinge_span* nexttspan;
		nexttspan = first_hinge_span[nextsection];	//point nexttspan at first
		while(nexttspan->nexthinge != NULL){		//loop through all of nextsection's outorders.
			nexttspan = nexttspan->nexthinge;
			recursivegenerate(curspan+1, nexttspan, nextsection, x, y);
		}
	}

	//reset for different intermediate 2-span
	curlength = curlength - mem1;
	curlength = curlength - mem2;



}




















