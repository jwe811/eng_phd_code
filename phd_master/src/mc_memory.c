#include "mc_globals.h"
#include "mc_runtime.h"

struct hinge_span* newhinge(unsigned int numberofwalks)		/* returns a pointer to a hinge_span structure */
/* this is a memory allocation function */
{
	struct hinge_span *nextnewhinge;

	nextnewhinge = (struct hinge_span *)mc_xcalloc(1, sizeof(struct hinge_span), "hinge span");
	int i,j;
	nextnewhinge->hedges = (unsigned int*)mc_xcalloc(vM * vL, sizeof(unsigned int), "hinge edge list");
	nextnewhinge->walks = (int**)mc_xcalloc(vM * vL + 1, sizeof(int*), "hinge walks");
	
	for(i=0; i<=vM*vL-1; i++){
		nextnewhinge->hedges[i]=0;
	}

	//allocate space for start, end, and walks
	for(i=0; i<=2; i++){
		nextnewhinge->start[i] = int_vecalloc(0,numberofwalks-1);
		nextnewhinge->end[i] = int_vecalloc(0,numberofwalks-1);
	}
	for(i=0; i<=vM*vL; i++){
		nextnewhinge->walks[i] = int_vecalloc(0,numberofwalks-1);
	}

	//initialize start, end, and walks to all zeros
	for(i=0; i<=numberofwalks-1; i++){
		for(j=0; j<=2; j++){
			nextnewhinge->start[j][i]=0;
			nextnewhinge->end[j][i]=0;
		}
		for(j=0; j<=vM*vL; j++){
			nextnewhinge->walks[j][i]=0;
		}
	}

	return nextnewhinge;
}


/***************************************************************************/
struct endhinge* newendhinge(unsigned int numberofwalks)		/* returns a pointer to an endhinge structure */
/* this is a memory allocation function */
{
	struct endhinge *nextnewendhinge;

	nextnewendhinge = (struct endhinge *)mc_xcalloc(1, sizeof(struct endhinge), "endhinge");
	int i,j;
	nextnewendhinge->hedges = (unsigned int*)mc_xcalloc(vM * vL, sizeof(unsigned int), "endhinge edge list");
	nextnewendhinge->walks = (int**)mc_xcalloc(vM * vL + 1, sizeof(int*), "endhinge walks");

	for(i=0; i<=vM*vL-1; i++){
		nextnewendhinge->hedges[i]=0;
	}

	//allocate space for start, end, and walks
	for(i=0; i<=2; i++){
		nextnewendhinge->start[i] = int_vecalloc(0,numberofwalks-1);
		nextnewendhinge->end[i] = int_vecalloc(0,numberofwalks-1);
	}
	for(i=0; i<=vM*vL; i++){
		nextnewendhinge->walks[i] = int_vecalloc(0,numberofwalks-1);
	}

	//initialize start, end, and walks to all zeros
	for(i=0; i<=numberofwalks-1; i++){
		for(j=0; j<=2; j++){
			nextnewendhinge->start[j][i]=0;
			nextnewendhinge->end[j][i]=0;
		}
		for(j=0; j<=vM*vL; j++){
			nextnewendhinge->walks[j][i]=0;
		}
	}

	return nextnewendhinge;
}
