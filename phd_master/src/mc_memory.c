#include "mc_globals.h"

struct hinge_span* newhinge(unsigned int numberofwalks)		/* returns a pointer to a hinge_span structure */
/* this is a memory allocation function */
{
	struct hinge_span *nextnewhinge;

	nextnewhinge = ((struct hinge_span *) calloc(1, sizeof(struct hinge_span)));
	if (nextnewhinge == NULL) {
		fprintf(stderr, "unable to allocate memory");
		exit(1);
	}
	int i,j;
	nextnewhinge->hedges = (unsigned int*)calloc(vM * vL, sizeof(unsigned int));
	nextnewhinge->walks = (int**)calloc(vM * vL + 1, sizeof(int*));
	
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

	nextnewendhinge = ((struct endhinge *) calloc(1, sizeof(struct endhinge)));
	if (nextnewendhinge == NULL) {
		fprintf(stderr, "unable to allocate memory");
		exit(1);
	}
	int i,j;
	nextnewendhinge->hedges = (unsigned int*)calloc(vM * vL, sizeof(unsigned int));
	nextnewendhinge->walks = (int**)calloc(vM * vL + 1, sizeof(int*));

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
