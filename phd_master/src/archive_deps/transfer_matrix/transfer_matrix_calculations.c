/* 
 * transfer_matrix_calculations.c
 *
 * Consolidated Master File for Transfer Matrix calculations in LxM tubes.
 * Features: State Generation, Numerical Solving, Statistical Analysis,
 * and Eigenvector File Export for Monte Carlo Sampling.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

/* ==========================================================================
   CONFIGURATION TOGGLES
   ========================================================================== */

#define MODE_STANDARD    0
#define MODE_HAMILTONIAN 1
#define MODE_2SAP        2
#define MODE_2SAP_HAM    3

#ifndef SYSTEM_MODE
#define SYSTEM_MODE      MODE_STANDARD
#endif

#define LAT_L 1
#define LAT_M 1

/* Compatibility macros */
#define L LAT_L
#define M LAT_M
#define force 0
#define SS

#if SYSTEM_MODE == MODE_2SAP || SYSTEM_MODE == MODE_2SAP_HAM
    #define NUM_POLYS 2
#else
    #define NUM_POLYS 1
#endif

#if SYSTEM_MODE == MODE_HAMILTONIAN || SYSTEM_MODE == MODE_2SAP_HAM
    #define HAM_CHECK 1
#else
    #define HAM_CHECK 0
#endif

#define bits_16 unsigned short int
#define vM (LAT_M+1)
#define vL (LAT_L+1)

/* Memory Scaling */
#if LAT_M==1 && LAT_L==1
    #define max_sections 100
    #define max_keynum 100
#elif LAT_M==1 && LAT_L==2
    #define max_sections 1000
    #define max_keynum 5000
#elif LAT_M==1 && LAT_L==3
    #define max_sections 10000
    #define max_keynum 50000
#else
    #define max_sections 200000
    #define max_keynum 500000
#endif

#define vec_length (max_sections+1)
#define MAX_HEDGE ( (LAT_M+1)*LAT_L + LAT_M*(LAT_L+1) + 1 )

typedef double vec_ent;
typedef double mat_ent;

/* ==========================================================================
   STRUCTURES
   ========================================================================== */

struct hinge_span {
    unsigned long int inorder2;
    unsigned long int outorder;
    unsigned long int outorder2;
    unsigned int edgecount;
    unsigned int edgecount2;
    unsigned int hedges[MAX_HEDGE];
    struct hinge_span *nexthinge;
} *first_hinge_span[vec_length], *current_hinge_span[vec_length];

/* ==========================================================================
   GLOBAL VARIABLES
   ========================================================================== */

unsigned int ordertemplate_master[NUM_POLYS][2][vM][vL];
#define ordertemplate ordertemplate_master[0]
#define ordertemplate2 ordertemplate_master[1]

unsigned int reordertemplate_master[NUM_POLYS][2][vM][vL];
#define reordertemplate reordertemplate_master[0]
#define reordertemplate2 reordertemplate_master[1]

unsigned long int num_2_spans = 0;
unsigned long int valid_2_spans = 0;
unsigned short int hingestatus[vM][vL];
unsigned short int alreadyentered[NUM_POLYS][vM][vL];
unsigned short int colhingeedges[vM][vL];
unsigned short int rowhingeedges[vM][vL];
unsigned long int sectionkey[vec_length];
unsigned long int sectionkey2SAP[max_keynum+1][2];
unsigned long int dupcounter = 0;

unsigned long int num_outsections[max_keynum+1];
unsigned long int *tspans_outsection[max_keynum+1];
unsigned long int *tspans_edges[max_keynum+1];

vec_ent *L_Evector[2], *R_Evector[2];
double fval = 0;
double connectivity_inv;
unsigned long int actual_max_states = 0;

/* ==========================================================================
   UTILITY INCLUDES
   ========================================================================== */

#include "../utils/unsgnlong_vecalloc.c"
#include "../utils/vecalloc.c"
#include "../utils/matalloc.c"

#ifdef CS
    #include "../topology/LFlag_0.c"
#else
    #include "../topology/LFlag_norder.c"
    #if NUM_POLYS == 2
    #include "../topology/LFlag_norder2.c"
    #endif
    
    #if vM*vL < 7
        #include "../sections/Num_section_6V_nonordered.c"
        #if NUM_POLYS == 2
        #include "../sections/Num_section_6V_nonordered2.c"
        #endif
    #elif vM*vL < 9
        #include "../sections/Num_section_8V_nonordered.c"
        #if NUM_POLYS == 2
        #include "../sections/Num_section_8V_nonordered2.c"
        #endif
    #elif vM*vL < 11
        #include "../sections/Num_section_10V_nonordered.c"
        #if NUM_POLYS == 2
        #include "../sections/Num_section_10V_nonordered2.c"
        #endif
    #else
        #include "../sections/Num_section_12V_nonordered.c"
        #if NUM_POLYS == 2
        #include "../sections/Num_section_12V_nonordered2.c"
        #endif
    #endif
#endif

#define QM LAT_M
#undef M
#include "../utils/qcksrtII_int.c"
#define M QM
#include "../analysis/statistics/rtflsp.c"

/* ==========================================================================
   CORE ENGINE
   ========================================================================== */

struct hinge_span* newhinge(void) {
    struct hinge_span *n = (struct hinge_span *) calloc(1, sizeof(struct hinge_span));
    if(!n) exit(1);
    return n;
}

void fillreordertemplate(int poly, int ledges, int redges) {
    int i, j, ii, jj, s;
    int firstentry, secondentry, connectingedge;
    for(s=0; s<=1; s++) {
        for(i=0; i<=LAT_M; i++) for(j=0; j<=LAT_L; j++) reordertemplate_master[poly][s][i][j]=0;
        firstentry=1;
        for(i=0; i<=LAT_M; i++) {
            for(j=0; j<=LAT_L; j++) {
                if(ordertemplate_master[poly][s][i][j] % 2 != 0 && reordertemplate_master[poly][s][i][j]==0) {
                    reordertemplate_master[poly][s][i][j] = firstentry;
                    connectingedge = (ordertemplate_master[poly][s][i][j]==1) ? (s==0 ? ledges : redges) : (ordertemplate_master[poly][s][i][j]-1);
                    secondentry = (firstentry==1) ? (s==0 ? ledges : redges) : (firstentry+1);
                    for(ii=0; ii<=LAT_M; ii++) for(jj=0; jj<=LAT_L; jj++) if(ordertemplate_master[poly][s][ii][jj]==connectingedge) {
                        reordertemplate_master[poly][s][ii][jj] = secondentry;
                        firstentry = (firstentry==1) ? 2 : (firstentry+2);
                        break;
                    }
                } else if(ordertemplate_master[poly][s][i][j]>0 && reordertemplate_master[poly][s][i][j]==0) {
                    reordertemplate_master[poly][s][i][j] = firstentry;
                    connectingedge = (s==0 && ordertemplate_master[poly][s][i][j]==ledges) || (s==1 && ordertemplate_master[poly][s][i][j]==redges) ? 1 : (ordertemplate_master[poly][s][i][j]+1);
                    secondentry = (firstentry==1) ? (s==0 ? ledges : redges) : (firstentry+1);
                    for(ii=0; ii<=LAT_M; ii++) for(jj=0; jj<=LAT_L; jj++) if(ordertemplate_master[poly][s][ii][jj]==connectingedge) {
                        reordertemplate_master[poly][s][ii][jj] = secondentry;
                        firstentry = (firstentry==1) ? 2 : (firstentry+2);
                        break;
                    }
                }
            }
        }
    }
}

unsigned long int get_section_num(int poly, int side) {
    #if NUM_POLYS == 2
    if(poly == 1) return num_section_norder2(side);
    #endif
    return num_section_norder(side);
}

void recordtemplate(int (*pON)[NUM_POLYS*3]) {
    unsigned long int inNum = get_section_num(0, 0);
    int i, j, found=0;
    for(i=1; i<vec_length && sectionkey[i]!=0; i++) if(sectionkey[i]==inNum) { found=1; inNum=i; break; }
    if(!found && i<vec_length) { sectionkey[i]=inNum; inNum=i; }

    unsigned long int outNum = get_section_num(0, 1);
    found=0;
    for(i=1; i<vec_length && sectionkey[i]!=0; i++) if(sectionkey[i]==outNum) { found=1; outNum=i; break; }
    if(!found && i<vec_length) { sectionkey[i]=outNum; outNum=i; }

    unsigned long int inNum2=0, outNum2=0;
    #if NUM_POLYS == 2
    inNum2 = get_section_num(1, 0);
    for(i=1; i<vec_length && sectionkey[i]!=0; i++) if(sectionkey[i]==inNum2) { found=1; inNum2=i; break; }
    if(!found && i<vec_length) { sectionkey[i]=inNum2; inNum2=i; }
    outNum2 = get_section_num(1, 1);
    for(i=1; i<vec_length && sectionkey[i]!=0; i++) if(sectionkey[i]==outNum2) { found=1; outNum2=i; break; }
    if(!found && i<vec_length) { sectionkey[i]=outNum2; outNum2=i; }
    #endif

    unsigned int n_edges = (*pON)[0] + (*pON)[2] - 1;
    unsigned int n_edges2 = (NUM_POLYS==2) ? ((*pON)[3] + (*pON)[5] - 1) : 0;
    
    unsigned int temp_hedges[MAX_HEDGE];
    for(i=0; i<MAX_HEDGE; i++) temp_hedges[i]=0;
    int idx=0, edgenum=1;
    for(i=0; i<=LAT_M; i++) for(j=0; j<LAT_L; j++) { if(rowhingeedges[i][j]) temp_hedges[idx]=edgenum; idx++; edgenum++; }
    for(i=0; i<LAT_M; i++) for(j=0; j<=LAT_L; j++) { if(colhingeedges[i][j]) temp_hedges[idx]=edgenum; idx++; edgenum++; }

    int duplicate=0;
    struct hinge_span* test = first_hinge_span[inNum];
    while(test->nexthinge != NULL) {
        test = test->nexthinge;
        if(test->outorder==outNum && test->inorder2==inNum2 && test->outorder2==outNum2) {
            if(test->edgecount==n_edges && test->edgecount2==n_edges2) {
                int match=1;
                for(i=0; i<MAX_HEDGE; i++) if(temp_hedges[i]!=test->hedges[i]) { match=0; break; }
                if(match) { duplicate=1; dupcounter++; break; }
            }
        }
    }
    if(!duplicate) {
        (*current_hinge_span[inNum]).nexthinge = newhinge();
        current_hinge_span[inNum] = (*current_hinge_span[inNum]).nexthinge;
        (*current_hinge_span[inNum]).inorder2 = inNum2;
        (*current_hinge_span[inNum]).outorder = outNum;
        (*current_hinge_span[inNum]).outorder2 = outNum2;
        (*current_hinge_span[inNum]).edgecount = n_edges;
        (*current_hinge_span[inNum]).edgecount2 = n_edges2;
        for(i=0; i<MAX_HEDGE; i++) (*current_hinge_span[inNum]).hedges[i]=temp_hedges[i];
    }
}

/* NAVIGATION */
void leavehinge(int i, int j, int side, int poly, int (*pON)[NUM_POLYS*3]);
void rowedges(int i, int j, int poly, int (*pON)[NUM_POLYS*3]);
void coledges(int i, int j, int poly, int (*pON)[NUM_POLYS*3]);

#if NUM_POLYS == 2
void enterhinge2(int i, int j, int side, int (*pON)[6]);
void leavehinge2(int i, int j, int side, int (*pON)[6]);
void rowedges2(int i, int j, int (*pON)[6]);
void coledges2(int i, int j, int (*pON)[6]);

void start_poly2(int (*pON)[6]) {
    for (int i=0; i<=LAT_M; i++) {
        for (int j=0; j<=LAT_L; j++) {
            if(!hingestatus[i][j]) enterhinge2(i, j, 0, pON);
            alreadyentered[1][i][j]=1;
        }
    }
    for (int i=0; i<=LAT_M; i++) for (int j=0; j<=LAT_L; j++) alreadyentered[1][i][j]=0;
}

void enterhinge2(int i, int j, int side, int (*pON)[6]) {
    if(alreadyentered[1][i][j]==1 && side==0) return;
    hingestatus[i][j]=1; ordertemplate_master[1][side][i][j] = (*pON)[3 + side]; (*pON)[3 + side]++;
    leavehinge2(i, j, side ^ 1, pON); rowedges2(i, j, pON); coledges2(i, j, pON);
    (*pON)[3 + side]--; ordertemplate_master[1][side][i][j] = 0; hingestatus[i][j]=0;
}

void leavehinge2(int i, int j, int side, int (*pON)[6]) {
    if(alreadyentered[1][i][j]==1 && side==0) return;
    ordertemplate_master[1][side][i][j] = (*pON)[3 + side]; (*pON)[3 + side]++;
    if (side == 0 && (*pON)[4] > 1) {
        if (LFlag2(pON)) {
            #if HAM_CHECK
            int isHam=1;
            for (int ii=0; ii<=LAT_M; ii++) for (int jj=0; jj<=LAT_L; jj++) if (!hingestatus[ii][jj]) { isHam=0; break; }
            if(isHam) { fillreordertemplate(0, (*pON)[0]-1, (*pON)[1]-1); fillreordertemplate(1, (*pON)[3]-1, (*pON)[4]-1); recordtemplate(pON); valid_2_spans++; }
            #else
            fillreordertemplate(0, (*pON)[0]-1, (*pON)[1]-1); fillreordertemplate(1, (*pON)[3]-1, (*pON)[4]-1); recordtemplate(pON); valid_2_spans++;
            #endif
        }
    }
    for (int ii=0; ii<=LAT_M; ii++) for (int jj=0; jj<=LAT_L; jj++) if (!hingestatus[ii][jj]) enterhinge2(ii, jj, side, pON);
    (*pON)[3 + side]--; ordertemplate_master[1][side][i][j] = 0;
}

void rowedges2(int i, int j, int (*pON)[6]) {
    if(j>0 && !hingestatus[i][j-1]) { hingestatus[i][j-1]=1; (*pON)[5]++; rowhingeedges[i][j-1]=1; leavehinge2(i, j-1, 0, pON); leavehinge2(i, j-1, 1, pON); rowedges2(i, j-1, pON); coledges2(i, j-1, pON); (*pON)[5]--; hingestatus[i][j-1]=0; rowhingeedges[i][j-1]=0; }
    if(j<LAT_L && !hingestatus[i][j+1]) { hingestatus[i][j+1]=1; (*pON)[5]++; rowhingeedges[i][j]=1; leavehinge2(i, j+1, 0, pON); leavehinge2(i, j+1, 1, pON); rowedges2(i, j+1, pON); coledges2(i, j+1, pON); (*pON)[5]--; hingestatus[i][j+1]=0; rowhingeedges[i][j]=0; }
}

void coledges2(int i, int j, int (*pON)[6]) {
    if(i>0 && !hingestatus[i-1][j]) { hingestatus[i-1][j]=1; (*pON)[5]++; colhingeedges[i-1][j]=1; leavehinge2(i-1, j, 0, pON); leavehinge2(i-1, j, 1, pON); rowedges2(i-1, j, pON); coledges2(i-1, j, pON); (*pON)[5]--; hingestatus[i-1][j]=0; colhingeedges[i-1][j]=0; }
    if(i<LAT_M && !hingestatus[i+1][j]) { hingestatus[i+1][j]=1; (*pON)[5]++; colhingeedges[i][j]=1; leavehinge2(i+1, j, 0, pON); leavehinge2(i+1, j, 1, pON); rowedges2(i+1, j, pON); coledges2(i+1, j, pON); (*pON)[5]--; hingestatus[i+1][j]=0; colhingeedges[i][j]=0; }
}
#endif

void enterhinge(int i, int j, int side, int poly, int (*pON)[NUM_POLYS*3]) {
    if(alreadyentered[poly][i][j]==1 && side==0) return;
    hingestatus[i][j]=1; ordertemplate_master[poly][side][i][j] = (*pON)[poly*3 + side]; (*pON)[poly*3 + side]++;
    leavehinge(i, j, side ^ 1, poly, pON); rowedges(i, j, poly, pON); coledges(i, j, poly, pON);
    (*pON)[poly*3 + side]--; ordertemplate_master[poly][side][i][j] = 0; hingestatus[i][j]=0;
}

void leavehinge(int i, int j, int side, int poly, int (*pON)[NUM_POLYS*3]) {
    if(alreadyentered[poly][i][j]==1 && side==0) return;
    ordertemplate_master[poly][side][i][j] = (*pON)[poly*3 + side]; (*pON)[poly*3 + side]++;
    if (side == 0 && (*pON)[poly*3 + 1] > 1) {
        #if NUM_POLYS == 2
        if ((*pON)[0]-1 <= vM * vL - 2) start_poly2(pON);
        #else
        if (LFlag(pON)) {
            #if HAM_CHECK
            int isHam=1;
            for (int ii=0; ii<=LAT_M; ii++) for (int jj=0; jj<=LAT_L; jj++) if (!hingestatus[ii][jj]) { isHam=0; break; }
            if(isHam) { fillreordertemplate(0, (*pON)[0]-1, (*pON)[1]-1); recordtemplate(pON); valid_2_spans++; }
            #else
            fillreordertemplate(0, (*pON)[0]-1, (*pON)[1]-1); recordtemplate(pON); valid_2_spans++;
            #endif
        }
        #endif
        num_2_spans++;
    }
    for (int ii=0; ii<=LAT_M; ii++) for (int jj=0; jj<=LAT_L; jj++) if (!hingestatus[ii][jj]) enterhinge(ii, jj, side, poly, pON);
    (*pON)[poly*3 + side]--; ordertemplate_master[poly][side][i][j] = 0;
}

void rowedges(int i, int j, int poly, int (*pON)[NUM_POLYS*3]) {
    if(j>0 && !hingestatus[i][j-1]) { hingestatus[i][j-1]=1; (*pON)[poly*3 + 2]++; rowhingeedges[i][j-1]=1; leavehinge(i, j-1, 0, poly, pON); leavehinge(i, j-1, 1, poly, pON); rowedges(i, j-1, poly, pON); coledges(i, j-1, poly, pON); (*pON)[poly*3 + 2]--; hingestatus[i][j-1]=0; rowhingeedges[i][j-1]=0; }
    if(j<LAT_L && !hingestatus[i][j+1]) { hingestatus[i][j+1]=1; (*pON)[poly*3 + 2]++; rowhingeedges[i][j]=1; leavehinge(i, j+1, 0, poly, pON); leavehinge(i, j+1, 1, poly, pON); rowedges(i, j+1, poly, pON); coledges(i, j+1, poly, pON); (*pON)[poly*3 + 2]--; hingestatus[i][j+1]=0; rowhingeedges[i][j]=0; }
}

void coledges(int i, int j, int poly, int (*pON)[NUM_POLYS*3]) {
    if(i>0 && !hingestatus[i-1][j]) { hingestatus[i-1][j]=1; (*pON)[poly*3 + 2]++; colhingeedges[i-1][j]=1; leavehinge(i-1, j, 0, poly, pON); leavehinge(i-1, j, 1, poly, pON); rowedges(i-1, j, poly, pON); coledges(i-1, j, poly, pON); (*pON)[poly*3 + 2]--; hingestatus[i-1][j]=0; colhingeedges[i-1][j]=0; }
    if(i<LAT_M && !hingestatus[i+1][j]) { hingestatus[i+1][j]=1; (*pON)[poly*3 + 2]++; colhingeedges[i][j]=1; leavehinge(i+1, j, 0, poly, pON); leavehinge(i+1, j, 1, poly, pON); rowedges(i+1, j, poly, pON); coledges(i+1, j, poly, pON); (*pON)[poly*3 + 2]--; hingestatus[i+1][j]=0; colhingeedges[i][j]=0; }
}

void conv_to_array(void) {
    unsigned long int i, j, s, key=0, found;
    struct hinge_span *curr, *tofree;
    #if NUM_POLYS == 2
    for(s=1; s<=max_sections; s++) {
        curr = first_hinge_span[s];
        while(curr->nexthinge) {
            curr = curr->nexthinge;
            found=0;
            for(j=1; j<=key && !found; j++) if(sectionkey2SAP[j][0]==s && sectionkey2SAP[j][1]==curr->inorder2) found=1;
            if(!found && key < max_keynum) { key++; sectionkey2SAP[key][0]=s; sectionkey2SAP[key][1]=curr->inorder2; }
        }
    }
    actual_max_states = key;
    #else
    for(s=1; s<=max_sections; s++) if(sectionkey[s]!=0) actual_max_states = s;
    #endif
    for(s=1; s<=actual_max_states; s++) {
        unsigned long int n = 0;
        #if NUM_POLYS == 2
        unsigned long int s1 = sectionkey2SAP[s][0], s2_in = sectionkey2SAP[s][1];
        curr = first_hinge_span[s1];
        while(curr->nexthinge) { curr = curr->nexthinge; if(curr->inorder2 == s2_in) n++; }
        #else
        curr = first_hinge_span[s];
        while(curr->nexthinge) { curr = curr->nexthinge; n++; }
        #endif
        num_outsections[s] = n;
        if(n > 0) { tspans_outsection[s]=unsgnlong_vecalloc(1,n); tspans_edges[s]=unsgnlong_vecalloc(1,n); }
    }
    unsigned long int curIdx[max_keynum+1];
    for(i=1; i<=actual_max_states; i++) curIdx[i]=1;
    for(s=1; s<=max_sections; s++) {
        curr = first_hinge_span[s];
        while(curr->nexthinge) {
            tofree = curr; curr = curr->nexthinge; free(tofree);
            unsigned long int srcKey = s, destKey=0;
            #if NUM_POLYS == 2
            for(j=1; j<=actual_max_states; j++) if(sectionkey2SAP[j][0]==s && sectionkey2SAP[j][1]==curr->inorder2) { srcKey=j; break; }
            for(j=1; j<=actual_max_states; j++) if(sectionkey2SAP[j][0]==curr->outorder && sectionkey2SAP[j][1]==curr->outorder2) { destKey=j; break; }
            #else
            destKey = curr->outorder;
            #endif
            if(srcKey && destKey) {
                tspans_outsection[srcKey][curIdx[srcKey]] = destKey;
                tspans_edges[srcKey][curIdx[srcKey]] = curr->edgecount + curr->edgecount2;
                curIdx[srcKey]++;
            }
        }
    }
    for(i=1; i<=actual_max_states; i++) qcksrtII(num_outsections[i], tspans_outsection[i], tspans_edges[i]);
}

/* ==========================================================================
   SOLVER & STATS
   ========================================================================== */

double max_eval_LRvec(double fugacity) {
    int i, k; unsigned long int in_s, out_s, nth; double L0, L1, R0, R1, pwf;
    for(i=1; i<=actual_max_states; i++) { L_Evector[0][i]=1; R_Evector[0][i]=1; }
    for(k=1; k<=500; k++) {
        L0=0.1; R0=0.1;
        for(i=1; i<=actual_max_states; i++) { L_Evector[1][i]=0; R_Evector[1][i]=0; }
        for(in_s=1; in_s<=actual_max_states; in_s++) {
            for(nth=1; nth<=num_outsections[in_s]; nth++) {
                pwf = pow(fugacity, tspans_edges[in_s][nth]);
                out_s = tspans_outsection[in_s][nth];
                L_Evector[1][out_s] += L_Evector[0][in_s]*pwf;
                R_Evector[1][in_s] += R_Evector[0][out_s]*pwf;
            }
        }
        for(i=1; i<=actual_max_states; i++) { if(L_Evector[1][i]>L0) L0=L_Evector[1][i]; if(R_Evector[1][i]>R0) R0=R_Evector[1][i]; }
        for(i=1; i<=actual_max_states; i++) { L_Evector[1][i]/=L0; R_Evector[1][i]/=R0; }
        L1=0.1; R1=0.1;
        for(i=1; i<=actual_max_states; i++) { L_Evector[0][i]=0; R_Evector[0][i]=0; }
        for(in_s=1; in_s<=actual_max_states; in_s++) {
            for(nth=1; nth<=num_outsections[in_s]; nth++) {
                pwf = pow(fugacity, tspans_edges[in_s][nth]);
                out_s = tspans_outsection[in_s][nth];
                L_Evector[0][out_s] += L_Evector[1][in_s]*pwf;
                R_Evector[0][in_s] += R_Evector[1][out_s]*pwf;
            }
        }
        for(i=1; i<=actual_max_states; i++) { if(L_Evector[0][i]>L1) L1=L_Evector[0][i]; if(R_Evector[0][i]>R1) R1=R_Evector[0][i]; }
        for(i=1; i<=actual_max_states; i++) { L_Evector[0][i]/=L1; R_Evector[0][i]/=R1; }
        if(fabs(L1-L0)<1e-8 && fabs(R1-R0)<1e-8) return R1 - 1.0;
    }
    exit(1);
}

double get_Beta(double kappa) {
    double B=0, pwf; unsigned long int i, j, out;
    for(i=1; i<=actual_max_states; i++) for(j=1; j<=num_outsections[i]; j++) {
        pwf = pow(kappa, tspans_edges[i][j]); out = tspans_outsection[i][j];
        B += L_Evector[0][i] * tspans_edges[i][j] * pwf * R_Evector[0][out];
    }
    return B;
}

void export_eigenvectors(void) {
    char filename[100], filename2[100];
    sprintf(filename, "L_Evector_TS_L%dM%d.txt", L, M);
    sprintf(filename2, "R_Evector_TS_L%dM%d.txt", L, M);
    FILE *fp = fopen(filename, "w"), *fp2 = fopen(filename2, "w");
    if(!fp || !fp2) { printf("Error opening eigenvector files for writing!\n"); return; }
    for(unsigned long int i=1; i<=actual_max_states; i++) {
        fprintf(fp, "%.15f\n", L_Evector[0][i]);
        fprintf(fp2, "%.15f\n", R_Evector[0][i]);
    }
    fclose(fp); fclose(fp2);
    printf("Eigenvectors exported to %s and %s\n", filename, filename2);
}

int main(void) {
    int ordNum[NUM_POLYS*3], i, j;
    printf("Master TMcalc Engine (v9) - Mode: %d, Lattice: %dx%d\n", SYSTEM_MODE, LAT_L, LAT_M);
    for (int p=0; p<NUM_POLYS; p++) { ordNum[p*3+0]=1; ordNum[p*3+1]=1; ordNum[p*3+2]=0; }
    for (i=1; i<=max_sections; i++) { current_hinge_span[i]=newhinge(); first_hinge_span[i]=current_hinge_span[i]; sectionkey[i]=0; }
    for (i=0; i<=LAT_M; i++) for (j=0; j<=LAT_L; j++) { if(!(i==LAT_M && j==LAT_L)) enterhinge(i, j, 0, 0, &ordNum); alreadyentered[0][i][j]=1; }
    conv_to_array();
    L_Evector[0]=vecalloc(0, actual_max_states); L_Evector[1]=vecalloc(0, actual_max_states);
    R_Evector[0]=vecalloc(0, actual_max_states); R_Evector[1]=vecalloc(0, actual_max_states);
    connectivity_inv = rtflsp(&max_eval_LRvec, 0.3, 0.7, 1e-8);
    double Beta = get_Beta(connectivity_inv), LmultR = 0;
    for(i=1; i<=actual_max_states; i++) LmultR += L_Evector[0][i]*R_Evector[0][i];
    double Alpha = pow(connectivity_inv, 4) * (L_Evector[0][1]*R_Evector[0][1]) / (LmultR * Beta);
    printf("\n--- RESEARCH RESULTS ---\nConnective Constant: %f\nAlpha (Entanglement): %e\nBeta (Normalization): %e\n", 1.0/connectivity_inv, Alpha, Beta);
    export_eigenvectors();
    return 0;
}
