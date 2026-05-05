/*
 * Unified transfer-matrix driver for self-avoiding polygon models on an
 * L-by-M tube. This file owns the combinatorial state generation: it walks all
 * possible two-span hinge configurations, assigns each boundary connectivity
 * pattern a compact section id, converts the resulting graph to CSR form, then
 * hands the numerical eigenvector solve to tm_spectral.c.
 *
 * Modes:
 *   0: one SAP, ordinary
 *   1: one SAP, Hamiltonian
 *   2: two SAPs on the same lattice, ordinary
 *   3: two SAPs on the same lattice, Hamiltonian
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <stdarg.h>
#include <omp.h>
#include <getopt.h>
#include <sys/stat.h>
#include "tm_spectral.h"

/* Forward Declarations */
struct hinge_span;
void enterhinge(int i, int j, int side, int poly, int (*pON)[]);
void leavehinge(int i, int j, int side, int poly, int (*pON)[]);
void rowedges(int i, int j, int poly, int (*pON)[]);
void coledges(int i, int j, int poly, int (*pON)[]);
void start_poly2(int (*pON)[]);
unsigned short int LFlag(int (*pON)[]);
unsigned short int LFlag2(int (*pON)[]);
void export_matrix(void);

/* ==========================================================================
   CONFIGURATION VARIABLES
   ========================================================================== */

int mode = 0;
int lat_L = 1;
int lat_M = 1;

#define L lat_L
#define M lat_M
#define vM (lat_M+1)
#define vL (lat_L+1)
#define force 0
#define SS

int num_polys = 1;
int ham_check = 0;
int damping_enabled = 0;
double convergence_threshold = 1e-8;
double manual_x = -1.0;
int sampling_export = 0;

/* Capacity limits are defaults for the generated graph. Large lattices can
   override them with -S and -K rather than recompiling archival constants. */
unsigned long int actual_max_sections = 0;
unsigned long int actual_max_keynum = 0;
unsigned long int user_max_sections = 0;
unsigned long int user_max_keynum = 0;

unsigned long int *sectionkey = NULL;
unsigned long int (*sectionkey2SAP)[2] = NULL;
unsigned long int *num_outsections = NULL;
unsigned long int *csr_out_states = NULL;
unsigned long int *csr_edges = NULL;
unsigned long int *csr_row_ptr = NULL;
TmSpectralProblem spectral_problem;
unsigned int *section_template_scratch = NULL;
size_t section_template_len = 0;

#define MAX_HEDGE ( (lat_M+1)*lat_L + lat_M*(lat_L+1) + 1 )

typedef double vec_ent;
typedef double mat_ent;

/* ==========================================================================
   STRUCTURES
   ========================================================================== */

struct hinge_span {
    /* For one SAP, a transition is "input section -> output section".
       For 2SAP, a state is the pair (poly A section, poly B section), so
       inorder2/outorder2 carry polygon B's section ids. */
    unsigned long int inorder2;
    unsigned long int outorder;
    unsigned long int outorder2;
    /* edgecount values are the powers of fugacity assigned to this transition. */
    unsigned int edgecount;
    unsigned int edgecount2;
    /* Fast duplicate prefilter; hedges[] remains the exact identity check. */
    unsigned long int hedge_hash;
    unsigned int *hedges;
    struct hinge_span *nexthinge;
};

struct hinge_span **first_hinge_span = NULL;
struct hinge_span **current_hinge_span = NULL;

/* ==========================================================================
   GLOBAL VARIABLES (Dynamic)

   ordertemplate stores the raw discovery numbering for boundary hinges.
   reordertemplate stores the canonical numbering used as the section key.
   Both are indexed [polygon][side][row][column], where side 0 is the left
   boundary of a two-span and side 1 is the right boundary.
   ========================================================================== */

unsigned int ****ordertemplate_master = NULL;
#define ordertemplate ordertemplate_master[0]
#define ordertemplate2 ordertemplate_master[1]

unsigned int ****reordertemplate_master = NULL;
#define reordertemplate reordertemplate_master[0]
#define reordertemplate2 reordertemplate_master[1]

unsigned long int num_2_spans = 0;
unsigned long int valid_2_spans = 0;
unsigned short int **hingestatus = NULL;
unsigned short int ***alreadyentered = NULL;
unsigned short int **colhingeedges = NULL;
unsigned short int **rowhingeedges = NULL;
unsigned long int dupcounter = 0;

vec_ent *L_Evector[2], *R_Evector[2];
double connectivity_inv;
unsigned long int actual_max_states = 0;

#include "tm_runtime.h"

/* ==========================================================================
   UTILITY INCLUDES
   ========================================================================== */

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "../deps/utils/unsgnlong_vecalloc.c"
#include "../deps/utils/vecalloc.c"
#include "../deps/utils/matalloc.c"

#include "../deps/topology/LFlag_norder.c"
#include "../deps/topology/LFlag_norder2.c"
#pragma GCC diagnostic pop

static unsigned long int hash_template_entry(unsigned long int hash, unsigned int value)
{
    hash ^= value;
    hash *= 1099511628211ULL;
    return hash;
}

/* FNV-1a hash of the canonical boundary template. Hashes are only a lookup
   accelerator: tm_runtime.h also stores the full template and compares it on
   collision, so section identity is exact rather than probabilistic. */
unsigned long int get_section_hash(int poly, int side) {
    unsigned long int h = 14695981039346656037ULL;
    for(int i=0; i<=lat_M; i++) {
        for(int j=0; j<=lat_L; j++) {
            h = hash_template_entry(h, reordertemplate_master[poly][side][i][j]);
        }
    }
    return h;
}

static unsigned long int copy_and_hash_section_template(int poly, int side, unsigned int *dest) {
    unsigned long int h = 14695981039346656037ULL;
    size_t idx = 0;
    for (int i = 0; i <= lat_M; i++) {
        for (int j = 0; j <= lat_L; j++) {
            unsigned int value = reordertemplate_master[poly][side][i][j];
            dest[idx++] = value;
            h = hash_template_entry(h, value);
        }
    }
    return h;
}

/* Get the stable section id for the current canonical template. The scratch
   buffer avoids repeated allocations during the recursive state search. */
static unsigned long int get_or_add_current_section(int poly, int side) {
    unsigned long int section_num;
    unsigned long int section_hash;

    section_hash = copy_and_hash_section_template(poly, side, section_template_scratch);
    section_num = get_or_add_section(section_hash, section_template_scratch, section_template_len);
    return section_num;
}

#define QM lat_M
#undef M
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "../deps/utils/qcksrtII_int.c"
#pragma GCC diagnostic pop
#define M QM
#include "../deps/numerics/rtflsp.c"

/* ==========================================================================
   CORE ENGINE
   ========================================================================== */

struct hinge_span* newhinge(void) {
    return (struct hinge_span *)xcalloc(1, sizeof(struct hinge_span), "hinge span");
}

void fillreordertemplate(int poly, int ledges, int redges) {
    int i, j, ii, jj, s;
    unsigned int firstentry, secondentry, connectingedge;
    for(s=0; s<=1; s++) {
        for(i=0; i<=lat_M; i++) for(j=0; j<=lat_L; j++) reordertemplate_master[poly][s][i][j]=0;
        firstentry=1;
        for(i=0; i<=lat_M; i++) {
            for(j=0; j<=lat_L; j++) {
                if(ordertemplate_master[poly][s][i][j] % 2 != 0 && reordertemplate_master[poly][s][i][j]==0) {
                    reordertemplate_master[poly][s][i][j] = firstentry;
                    connectingedge = (ordertemplate_master[poly][s][i][j]==1) ? (unsigned int)(s==0 ? ledges : redges) : (ordertemplate_master[poly][s][i][j]-1);
                    secondentry = (firstentry==1) ? (unsigned int)(s==0 ? ledges : redges) : (firstentry+1);
                    for(ii=0; ii<=lat_M; ii++) for(jj=0; jj<=lat_L; jj++) if(ordertemplate_master[poly][s][ii][jj]==connectingedge) {
                        reordertemplate_master[poly][s][ii][jj] = secondentry;
                        firstentry = (firstentry==1) ? 2 : (firstentry+2);
                        break;
                    }
                } else if(ordertemplate_master[poly][s][i][j]>0 && reordertemplate_master[poly][s][i][j]==0) {
                    reordertemplate_master[poly][s][i][j] = firstentry;
                    connectingedge = (s==0 && ordertemplate_master[poly][s][i][j]==(unsigned int)ledges) || (s==1 && ordertemplate_master[poly][s][i][j]==(unsigned int)redges) ? 1 : (ordertemplate_master[poly][s][i][j]+1);
                    secondentry = (firstentry==1) ? (unsigned int)(s==0 ? ledges : redges) : (firstentry+1);
                    for(ii=0; ii<=lat_M; ii++) for(jj=0; jj<=lat_L; jj++) if(ordertemplate_master[poly][s][ii][jj]==connectingedge) {
                        reordertemplate_master[poly][s][ii][jj] = secondentry;
                        firstentry = (firstentry==1) ? 2 : (firstentry+2);
                        break;
                    }
                }
            }
        }
    }
}

/* Store one valid two-span transition unless an identical transition was
   already found. The linked-list form is convenient while recursively
   discovering states; conv_to_array() later compacts it to CSR. */
void recordtemplate(int (*pON)[]) {
    unsigned long int inNum = get_or_add_current_section(0, 0);
    unsigned long int outNum = get_or_add_current_section(0, 1);
    int i, j;

    unsigned long int inNum2=0, outNum2=0;
    if (num_polys == 2) {
        inNum2 = get_or_add_current_section(1, 0);
        outNum2 = get_or_add_current_section(1, 1);
    }

    int *pON_raw = (int*)pON;
    unsigned int n_edges = pON_raw[0] + pON_raw[2] - 1;
    unsigned int n_edges2 = (num_polys==2) ? (pON_raw[3] + pON_raw[5] - 1) : 0;
    
    /* hedges[] is positional: every possible hinge edge consumes one slot,
       and unoccupied edges leave a zero. This preserves the archival identity
       test for duplicate transitions. */
    unsigned int temp_hedges[MAX_HEDGE];
    for(i=0; i<MAX_HEDGE; i++) temp_hedges[i]=0;
    int idx=0, edgenum=1;
    for(i=0; i<=lat_M; i++) for(j=0; j<lat_L; j++) {
        if (idx >= MAX_HEDGE) {
            fprintf(stderr, "Fatal: hinge edge list overflow while recording template for L=%d M=%d\n", lat_L, lat_M);
            exit(EXIT_FAILURE);
        }
        if(rowhingeedges[i][j]) temp_hedges[idx]=edgenum;
        idx++;
        edgenum++;
    }
    for(i=0; i<lat_M; i++) for(j=0; j<=lat_L; j++) {
        if (idx >= MAX_HEDGE) {
            fprintf(stderr, "Fatal: hinge edge list overflow while recording template for L=%d M=%d\n", lat_L, lat_M);
            exit(EXIT_FAILURE);
        }
        if(colhingeedges[i][j]) temp_hedges[idx]=edgenum;
        idx++;
        edgenum++;
    }
    unsigned long int hedge_hash = 14695981039346656037ULL;
    for(i=0; i<MAX_HEDGE; i++) hedge_hash = hash_template_entry(hedge_hash, temp_hedges[i]);

    int duplicate=0;
    struct hinge_span* test = first_hinge_span[inNum];
    while(test->nexthinge != NULL) {
        test = test->nexthinge;
        if(test->outorder==outNum && test->inorder2==inNum2 && test->outorder2==outNum2) {
            if(test->edgecount==n_edges && test->edgecount2==n_edges2 && test->hedge_hash==hedge_hash) {
                int match=1;
                for(i=0; i<MAX_HEDGE; i++) if(temp_hedges[i]!=test->hedges[i]) { match=0; break; }
                if(match) { duplicate=1; dupcounter++; break; }
            }
        }
    }
    if(!duplicate) {
        (*current_hinge_span[inNum]).nexthinge = (struct hinge_span*)xcalloc(1, sizeof(struct hinge_span), "hinge span link");
        current_hinge_span[inNum] = (*current_hinge_span[inNum]).nexthinge;
        (*current_hinge_span[inNum]).inorder2 = inNum2;
        (*current_hinge_span[inNum]).outorder = outNum;
        (*current_hinge_span[inNum]).outorder2 = outNum2;
        (*current_hinge_span[inNum]).edgecount = n_edges;
        (*current_hinge_span[inNum]).edgecount2 = n_edges2;
        (*current_hinge_span[inNum]).hedge_hash = hedge_hash;
        (*current_hinge_span[inNum]).hedges = (unsigned int*)xcalloc(MAX_HEDGE, sizeof(unsigned int), "hinge edge list");
        for(i=0; i<MAX_HEDGE; i++) (*current_hinge_span[inNum]).hedges[i]=temp_hedges[i];
    }
}

/* 2SAP discovery nests polygon B inside each partial polygon A discovery. Both
   polygons share hingestatus/edge occupancy, but use independent templates. */
void start_poly2(int (*pON)[]) {
    for (int i=0; i<=lat_M; i++) {
        for (int j=0; j<=lat_L; j++) {
            if(!hingestatus[i][j]) enterhinge(i, j, 0, 1, pON);
            alreadyentered[1][i][j]=1;
        }
    }
    for (int i=0; i<=lat_M; i++) for (int j=0; j<=lat_L; j++) alreadyentered[1][i][j]=0;
}

/* Enter/leave/rowedges/coledges are the recursive state generator. pON holds
   three counters per polygon: left boundary order, right boundary order, and
   hinge-edge count. */
void enterhinge(int i, int j, int side, int poly, int (*pON)[]) {
    if(alreadyentered[poly][i][j]==1 && side==0) return;
    int *pON_raw = (int*)pON;
    hingestatus[i][j]=1; ordertemplate_master[poly][side][i][j] = pON_raw[poly*3 + side]; pON_raw[poly*3 + side]++;
    leavehinge(i, j, side ^ 1, poly, pON); rowedges(i, j, poly, pON); coledges(i, j, poly, pON);
    pON_raw[poly*3 + side]--; ordertemplate_master[poly][side][i][j] = 0; hingestatus[i][j]=0;
}

void leavehinge(int i, int j, int side, int poly, int (*pON)[]) {
    if(alreadyentered[poly][i][j]==1 && side==0) return;
    int *pON_raw = (int*)pON;
    ordertemplate_master[poly][side][i][j] = pON_raw[poly*3 + side]; pON_raw[poly*3 + side]++;
    if (side == 0 && pON_raw[poly*3 + 1] > 1) {
        if (num_polys == 2 && poly == 0) {
            if (pON_raw[0]-1 <= vM * vL - 2) start_poly2(pON);
        } else {
            if (num_polys == 2 ? LFlag2(pON) : LFlag(pON)) {
                if (ham_check) {
                    int isHam=1;
                    for (int ii=0; ii<=lat_M; ii++) for (int jj=0; jj<=lat_L; jj++) if (!hingestatus[ii][jj]) { isHam=0; break; }
                    if(isHam) { fillreordertemplate(0, pON_raw[0]-1, pON_raw[1]-1); if(num_polys==2) fillreordertemplate(1, pON_raw[3]-1, pON_raw[4]-1); recordtemplate(pON); valid_2_spans++; }
                } else {
                    fillreordertemplate(0, pON_raw[0]-1, pON_raw[1]-1); if(num_polys==2) fillreordertemplate(1, pON_raw[3]-1, pON_raw[4]-1); recordtemplate(pON); valid_2_spans++;
                }
            }
        }
        num_2_spans++;
    }
    for (int ii=0; ii<=lat_M; ii++) for (int jj=0; jj<=lat_L; jj++) if (!hingestatus[ii][jj]) enterhinge(ii, jj, side, poly, pON);
    pON_raw[poly*3 + side]--; ordertemplate_master[poly][side][i][j] = 0;
}

void rowedges(int i, int j, int poly, int (*pON)[]) {
    int *pON_raw = (int*)pON;
    if(j>0 && !hingestatus[i][j-1]) { hingestatus[i][j-1]=1; pON_raw[poly*3 + 2]++; rowhingeedges[i][j-1]=1; leavehinge(i, j-1, 0, poly, pON); leavehinge(i, j-1, 1, poly, pON); rowedges(i, j-1, poly, pON); coledges(i, j-1, poly, pON); pON_raw[poly*3 + 2]--; hingestatus[i][j-1]=0; rowhingeedges[i][j-1]=0; }
    if(j<lat_L && !hingestatus[i][j+1]) { hingestatus[i][j+1]=1; pON_raw[poly*3 + 2]++; rowhingeedges[i][j]=1; leavehinge(i, j+1, 0, poly, pON); leavehinge(i, j+1, 1, poly, pON); rowedges(i, j+1, poly, pON); coledges(i, j+1, poly, pON); pON_raw[poly*3 + 2]--; hingestatus[i][j+1]=0; rowhingeedges[i][j]=0; }
}

void coledges(int i, int j, int poly, int (*pON)[]) {
    int *pON_raw = (int*)pON;
    if(i>0 && !hingestatus[i-1][j]) { hingestatus[i-1][j]=1; pON_raw[poly*3 + 2]++; colhingeedges[i-1][j]=1; leavehinge(i-1, j, 0, poly, pON); leavehinge(i-1, j, 1, poly, pON); rowedges(i-1, j, poly, pON); coledges(i-1, j, poly, pON); pON_raw[poly*3 + 2]--; hingestatus[i-1][j]=0; colhingeedges[i-1][j]=0; }
    if(i<lat_M && !hingestatus[i+1][j]) { hingestatus[i+1][j]=1; pON_raw[poly*3 + 2]++; colhingeedges[i][j]=1; leavehinge(i+1, j, 0, poly, pON); leavehinge(i+1, j, 1, poly, pON); rowedges(i+1, j, poly, pON); coledges(i+1, j, poly, pON); pON_raw[poly*3 + 2]--; hingestatus[i+1][j]=0; colhingeedges[i][j]=0; }
}

void allocate_resources(void) {
    /* These defaults are deliberately generous for audited small lattices.
       Researchers can raise them at runtime when exploring larger instances. */
    if (user_max_sections > 0) actual_max_sections = user_max_sections;
    else if (lat_L * lat_M == 1) actual_max_sections = 1000;
    else if (lat_L * lat_M == 2) actual_max_sections = 10000;
    else if (lat_L * lat_M == 3) actual_max_sections = 50000;
    else actual_max_sections = 500000;

    if (user_max_keynum > 0) actual_max_keynum = user_max_keynum;
    else if (lat_L * lat_M == 1) actual_max_keynum = 1000;
    else if (lat_L * lat_M == 2) actual_max_keynum = 50000;
    else if (lat_L * lat_M == 3) actual_max_keynum = 200000;
    else actual_max_keynum = 2000000;

    sectionkey = (unsigned long int*)xcalloc(actual_max_sections + 1, sizeof(unsigned long int), "section keys");
    sectionkey2SAP = (unsigned long int(*)[2])xcalloc(actual_max_keynum + 1, sizeof(unsigned long int[2]), "2SAP section keys");
    num_outsections = (unsigned long int*)xcalloc(actual_max_keynum + 1, sizeof(unsigned long int), "outsection counts");
    csr_row_ptr = (unsigned long int*)xcalloc(actual_max_keynum + 2, sizeof(unsigned long int), "CSR row pointer");
    first_hinge_span = (struct hinge_span**)xcalloc(actual_max_sections + 1, sizeof(struct hinge_span*), "first hinge span table");
    current_hinge_span = (struct hinge_span**)xcalloc(actual_max_sections + 1, sizeof(struct hinge_span*), "current hinge span table");
    section_template_len = (size_t)vM * (size_t)vL;
    section_template_scratch = (unsigned int *)xcalloc(section_template_len, sizeof(*section_template_scratch), "section template scratch");

    int p, s, i;
    ordertemplate_master = (unsigned int ****)xcalloc(num_polys, sizeof(unsigned int ***), "order templates");
    reordertemplate_master = (unsigned int ****)xcalloc(num_polys, sizeof(unsigned int ***), "reorder templates");
    for(p=0; p<num_polys; p++) {
        ordertemplate_master[p] = (unsigned int ***)xcalloc(2, sizeof(unsigned int **), "order template sides");
        reordertemplate_master[p] = (unsigned int ***)xcalloc(2, sizeof(unsigned int **), "reorder template sides");
        for(s=0; s<2; s++) {
            ordertemplate_master[p][s] = (unsigned int **)xcalloc(vM, sizeof(unsigned int *), "order template rows");
            reordertemplate_master[p][s] = (unsigned int **)xcalloc(vM, sizeof(unsigned int *), "reorder template rows");
            for(i=0; i<vM; i++) {
                ordertemplate_master[p][s][i] = (unsigned int *)xcalloc(vL, sizeof(unsigned int), "order template row");
                reordertemplate_master[p][s][i] = (unsigned int *)xcalloc(vL, sizeof(unsigned int), "reorder template row");
            }
        }
    }
    hingestatus = (unsigned short int **)xcalloc(vM, sizeof(unsigned short int *), "hinge status rows");
    colhingeedges = (unsigned short int **)xcalloc(vM, sizeof(unsigned short int *), "column hinge edge rows");
    rowhingeedges = (unsigned short int **)xcalloc(vM, sizeof(unsigned short int *), "row hinge edge rows");
    for(i=0; i<vM; i++) {
        hingestatus[i] = (unsigned short int *)xcalloc(vL, sizeof(unsigned short int), "hinge status row");
        colhingeedges[i] = (unsigned short int *)xcalloc(vL, sizeof(unsigned short int), "column hinge edge row");
        rowhingeedges[i] = (unsigned short int *)xcalloc(vL, sizeof(unsigned short int), "row hinge edge row");
    }
    alreadyentered = (unsigned short int ***)xcalloc(num_polys, sizeof(unsigned short int **), "already-entered tables");
    for(p=0; p<num_polys; p++) {
        alreadyentered[p] = (unsigned short int **)xcalloc(vM, sizeof(unsigned short int *), "already-entered rows");
        for(i=0; i<vM; i++) alreadyentered[p][i] = (unsigned short int *)xcalloc(vL, sizeof(unsigned short int), "already-entered row");
    }
}

#define HASH_SIZE 131072
struct hash_node {
    unsigned long int s1, s2;
    unsigned long int key;
    struct hash_node *next;
} *hash_table[HASH_SIZE];

/* In 2SAP mode the spectral matrix state is a pair of one-polygon sections.
   This map assigns a dense state id to each pair actually encountered. */
unsigned long int get_composite_key(unsigned long int s1, unsigned long int s2, unsigned long int *pKey) {
    unsigned long int h = (s1 * 31337 + s2) % HASH_SIZE;
    struct hash_node *n = hash_table[h];
    while(n) { if(n->s1 == s1 && n->s2 == s2) return n->key; n = n->next; }
    (*pKey)++;
    if (*pKey > actual_max_keynum) { printf("Fatal: actual_max_keynum overflow (%lu)\n", actual_max_keynum); exit(1); }
    n = (struct hash_node*)xmalloc(sizeof(struct hash_node), "composite hash node");
    n->s1 = s1; n->s2 = s2; n->key = *pKey; n->next = hash_table[h];
    hash_table[h] = n;
    sectionkey2SAP[*pKey][0] = s1; sectionkey2SAP[*pKey][1] = s2;
    return *pKey;
}

/* Convert the generated transition lists to compressed sparse row storage.
   Arrays are 1-based to match the archival Numerical Recipes-style helpers:
   row s occupies positions csr_row_ptr[s] + 1 through csr_row_ptr[s + 1]. */
void conv_to_array(void) {
    unsigned long int i, s, key=0;
    struct hinge_span *curr;
    for(i=0; i<HASH_SIZE; i++) hash_table[i]=NULL;

    if (num_polys == 2) {
        for(s=1; s<=actual_max_sections; s++) {
            if (!first_hinge_span[s]) continue;
            curr = first_hinge_span[s];
            while(curr->nexthinge) {
                curr = curr->nexthinge;
                get_composite_key(s, curr->inorder2, &key);
                get_composite_key(curr->outorder, curr->outorder2, &key);
            }
        }
        actual_max_states = key;
    } else {
        for(s=1; s<=actual_max_sections; s++) if(sectionkey[s]!=0) actual_max_states = s;
    }

    unsigned long int total_transitions = 0;
    for(s=1; s<=actual_max_states; s++) {
        unsigned long int n = 0;
        if (num_polys == 2) {
            unsigned long int s1 = sectionkey2SAP[s][0], s2_in = sectionkey2SAP[s][1];
            curr = first_hinge_span[s1];
            while(curr->nexthinge) { curr = curr->nexthinge; if(curr->inorder2 == s2_in) n++; }
        } else {
            curr = first_hinge_span[s];
            while(curr->nexthinge) { curr = curr->nexthinge; n++; }
        }
        num_outsections[s] = n;
        csr_row_ptr[s+1] = csr_row_ptr[s] + n;
        total_transitions += n;
    }
    csr_out_states = (unsigned long int*)xcalloc(total_transitions + 1, sizeof(unsigned long int), "CSR out states");
    csr_edges = (unsigned long int*)xcalloc(total_transitions + 1, sizeof(unsigned long int), "CSR edge weights");

    unsigned long int *curIdx = (unsigned long int*)xcalloc(actual_max_keynum + 1, sizeof(unsigned long int), "CSR write cursors");
    for(i=1; i<=actual_max_states; i++) curIdx[i]=1;
    for(s=1; s<=actual_max_sections; s++) {
        if (!first_hinge_span[s]) continue;
        curr = first_hinge_span[s]->nexthinge;
        while(curr) {
            struct hinge_span *next = curr->nexthinge;
            unsigned long int srcKey=0, destKey=0;
            if (num_polys == 2) {
                unsigned long int h = (s * 31337 + curr->inorder2) % HASH_SIZE;
                struct hash_node *n = hash_table[h];
                while(n) { if(n->s1 == s && n->s2 == curr->inorder2) { srcKey = n->key; break; } n = n->next; }
                h = (curr->outorder * 31337 + curr->outorder2) % HASH_SIZE;
                n = hash_table[h];
                while(n) { if(n->s1 == curr->outorder && n->s2 == curr->outorder2) { destKey = n->key; break; } n = n->next; }
            } else {
                srcKey = s;
                destKey = curr->outorder;
            }
            if(!srcKey || !destKey) {
                fprintf(stderr, "Fatal: missing composite CSR key while converting transition (%lu,%lu)->(%lu,%lu)\n",
                    s, curr->inorder2, curr->outorder, curr->outorder2);
                exit(EXIT_FAILURE);
            }
            {
                unsigned long int offset = csr_row_ptr[srcKey] + (curIdx[srcKey] - 1);
                if (offset + 1 > csr_row_ptr[srcKey + 1]) {
                    fprintf(stderr, "Fatal: CSR row cursor overflow while writing state %lu\n", srcKey);
                    exit(EXIT_FAILURE);
                }
                csr_out_states[offset + 1] = destKey;
                csr_edges[offset + 1] = curr->edgecount + curr->edgecount2;
                curIdx[srcKey]++;
            }
            free(curr->hedges);
            free(curr);
            curr = next;
        }
        free(first_hinge_span[s]);
        first_hinge_span[s] = NULL;
        current_hinge_span[s] = NULL;
    }
    free(curIdx);
    for(i=1; i<=actual_max_states; i++) if(num_outsections[i] > 0) qcksrtII(num_outsections[i], &csr_out_states[csr_row_ptr[i]], &csr_edges[csr_row_ptr[i]]);
}

double max_eval_LRvec(double fugacity) {
    return tm_spectral_max_eval(&spectral_problem, fugacity);
}

double get_Beta(double kappa) {
    return tm_spectral_beta(&spectral_problem, kappa);
}

static void find_root_bracket(double (*func)(double), double *x_low, double *x_high)
{
    /* The legacy solver expects a bracketed root. Start with the archival
       interval and expand only when a new lattice/mode falls outside it. */
    const int max_expansions = 24;
    double low = 0.3;
    double high = 0.8;
    double f_low = func(low);
    double f_high = func(high);

    for (int step = 0; step < max_expansions && f_low * f_high > 0.0; step++) {
        if (fabs(f_low) < fabs(f_high)) {
            low *= 0.5;
            f_low = func(low);
        } else {
            high *= 1.5;
            f_high = func(high);
        }
        if (!isfinite(f_low) || !isfinite(f_high)) {
            break;
        }
    }

    if (!isfinite(f_low) || !isfinite(f_high) || f_low * f_high > 0.0) {
        fprintf(stderr, "Fatal: unable to bracket spectral root; f(%g)=%g, f(%g)=%g\n",
            low, f_low, high, f_high);
        exit(EXIT_FAILURE);
    }

    *x_low = low;
    *x_high = high;
}

void export_eigenvectors(void) {
    const char *mode_suffixes[] = {"std", "ham", "2sap", "2sap_ham"};
    const char *suffix = (mode >= 0 && mode <= 3) ? mode_suffixes[mode] : "unk";

    ensure_directory("data");
    ensure_directory("data/TransferMatrix");
    ensure_directory("data/TransferMatrix/TMresults");
    
    if (sampling_export) {
        char fn[128];
        const char *mc_prefix = (mode == 2 || mode == 3) ? "2SAP_" : "";
        const char *mc_ham = (mode == 1 || mode == 3) ? "Ham" : "";
        checked_snprintf(fn, sizeof(fn), "data/TransferMatrix/TMresults/%sR_Evector%s_TS_L%dM%d.txt", mc_prefix, mc_ham, lat_L, lat_M);
        
        FILE *fp = xfopen(fn, "w");
        
        unsigned long int i, j;
        for (i = 1; i <= actual_max_states; i++) {
            unsigned long int row_start = csr_row_ptr[i];
            unsigned long int row_end = csr_row_ptr[i+1];
            for (j = row_start + 1; j <= row_end; j++) {
                /* The Monte Carlo sampler is transition-indexed: entry k is
                   the right eigenvector value of the state reached by
                   transition k, not of the source state. */
                fprintf(fp, "%.15f\n", R_Evector[0][csr_out_states[j]]);
            }
        }
        fclose(fp);
        printf("[Phase 3] Transition-indexed eigenvector exported to %s\n", fn);
        return;
    }

    char fn1[100], fn2[100];
    checked_snprintf(fn1, sizeof(fn1), "data/TransferMatrix/TMresults/L_Evector_L%dM%d_%s.txt", lat_L, lat_M, suffix);
    checked_snprintf(fn2, sizeof(fn2), "data/TransferMatrix/TMresults/R_Evector_L%dM%d_%s.txt", lat_L, lat_M, suffix);
    FILE *fp1 = xfopen(fn1, "w"), *fp2 = xfopen(fn2, "w");
    for(unsigned long int i=1; i<=actual_max_states; i++) { fprintf(fp1, "%.15f\n", L_Evector[0][i]); fprintf(fp2, "%.15f\n", R_Evector[0][i]); }
    fclose(fp1); fclose(fp2);
}

void export_matrix(void) {
    const char *mode_suffixes[] = {"std", "ham", "2sap", "2sap_ham"};
    const char *suffix = (mode >= 0 && mode <= 3) ? mode_suffixes[mode] : "unk";
    char fn[100];
    ensure_directory("data");
    ensure_directory("data/TransferMatrix");
    ensure_directory("data/TransferMatrix/TMresults");
    checked_snprintf(fn, sizeof(fn), "data/TransferMatrix/TMresults/CSR_L%dM%d_%s.bin", lat_L, lat_M, suffix);
    FILE *fp = xfopen(fn, "wb");
    unsigned long int total_transitions = csr_row_ptr[actual_max_states + 1];
    checked_fwrite(&actual_max_states, sizeof(unsigned long int), 1, fp, "state count");
    checked_fwrite(&total_transitions, sizeof(unsigned long int), 1, fp, "transition count");
    checked_fwrite(csr_row_ptr, sizeof(unsigned long int), actual_max_states + 2, fp, "CSR row pointer");
    checked_fwrite(csr_out_states, sizeof(unsigned long int), total_transitions + 1, fp, "CSR out states");
    checked_fwrite(csr_edges, sizeof(unsigned long int), total_transitions + 1, fp, "CSR edge weights");
    fclose(fp);
}

void cleanup_resources(void) {
    int p, s, i;
    for(p=0; p<num_polys; p++) {
        for(s=0; s<2; s++) {
            for(i=0; i<vM; i++) {
                free(ordertemplate_master[p][s][i]);
                free(reordertemplate_master[p][s][i]);
            }
            free(ordertemplate_master[p][s]);
            free(reordertemplate_master[p][s]);
        }
        free(ordertemplate_master[p]);
        free(reordertemplate_master[p]);
        for(i=0; i<vM; i++) free(alreadyentered[p][i]);
        free(alreadyentered[p]);
    }
    free(ordertemplate_master); free(reordertemplate_master); free(alreadyentered);
    for(i=0; i<vM; i++) { free(hingestatus[i]); free(colhingeedges[i]); free(rowhingeedges[i]); }
    free(hingestatus); free(colhingeedges); free(rowhingeedges);
    free(sectionkey); free(sectionkey2SAP); free(num_outsections); free(csr_row_ptr);
    free(section_template_scratch);
    free(first_hinge_span); free(current_hinge_span);
    free(csr_out_states); free(csr_edges);
    free(L_Evector[0]); free(L_Evector[1]); free(R_Evector[0]); free(R_Evector[1]);
    for(i=0; i<HASH_SIZE; i++) {
        struct hash_node *curr = hash_table[i], *tmp;
        while(curr) { tmp = curr; curr = curr->next; free(tmp); }
    }
    for(i=0; i<SECTION_HASH_SIZE; i++) {
        struct section_hash_node *curr = section_hash_table[i], *tmp;
        while(curr) { tmp = curr; curr = curr->next; free(tmp->template_values); free(tmp); }
    }
}

void print_usage(char *prog) { 
    printf("Usage: %s [-L lat_L] [-M lat_M] [-m mode] [-S max_sec] [-K max_key] [-c conv] [-d] [-x x_val] [-E]\n", prog);
    printf("  -x x_val : Manually set fugacity (skips root finding, e.g., -x 1.0 for sampling)\n");
    printf("  -E       : Export transition-indexed eigenvectors for Monte Carlo sampling\n");
}

int main(int argc, char **argv) {
    int opt;
    while ((opt = getopt(argc, argv, "L:M:m:S:K:c:dx:Eh")) != -1) {
        char *endptr;
        switch (opt) {
            case 'L': 
                lat_L = (int)strtol(optarg, &endptr, 10);
                if (*endptr != '\0') { fprintf(stderr, "Error: Lattice Width (-L) must be an integer (received '%s')\n", optarg); exit(1); }
                break;
            case 'M': 
                lat_M = (int)strtol(optarg, &endptr, 10);
                if (*endptr != '\0') { fprintf(stderr, "Error: Lattice Height (-M) must be an integer (received '%s')\n", optarg); exit(1); }
                break;
            case 'm': 
                mode = (int)strtol(optarg, &endptr, 10);
                if (*endptr != '\0') { fprintf(stderr, "Error: Mode (-m) must be an integer (received '%s')\n", optarg); exit(1); }
                break;
            case 'S': 
                user_max_sections = strtoul(optarg, &endptr, 10);
                if (*endptr != '\0') { fprintf(stderr, "Error: Max Sections (-S) must be an integer (received '%s')\n", optarg); exit(1); }
                break;
            case 'K': 
                user_max_keynum = strtoul(optarg, &endptr, 10);
                if (*endptr != '\0') { fprintf(stderr, "Error: Max Keynum (-K) must be an integer (received '%s')\n", optarg); exit(1); }
                break;
            case 'c': 
                convergence_threshold = strtod(optarg, &endptr);
                if (*endptr != '\0') { fprintf(stderr, "Error: Convergence (-c) must be a number (received '%s')\n", optarg); exit(1); }
                break;
            case 'd': damping_enabled = 1; break;
            case 'x':
                manual_x = strtod(optarg, &endptr);
                if (*endptr != '\0') { fprintf(stderr, "Error: Fugacity (-x) must be a number (received '%s')\n", optarg); exit(1); }
                break;
            case 'E': sampling_export = 1; break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }
    if (mode == 2 || mode == 3) num_polys = 2;
    if (mode == 1 || mode == 3) { ham_check = 1; damping_enabled = 1; }

    /* Validate before allocation so bad command lines fail without creating
       partial output files. */
    int validation_failed = 0;
    if (lat_L < 0) { fprintf(stderr, "Error: Lattice Width (-L) must be >= 0 (received %d)\n", lat_L); validation_failed = 1; }
    if (lat_M < 0) { fprintf(stderr, "Error: Lattice Height (-M) must be >= 0 (received %d)\n", lat_M); validation_failed = 1; }
    if (lat_L == 0 && lat_M == 0) { fprintf(stderr, "Error: At least one dimension (L or M) must be > 0.\n"); validation_failed = 1; }
    if (mode < 0 || mode > 3) { fprintf(stderr, "Error: Invalid Simulation Mode (-m %d). Valid modes are 0, 1, 2, 3.\n", mode); validation_failed = 1; }
    if (convergence_threshold <= 0) { fprintf(stderr, "Error: Convergence threshold (-c) must be positive.\n"); validation_failed = 1; }
    
    if (validation_failed) {
        printf("\n");
        print_usage(argv[0]);
        exit(1);
    }

    double t_start = omp_get_wtime();
    int ordNum[6], i, j;
    unsigned long int state_idx;
    printf("[Phase 1] Initializing resources and generating state space...\n");
    allocate_resources();
    /* ordNum is laid out as [poly0 left, poly0 right, poly0 hinge_edges,
       poly1 left, poly1 right, poly1 hinge_edges]. */
    for (int p=0; p<num_polys; p++) { ordNum[p*3+0]=1; ordNum[p*3+1]=1; ordNum[p*3+2]=0; }
    for (state_idx=1; state_idx<=actual_max_sections; state_idx++) { first_hinge_span[state_idx]=newhinge(); current_hinge_span[state_idx]=first_hinge_span[state_idx]; sectionkey[state_idx]=0; }
    /* Try every lattice vertex except the upper-right sentinel as the first
       hinge. alreadyentered prevents generating the same state from an earlier
       starting vertex. */
    for (i=0; i<=lat_M; i++) for (j=0; j<=lat_L; j++) { if(!(i==lat_M && j==lat_L)) enterhinge(i, j, 0, 0, (int(*)[])&ordNum); alreadyentered[0][i][j]=1; }
    conv_to_array();
    double t_gen = omp_get_wtime();
    printf("[Phase 2] State space complete (%lu states). Solving eigenvalue problem...\n", actual_max_states);

    L_Evector[0]=vecalloc(0, actual_max_states); L_Evector[1]=vecalloc(0, actual_max_states);
    R_Evector[0]=vecalloc(0, actual_max_states); R_Evector[1]=vecalloc(0, actual_max_states);

    memset(&spectral_problem, 0, sizeof(spectral_problem));
    spectral_problem.states = actual_max_states;
    spectral_problem.row_ptr = csr_row_ptr;
    spectral_problem.out_states = csr_out_states;
    spectral_problem.edges = csr_edges;
    spectral_problem.left[0] = L_Evector[0];
    spectral_problem.left[1] = L_Evector[1];
    spectral_problem.right[0] = R_Evector[0];
    spectral_problem.right[1] = R_Evector[1];
    spectral_problem.convergence_threshold = convergence_threshold;
    spectral_problem.damping_enabled = damping_enabled;
    tm_spectral_prepare_transpose(&spectral_problem);
    
    if (manual_x > 0) {
        connectivity_inv = manual_x;
        double max_Evalue = max_eval_LRvec(connectivity_inv) + 1.0;
        printf("[Phase 3] Solved at manual fugacity x=%.6f. Dominant Eigenvalue=%.15f\n", connectivity_inv, max_Evalue);
    } else {
        double bracket_low, bracket_high;
        find_root_bracket(&max_eval_LRvec, &bracket_low, &bracket_high);
        connectivity_inv = rtflsp(&max_eval_LRvec, bracket_low, bracket_high, 1e-8);
    }
    printf("[Phase 3] Convergence reached. Aligning with thesis specification...\n");
    
    double LmultR = 0;
    for(state_idx=1; state_idx<=actual_max_states; state_idx++) LmultR += L_Evector[0][state_idx]*R_Evector[0][state_idx];

    /* Thesis normalization: scale Beta by the left/right inner product, then
       derive Alpha from the critical fugacity. */
    double Beta = get_Beta(connectivity_inv) / LmultR;
    double Alpha = connectivity_inv / Beta;
    double t_solve = omp_get_wtime();

    const char *mode_names[] = {"Standard", "Hamiltonian", "2SAP", "2SAP-Hamiltonian"};
    const char *mode_name = (mode >= 0 && mode <= 3) ? mode_names[mode] : "Unknown";

    printf("\n--- RESEARCH RESULTS (%dx%d %s) ---\n", lat_L, lat_M, mode_name);
    printf("Critical Fugacity (x_0):  %f\n", connectivity_inv);
    printf("Connective Constant (mu): %f\n", 1.0/connectivity_inv);
    printf("Alpha (Amplitude):        %e\n", Alpha);
    printf("Beta (Growth Parameter):  %e\n", Beta);

    unsigned long int total_trans = csr_row_ptr[actual_max_states+1];
    printf("\n--- PERFORMANCE TELEMETRY ---\n");
    printf("State Generation: %.3f s\n", t_gen - t_start);
    printf("Matrix Solving:   %.3f s (Total: %.3f s)\n", t_solve - t_gen, t_solve - t_start);
    printf("Transitions:      %lu (TPS: %.2e)\n", total_trans, (double)total_trans / (t_solve - t_gen));

    export_eigenvectors();
    export_matrix();
    tm_spectral_free(&spectral_problem);
    cleanup_resources();
    return 0;
}
