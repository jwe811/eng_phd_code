#include "mc_globals.h"



struct hinge_span **first_hinge_span, **current_hinge_span;
struct endhinge **firstendhinge;
struct endhinge **currentendhinge;

unsigned int ***ordertemplate;
unsigned long int num_2_spans = 1;
unsigned long int valid_2_spans = 0;

unsigned long int *num_outsections;
unsigned long int **t_outsection;
unsigned long int **t_nrr;
unsigned int **t_num_walks;
int ****t_start;
int ****t_end;
int ****t_walks;
unsigned int **tspans_edges;
double *L_Evector[2];
double *R_Evector_solve[2];

unsigned long int *num_left_endhinges;
unsigned int **Lend_num_walks;
int ****Lend_start;
int ****Lend_end;
int ****Lend_walks;

unsigned long int *num_right_endhinges;
unsigned int **Rend_num_walks;
int ****Rend_start;
int ****Rend_end;
int ****Rend_walks;

unsigned int ***reordertemplate;
unsigned short int **hingestatus;
unsigned short int **alreadyentered;
unsigned short int **colhingeedges;
unsigned short int **rowhingeedges;
unsigned long int *sectionkey;
unsigned long int dupcounter=0;

int **curstart;
int **curend;
int **curwalks;
unsigned int num_walks=0;

int EndOrdNum[2];
int **endtemplate;
int **rendtemplate;
int **endtemplate2;
int **rendtemplate2;

unsigned long int tot_right_endhinges=0;
unsigned long int tot_left_endhinges=0;
unsigned long int num_duplicate_right_endhinges=0;
unsigned long int num_duplicate_left_endhinges=0;
unsigned long int num_tspans=0;
unsigned long int num_duplicate_tspans=0;

int **built_walks_start;
int **built_walks_end;
int **built_walks_direcs;
int num_built_walks;

char filename[100];
FILE* fp;
unsigned long int filetotal=0;
unsigned long int filenum=1;



unsigned int*** allocate_3d_uint(int d1, int d2, int d3) {
    unsigned int ***arr = (unsigned int***)malloc(d1 * sizeof(unsigned int**));
    for (int i = 0; i < d1; i++) {
        arr[i] = (unsigned int**)malloc(d2 * sizeof(unsigned int*));
        for (int j = 0; j < d2; j++) {
            arr[i][j] = (unsigned int*)calloc(d3, sizeof(unsigned int));
        }
    }
    return arr;
}

unsigned short int** allocate_2d_ushort(int d1, int d2) {
    unsigned short int **arr = (unsigned short int**)malloc(d1 * sizeof(unsigned short int*));
    for (int i = 0; i < d1; i++) {
        arr[i] = (unsigned short int*)calloc(d2, sizeof(unsigned short int));
    }
    return arr;
}

int** allocate_2d_int(int d1, int d2) {
    int **arr = (int**)malloc(d1 * sizeof(int*));
    for (int i = 0; i < d1; i++) {
        arr[i] = (int*)calloc(d2, sizeof(int));
    }
    return arr;
}

void allocate_globals() {
    ordertemplate = allocate_3d_uint(2, vM, vL);
    reordertemplate = allocate_3d_uint(2, vM, vL);
    hingestatus = allocate_2d_ushort(vM, vL);
    alreadyentered = allocate_2d_ushort(vM, vL);
    colhingeedges = allocate_2d_ushort(vM, vL); // Used as colhingeedges[M][vL], vM is enough
    rowhingeedges = allocate_2d_ushort(vM, vL);
    sectionkey = (unsigned long int*)calloc(vec_length, sizeof(unsigned long int));
    
    num_outsections = (unsigned long int*)calloc(vec_length, sizeof(unsigned long int));
    t_outsection = (unsigned long int**)calloc(vec_length, sizeof(unsigned long int*));
    t_nrr = (unsigned long int**)calloc(vec_length, sizeof(unsigned long int*));
    t_num_walks = (unsigned int**)calloc(vec_length, sizeof(unsigned int*));
    t_start = (int****)calloc(vec_length, sizeof(int***));
    t_end = (int****)calloc(vec_length, sizeof(int***));
    t_walks = (int****)calloc(vec_length, sizeof(int***));
    tspans_edges = (unsigned int**)calloc(vec_length, sizeof(unsigned int*));
    
    L_Evector[0] = (double*)calloc(max_tspans + 1, sizeof(double));
    L_Evector[1] = (double*)calloc(max_tspans + 1, sizeof(double));
    R_Evector_solve[0] = (double*)calloc(max_tspans + 1, sizeof(double));
    R_Evector_solve[1] = (double*)calloc(max_tspans + 1, sizeof(double));
    
    num_left_endhinges = (unsigned long int*)calloc(vec_length, sizeof(unsigned long int));
    Lend_num_walks = (unsigned int**)calloc(vec_length, sizeof(unsigned int*));
    Lend_start = (int****)calloc(vec_length, sizeof(int***));
    Lend_end = (int****)calloc(vec_length, sizeof(int***));
    Lend_walks = (int****)calloc(vec_length, sizeof(int***));
    
    num_right_endhinges = (unsigned long int*)calloc(vec_length, sizeof(unsigned long int));
    Rend_num_walks = (unsigned int**)calloc(vec_length, sizeof(unsigned int*));
    Rend_start = (int****)calloc(vec_length, sizeof(int***));
    Rend_end = (int****)calloc(vec_length, sizeof(int***));
    Rend_walks = (int****)calloc(vec_length, sizeof(int***));
    
    first_hinge_span = (struct hinge_span**)calloc(vec_length, sizeof(struct hinge_span*));
    current_hinge_span = (struct hinge_span**)calloc(vec_length, sizeof(struct hinge_span*));
    firstendhinge = (struct endhinge**)calloc(vec_length, sizeof(struct endhinge*));
    currentendhinge = (struct endhinge**)calloc(vec_length, sizeof(struct endhinge*));
    
    curstart = allocate_2d_int(3, vM * vL);
    curend = allocate_2d_int(3, vM * vL);
    curwalks = allocate_2d_int(vM * vL + 1, vM * vL);
    
    endtemplate = allocate_2d_int(vM, vL);
    rendtemplate = allocate_2d_int(vM, vL);
    endtemplate2 = allocate_2d_int(vM, vL);
    rendtemplate2 = allocate_2d_int(vM, vL);
    
    built_walks_start = allocate_2d_int(vM * vL / 2 + 1, 3);
    built_walks_end = allocate_2d_int(vM * vL / 2 + 1, 3);
    built_walks_direcs = allocate_2d_int(vM * vL / 2 + 1, vM * vL * (totalspan + 1) + 1);
}
