#include "mc_globals.h"
#include "mc_runtime.h"



struct hinge_span **first_hinge_span, **current_hinge_span;
struct endhinge **firstendhinge;
struct endhinge **currentendhinge;
PolygonState primary_polygon_state;

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
    unsigned int ***arr = (unsigned int***)mc_xmalloc(d1 * sizeof(unsigned int**), "3D unsigned int outer array");
    unsigned int **planes = (unsigned int**)mc_xmalloc((size_t)d1 * d2 * sizeof(unsigned int*), "3D unsigned int plane table");
    unsigned int *data = (unsigned int*)mc_xcalloc((size_t)d1 * d2 * d3, sizeof(unsigned int), "3D unsigned int data");

    for (int i = 0; i < d1; i++) {
        arr[i] = planes + (size_t)i * d2;
        for (int j = 0; j < d2; j++) {
            arr[i][j] = data + ((size_t)i * d2 + j) * d3;
        }
    }
    return arr;
}

unsigned short int** allocate_2d_ushort(int d1, int d2) {
    unsigned short int **arr = (unsigned short int**)mc_xmalloc(d1 * sizeof(unsigned short int*), "2D unsigned short row table");
    unsigned short int *data = (unsigned short int*)mc_xcalloc((size_t)d1 * d2, sizeof(unsigned short int), "2D unsigned short data");

    for (int i = 0; i < d1; i++) {
        arr[i] = data + (size_t)i * d2;
    }
    return arr;
}

int** allocate_2d_int(int d1, int d2) {
    int **arr = (int**)mc_xmalloc(d1 * sizeof(int*), "2D int row table");
    int *data = (int*)mc_xcalloc((size_t)d1 * d2, sizeof(int), "2D int data");

    for (int i = 0; i < d1; i++) {
        arr[i] = data + (size_t)i * d2;
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
    sectionkey = (unsigned long int*)mc_xcalloc(vec_length, sizeof(unsigned long int), "section keys");
    
    num_outsections = (unsigned long int*)mc_xcalloc(vec_length, sizeof(unsigned long int), "outsection counts");
    t_outsection = (unsigned long int**)mc_xcalloc(vec_length, sizeof(unsigned long int*), "tspan outsections");
    t_nrr = (unsigned long int**)mc_xcalloc(vec_length, sizeof(unsigned long int*), "tspan numbers");
    t_num_walks = (unsigned int**)mc_xcalloc(vec_length, sizeof(unsigned int*), "tspan walk counts");
    t_start = (int****)mc_xcalloc(vec_length, sizeof(int***), "tspan starts");
    t_end = (int****)mc_xcalloc(vec_length, sizeof(int***), "tspan ends");
    t_walks = (int****)mc_xcalloc(vec_length, sizeof(int***), "tspan walks");
    tspans_edges = (unsigned int**)mc_xcalloc(vec_length, sizeof(unsigned int*), "tspan edges");
    
    L_Evector[0] = (double*)mc_xcalloc(max_tspans + 1, sizeof(double), "left eigenvector");
    L_Evector[1] = (double*)mc_xcalloc(max_tspans + 1, sizeof(double), "left eigenvector scratch");
    R_Evector_solve[0] = (double*)mc_xcalloc(max_tspans + 1, sizeof(double), "right eigenvector");
    R_Evector_solve[1] = (double*)mc_xcalloc(max_tspans + 1, sizeof(double), "right eigenvector scratch");
    
    num_left_endhinges = (unsigned long int*)mc_xcalloc(vec_length, sizeof(unsigned long int), "left endhinge counts");
    Lend_num_walks = (unsigned int**)mc_xcalloc(vec_length, sizeof(unsigned int*), "left endhinge walk counts");
    Lend_start = (int****)mc_xcalloc(vec_length, sizeof(int***), "left endhinge starts");
    Lend_end = (int****)mc_xcalloc(vec_length, sizeof(int***), "left endhinge ends");
    Lend_walks = (int****)mc_xcalloc(vec_length, sizeof(int***), "left endhinge walks");
    
    num_right_endhinges = (unsigned long int*)mc_xcalloc(vec_length, sizeof(unsigned long int), "right endhinge counts");
    Rend_num_walks = (unsigned int**)mc_xcalloc(vec_length, sizeof(unsigned int*), "right endhinge walk counts");
    Rend_start = (int****)mc_xcalloc(vec_length, sizeof(int***), "right endhinge starts");
    Rend_end = (int****)mc_xcalloc(vec_length, sizeof(int***), "right endhinge ends");
    Rend_walks = (int****)mc_xcalloc(vec_length, sizeof(int***), "right endhinge walks");
    
    first_hinge_span = (struct hinge_span**)mc_xcalloc(vec_length, sizeof(struct hinge_span*), "first hinge span table");
    current_hinge_span = (struct hinge_span**)mc_xcalloc(vec_length, sizeof(struct hinge_span*), "current hinge span table");
    firstendhinge = (struct endhinge**)mc_xcalloc(vec_length, sizeof(struct endhinge*), "first endhinge table");
    currentendhinge = (struct endhinge**)mc_xcalloc(vec_length, sizeof(struct endhinge*), "current endhinge table");
    
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

    primary_polygon_state.ordertemplate = ordertemplate;
    primary_polygon_state.reordertemplate = reordertemplate;
    primary_polygon_state.hingestatus = hingestatus;
    primary_polygon_state.alreadyentered = alreadyentered;
    primary_polygon_state.colhingeedges = colhingeedges;
    primary_polygon_state.rowhingeedges = rowhingeedges;
    primary_polygon_state.curstart = curstart;
    primary_polygon_state.curend = curend;
    primary_polygon_state.curwalks = curwalks;
    primary_polygon_state.num_walks = &num_walks;
}
