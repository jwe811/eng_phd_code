import os
import re

IN_FILE = "src/MASTER_MCsample.c"
OUT_DIR = "src"

with open(IN_FILE, "r") as f:
    lines = f.readlines()

def get_lines(start, end):
    return "".join(lines[start-1:end])

globals_h = """#ifndef MC_GLOBALS_H
#define MC_GLOBALS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include "../include/marsaglia.h"

// Macros
#define vec_ent double
#define mat_ent double
#define newline printf("\\n")
#define VALID 1
#define NOTVALID !VALID

// Global variable declarations
extern int L;
extern int M;
extern int totalspan;
extern int samplesize;
extern int runnum;
extern unsigned int seednum;
extern int max_sections;
extern int max_tspans;
extern double dom_evalue;
extern int maxpolys;
extern int vec_length;
extern int vM;
extern int vL;
extern int CS_mode;
extern int mode;
extern int ham_check;
extern char* output_dir;
extern double fval;

// Structs
""" + get_lines(64, 104) + """
extern struct endhinge **firstendhinge;
extern struct endhinge **currentendhinge;

// Arrays
extern unsigned int ***ordertemplate;
extern unsigned long int num_2_spans;
extern unsigned long int valid_2_spans;

extern unsigned long int *num_outsections;
extern unsigned long int **t_outsection;
extern unsigned long int **t_nrr;
extern unsigned int **t_num_walks;
extern int ****t_start;
extern int ****t_end;
extern int ****t_walks;
extern unsigned int **tspans_edges;
extern double *L_Evector[2];
extern double *R_Evector_solve[2];

extern unsigned long int *num_left_endhinges;
extern unsigned int **Lend_num_walks;
extern int ****Lend_start;
extern int ****Lend_end;
extern int ****Lend_walks;

extern unsigned long int *num_right_endhinges;
extern unsigned int **Rend_num_walks;
extern int ****Rend_start;
extern int ****Rend_end;
extern int ****Rend_walks;

extern unsigned int ***reordertemplate;
extern unsigned short int **hingestatus;
extern unsigned short int **alreadyentered;
extern unsigned short int **colhingeedges;
extern unsigned short int **rowhingeedges;
extern unsigned long int *sectionkey;
extern unsigned long int dupcounter;

extern int **curstart;
extern int **curend;
extern int **curwalks;
extern unsigned int num_walks;

extern int EndOrdNum[2];
extern int **endtemplate;
extern int **rendtemplate;
extern int **endtemplate2;
extern int **rendtemplate2;

extern unsigned long int tot_right_endhinges;
extern unsigned long int tot_left_endhinges;
extern unsigned long int num_duplicate_right_endhinges;
extern unsigned long int num_duplicate_left_endhinges;
extern unsigned long int num_tspans;
extern unsigned long int num_duplicate_tspans;

extern int **built_walks_start;
extern int **built_walks_end;
extern int **built_walks_direcs;
extern int num_built_walks;

extern char filename[100];
extern FILE* fp;
extern unsigned long int filetotal;
extern unsigned long int filenum;

// Function prototypes
void set_system_params();
void allocate_globals();
int run_legacy_2sap_sampler(void);
void generate_evectors();
void conv_to_array(void);
void conv_endhinges_to_array(void);
void enterhinge(int i, int j, int side, int (*pointordNum)[3], int curlength);
void leavehinge(int i, int j, int side, int (*pointordNum)[3], int curlength);
void rowedges(int i, int j, int (*pointordNum)[3], int curlength);
void coledges(int i, int j, int (*pointordNum)[3], int curlength);
void recordtemplate( int (*pointordNum)[3]);
void enterendhinge(int i, int j, int curlength);
void leaveendhinge(int i, int j, int curlength);
void endrowedges(int i, int j, int curlength);
void endcoledges(int i, int j, int curlength);
void recordendtemplate();
struct hinge_span *newhinge(unsigned int numberofwalks);
struct endhinge *newendhinge(unsigned int numberofwalks);
void fillreordertemplate(int ledges, int redges);
void printordtemp();
void printreordtemp();
void fillrendtemplate();
void fillendtemplate2();
void fillrendtemplate2();
void printbuiltwalks();
void printbuiltwalks_all();
void add_to_built_walks(unsigned long int secnum, int nth_tspan);
void add_right_endhinge(unsigned long int secnum, int nth_endhinge);
int reverse_direc(int direc);
void check_reachability();
void printtofile();

unsigned long int num_section(int side);
unsigned long int num_section_norder(int side);
unsigned long int num_section_endhinge_norder(void);
unsigned long int num_section_endhinge_norder2(void);
unsigned short int LFlag( int (*)[3]);
unsigned short int LFlag_endhinge(int (*EndOrdNum)[2]);
int *int_vecalloc(int low, int high);
unsigned int *unsgn_vecalloc(int low, int high);
unsigned long int *unsgnlong_vecalloc(int low, int high);

#endif
"""

globals_c = """#include "mc_globals.h"

int L = 2;
int M = 1;
int totalspan = 100;
int samplesize = 100;
int runnum = 1;
unsigned int seednum = 227001;
int max_sections;
int max_tspans;
double dom_evalue;
int maxpolys = 10000;
int vec_length;
int vM;
int vL;
int CS_mode = 0;
int mode = 0;
int ham_check = 0;
char* output_dir = "data/SAPs";
double fval = 0.0;

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

""" + get_lines(610, 764)

legacy_c = """#include "mc_globals.h"

""" + get_lines(328, 569)

builder_c = """#include "mc_globals.h"

""" + get_lines(1663, 1956) + "\n" + get_lines(2117, 2323)

utils_c = """#include "mc_globals.h"

""" + get_lines(1327, 1657) + "\n" + get_lines(1963, 2109) + "\n" + get_lines(2330, 2544) + "\n" + get_lines(2640, 3593) + "\n" + get_lines(3693, 3724)

memory_c = """#include "mc_globals.h"

""" + get_lines(2553, 2635)

validation_c = """#include "mc_globals.h"

""" + get_lines(3595, 3691)

deps_c = """#include "mc_globals.h"

""" + get_lines(3728, 3892)

master_c = """#include "mc_globals.h"

""" + get_lines(571, 607) + "\n" + get_lines(772, 1319)

with open(f"{OUT_DIR}/mc_globals.h", "w") as f: f.write(globals_h)
with open(f"{OUT_DIR}/mc_globals.c", "w") as f: f.write(globals_c)
with open(f"{OUT_DIR}/mc_legacy.c", "w") as f: f.write(legacy_c)
with open(f"{OUT_DIR}/mc_builder.c", "w") as f: f.write(builder_c)
with open(f"{OUT_DIR}/mc_utils.c", "w") as f: f.write(utils_c)
with open(f"{OUT_DIR}/mc_memory.c", "w") as f: f.write(memory_c)
with open(f"{OUT_DIR}/mc_validation.c", "w") as f: f.write(validation_c)
with open(f"{OUT_DIR}/mc_deps.c", "w") as f: f.write(deps_c)
with open(f"{OUT_DIR}/mc_master.c", "w") as f: f.write(master_c)

print("Files generated successfully.")
