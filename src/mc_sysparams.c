#include "mc_globals.h"

int L = 2;
int M = 1;
int totalspan = 100;
int samplesize = 100;
int runnum = 1;
unsigned int seednum = 227001;
int max_sections;
unsigned long int max_keynum;
int max_tspans;
double dom_evalue;
int maxpolys = 10000;
int vec_length;
int vM;
int vL;
int CS_mode = 0;
int mode = 0;
int ham_check = 0;
char* output_dir = "data/MonteCarlo/SAPs";
double fval = 0.0;

void set_system_params() {
    vM = M + 1;
    vL = L + 1;
    if (M == 0) CS_mode = 1;
    
    if (ham_check) {
        if (M == 1 && L == 1) {
            max_sections = 8;
            max_tspans = 22;
            dom_evalue = 3.732050810014727;
        } else if (M == 1 && L == 2) {
            max_sections = 73;
            max_tspans = 649;
            dom_evalue = 14.076438172824952;
        } else if (M == 1 && L == 3) {
            max_sections = 742;
            max_tspans = 19554;
            dom_evalue = 49.643407510907970;
        } else if (M == 1 && L == 4) {
            max_sections = 9309;
            max_tspans = 728925;
            dom_evalue = 172.714480164060546;
        } else if (M == 1 && L == 5) {
            max_sections = 138038;
            max_tspans = 32294131;
            dom_evalue = 596.477267188525502;
        } else if (M == 2 && L == 2) {
            max_sections = 2619;
            max_tspans = 184574;
            dom_evalue = 104.488979151088131;
        } else if (M == 2 && L == 3) {
            max_sections = 138322;
            max_tspans = 39671908;
            dom_evalue = 880.559222610926781;
        } else {
            fprintf(stderr, "Unsupported Hamiltonian L and M values (%d, %d).\n", L, M);
            exit(1);
        }
    } else {
        if (M == 1 && L == 1) {
            max_sections = 8;
            max_tspans = 48;
            dom_evalue = 7.0;
        } else if (M == 1 && L == 2) {
            max_sections = 73;
            max_tspans = 1829;
            dom_evalue = 34.360180657516501;
        } else if (M == 1 && L == 3) {
            max_sections = 742;
            max_tspans = 70306;
            dom_evalue = 165.169003076277392;
        } else if (M == 1 && L == 4) {
            max_sections = 9309;
            max_tspans = 3165653;
            dom_evalue = 770.927136914359608;
        } else if (M == 1 && L == 5) {
            max_sections = 138038;
            max_tspans = 165637127;
            dom_evalue = 3562.911364811740896;
        } else if (M == 2 && L == 2) {
            max_sections = 2619;
            max_tspans = 513585;
            dom_evalue = 416.870158713340345;
        } else if (M == 2 && L == 3) {
            max_sections = 138322;
            max_tspans = 201423784;
            dom_evalue = 4945.864451078689854;
        } else {
            fprintf(stderr, "Unsupported Standard L and M values (%d, %d).\n", L, M);
            exit(1);
        }
    }
    vec_length = max_sections + 1;
    max_keynum = max_sections;
}
