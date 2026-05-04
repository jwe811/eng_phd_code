import os
import re

def process_file(src_path, dst_path, is_ham):
    with open(src_path, 'r') as f:
        content = f.read()

    # 1. Remove #define L, M blocks
    content = re.sub(r'//Automatically set max_sections.*?#define\s+vec_length\s+\(max_sections\+1\)', '', content, flags=re.DOTALL)
    
    # 2. Inject constants and dynamic CLI parser
    fixed_constants = """
#include <getopt.h>

extern int L;
extern int M;
extern int totalspan;
extern int samplesize;
extern int runnum;
extern unsigned int seednum;
extern int maxpolys;

extern int max_sections;
extern int max_tspans;
extern unsigned long int max_keynum;
extern double dom_evalue;
extern double fval;

extern void set_system_params();
extern double max_eval_LRvec(double fugacity);

double *MC_L_Evector[2];
double *MC_R_Evector[2];
unsigned long int **MC_tspans_edges;


#define MAX_vM 5
#define MAX_vL 5
#define MAX_vMvL 25
#define MAX_SPAN 20
#define MAX_KEYNUM_ARR 140000

#define vM (M+1)
#define vL (L+1)
"""
    content = content.replace('#include <stdio.h>', '#include <stdio.h>\n' + fixed_constants)

    # 3. Array bounds
    content = re.sub(r'\[vM\*vL\]', '[MAX_vMvL]', content)
    content = re.sub(r'\[vM\*vL\+1\]', '[MAX_vMvL+1]', content)
    content = re.sub(r'\[vM\*vL/2\]', '[MAX_vMvL/2]', content)
    content = re.sub(r'\[vM\*vL\*\(totalspan\+1\)\]', '[MAX_vMvL*(MAX_SPAN+1)]', content)
    content = re.sub(r'\[vM\]', '[MAX_vM]', content)
    content = re.sub(r'\[vL\]', '[MAX_vL]', content)
    content = re.sub(r'\[M\]', '[MAX_vM]', content)
    content = content.replace('[L]', '[MAX_vL]')
    content = content.replace('[max_keynum+1]', '[MAX_KEYNUM_ARR+1]')
    content = content.replace('[max_keynum]', '[MAX_KEYNUM_ARR]')
    content = content.replace('[vec_length]', '[MAX_KEYNUM_ARR+1]')

    content = re.sub(r'#\s*define\s+L\s+\d+', '', content)
    content = re.sub(r'#\s*define\s+M\s+\d+', '', content)
    content = re.sub(r'#\s*define\s+totalspan\s+\d+', '', content)
    content = re.sub(r'#\s*define\s+samplesize\s+\d+', '', content)
    content = re.sub(r'#\s*define\s+runnum\s+\d+', '', content)
    content = re.sub(r'#\s*define\s+seednum\s+\d+', '', content)
    content = re.sub(r'#\s*define\s+maxpolys\s+\d+', '', content)

    content = content.replace('}	*first_hinge_span[MAX_KEYNUM_ARR+1], *current_hinge_span[MAX_KEYNUM_ARR+1];', '} **first_hinge_span, **current_hinge_span;')
    content = content.replace('struct endhinge* firstendhinge[MAX_KEYNUM_ARR+1];', 'struct endhinge **firstendhinge;')
    content = content.replace('struct endhinge* currentendhinge[MAX_KEYNUM_ARR+1];', 'struct endhinge **currentendhinge;')
    content = content.replace('struct endhinge *firstendhinge[MAX_KEYNUM_ARR+1], *currentendhinge[MAX_KEYNUM_ARR+1];', 'struct endhinge **firstendhinge, **currentendhinge;')

    # 4. Topology Dispatcher
    # Replace the massive #if defined(CS) ... block at the end with a simple include to mc_deps_2sap.c
    content = re.sub(r'#if defined\(CS\)\n#include "\.\./topology/LFlag_0\.c".*?#include "\.\./utils/matalloc\.c"', '#include "mc_deps_2sap.c"', content, flags=re.DOTALL)
    
    # 5. Fix includes
    content = content.replace('#include "../sections/', '#include "archive_deps/sections/')
    content = content.replace('#include "../topology/', '#include "archive_deps/topology/')
    content = content.replace('#include "../utils/', '#include "archive_deps/utils/')
    content = content.replace('#include "../analysis/', '#include "archive_deps/analysis/')
    content = content.replace('#include "../../include/marsaglia.h"', '#include "../include/marsaglia.h"')
    content = content.replace(
        "ACTUAL 'VALID' SECTIONS FOR THE L=%d, M=%d COMPRESSED CASE IS: %d-%lu=%lu=%lu",
        "ACTUAL 'VALID' SECTIONS FOR THE L=%d, M=%d COMPRESSED CASE IS: %lu-%lu=%lu=%lu",
    )

    if is_ham:
        content += '\n#define tspans_nrr t_nrr\n#define tspans_edges MC_tspans_edges\n#define L_Evector MC_L_Evector\n#define R_Evector MC_R_Evector\n#define tspans_outsection t_outsection\n#include "archive_deps/transfer_matrix/pw_meth_ts_LRvec_fcheck_2SAP_HAM.c"\n#undef tspans_nrr\n#undef tspans_edges\n#undef L_Evector\n#undef R_Evector\n#undef tspans_outsection\n'
    else:
        content += '\n#define tspans_nrr t_nrr\n#define tspans_edges MC_tspans_edges\n#define L_Evector MC_L_Evector\n#define R_Evector MC_R_Evector\n#define tspans_outsection t_outsection\n#include "archive_deps/transfer_matrix/pw_meth_ts_LRvec_fcheck_2SAP.c"\n#undef tspans_nrr\n#undef tspans_edges\n#undef L_Evector\n#undef R_Evector\n#undef tspans_outsection\n'

    # 6. Main and Allocations
    cli_parser = """
    int opt;
    while ((opt = getopt(argc, argv, "L:M:s:n:r:S:")) != -1) {
        switch (opt) {
            case 'L': L = atoi(optarg); break;
            case 'M': M = atoi(optarg); break;
            case 's': totalspan = atoi(optarg); break;
            case 'n': samplesize = atoi(optarg); break;
            case 'r': runnum = atoi(optarg); break;
            case 'S': seednum = (unsigned int)atoi(optarg); break;
        }
    }
    set_system_params();
    """
    content = content.replace('int main(void)', 'int main(int argc, char *argv[])')
    
    dynamic_allocs = cli_parser + """
    first_hinge_span = (struct hinge_span**)malloc(sizeof(struct hinge_span*) * (max_sections + 1));
    current_hinge_span = (struct hinge_span**)malloc(sizeof(struct hinge_span*) * (max_sections + 1));
    firstendhinge = (struct endhinge**)malloc(sizeof(struct endhinge*) * (max_sections + 1));
    currentendhinge = (struct endhinge**)malloc(sizeof(struct endhinge*) * (max_sections + 1));
    MC_L_Evector[0] = (double*)calloc(max_tspans+1, sizeof(double));
    MC_L_Evector[1] = (double*)calloc(max_tspans+1, sizeof(double));
    MC_R_Evector[0] = (double*)calloc(max_tspans+1, sizeof(double));
    MC_R_Evector[1] = (double*)calloc(max_tspans+1, sizeof(double));
    MC_tspans_edges = (unsigned long int**)malloc(sizeof(unsigned long int*)*(max_keynum+1));
    for(int k=1; k<=max_keynum; k++){
        MC_tspans_edges[k] = unsgnlong_vecalloc(1, num_outsections[k]);
        for(int m=1; m<=num_outsections[k]; m++) MC_tspans_edges[k][m] = 0;
    }
"""
    content = re.sub(r'int\s+main\s*\([^)]*\)\s*\{', 'int main(int argc, char *argv[])\n{\n' + dynamic_allocs.replace('\\', '\\\\'), content)

    evector_alloc = """
    double calculated_dom_evalue = max_eval_LRvec(1.0) + 1.0;
    double* R_Evector_ptr = R_Evector;
"""
    content = content.replace('double* R_Evector;', 'double* R_Evector;\n' + evector_alloc)

    with open(dst_path, 'w') as f:
        f.write(content)

os.makedirs('generated', exist_ok=True)
process_file('src/archive_deps/monte_carlo/2SAP_MCsample.c', 'generated/mc_2sap.c', False)
process_file('src/archive_deps/monte_carlo/2SAP_MCsample_Ham.c', 'generated/mc_2sap_ham.c', True)
