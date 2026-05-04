import os
import re

def process_file(src_path, dst_path, is_ham):
    with open(src_path, 'r') as f:
        content = f.read()

    # 1. Remove #define L, M blocks
    content = re.sub(r'//Automatically set max_sections.*?#define\s+vec_length\s+\(max_sections\+1\)', '', content, flags=re.DOTALL)
    
    # 2. Inject constants and dynamic CLI parser
    fixed_constants = """
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

extern int L;
extern int M;
extern int totalspan;
extern int samplesize;
extern int runnum;
extern unsigned int seednum;
extern int maxpolys;
extern int mode;
extern int ham_check;

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

static void *generated_xcalloc(size_t count, size_t size, const char *label)
{
    void *ptr = calloc(count, size);
    if (ptr == NULL) {
        fprintf(stderr, "Fatal: unable to allocate %s (%zu x %zu bytes)\\n", label, count, size);
        exit(EXIT_FAILURE);
    }
    return ptr;
}


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
    content = content.replace(
        'struct endhinge *nextnewendhinge;\n\n\tnextnewendhinge =',
        'struct endhinge *nextnewendhinge;\n\tunsigned int max_walks = numberofwalks > numberofwalks2 ? numberofwalks : numberofwalks2;\n\n\tnextnewendhinge =',
    )
    content = content.replace('nextnewendhinge->start2[i] = int_vecalloc(0,numberofwalks2-1);', 'nextnewendhinge->start2[i] = int_vecalloc(0,max_walks-1);')
    content = content.replace('nextnewendhinge->end2[i] = int_vecalloc(0,numberofwalks2-1);', 'nextnewendhinge->end2[i] = int_vecalloc(0,max_walks-1);')
    content = content.replace('nextnewendhinge->walks2[i] = int_vecalloc(0,numberofwalks-1);', 'nextnewendhinge->walks2[i] = int_vecalloc(0,max_walks-1);')
    for table_name in ("Lend", "Rend"):
        for suffix in ("", "2"):
            content = content.replace(
                f'(int*)malloc({table_name}_num_walks{suffix}[section_num][i] * sizeof(int))',
                f'(int*)generated_xcalloc(MAX_vMvL + 1, sizeof(int), "{table_name} walk data{suffix}")',
            )

    # 4. Topology Dispatcher
    # Replace the massive #if defined(CS) ... block at the end with a simple include to mc_deps_2sap.c
    content = re.sub(r'#if defined\(CS\)\n#include "\.\./topology/LFlag_0\.c".*?#include "\.\./utils/matalloc\.c"', '#include "mc_deps_2sap.c"', content, flags=re.DOTALL)
    
    # 5. Fix includes
    content = content.replace('#include "../sections/', '#include "../deps/archive/sections/')
    content = content.replace('#include "../topology/', '#include "../deps/archive/topology/')
    content = content.replace('#include "../utils/', '#include "../deps/archive/utils/')
    content = content.replace('#include "../analysis/', '#include "../deps/archive/analysis/')
    content = content.replace('#include "../../include/marsaglia.h"', '#include "../include/marsaglia.h"')
    content = content.replace(
        "ACTUAL 'VALID' SECTIONS FOR THE L=%d, M=%d COMPRESSED CASE IS: %d-%lu=%lu=%lu",
        "ACTUAL 'VALID' SECTIONS FOR THE L=%d, M=%d COMPRESSED CASE IS: %lu-%lu=%lu=%lu",
    )
    content = content.replace('sprintf(filename2, "2SAP_R_Evector_TS_L%dM%d.txt", L, M);', 'sprintf(filename2, "data/MC_Evectors/2SAP_R_Evector_TS_L%dM%d.txt", L, M);')
    content = content.replace('sprintf(filename2, "2SAP_R_EvectorHam_TS_L%dM%d.txt", L, M);', 'sprintf(filename2, "data/MC_Evectors/2SAP_R_EvectorHam_TS_L%dM%d.txt", L, M);')

    if is_ham:
        content += '\n#define tspans_nrr t_nrr\n#define tspans_edges MC_tspans_edges\n#define L_Evector MC_L_Evector\n#define R_Evector MC_R_Evector\n#define tspans_outsection t_outsection\n#include "../deps/archive/transfer_matrix/pw_meth_ts_LRvec_fcheck_2SAP_HAM.c"\n#undef tspans_nrr\n#undef tspans_edges\n#undef L_Evector\n#undef R_Evector\n#undef tspans_outsection\n'
    else:
        content += '\n#define tspans_nrr t_nrr\n#define tspans_edges MC_tspans_edges\n#define L_Evector MC_L_Evector\n#define R_Evector MC_R_Evector\n#define tspans_outsection t_outsection\n#include "../deps/archive/transfer_matrix/pw_meth_ts_LRvec_fcheck_2SAP.c"\n#undef tspans_nrr\n#undef tspans_edges\n#undef L_Evector\n#undef R_Evector\n#undef tspans_outsection\n'

    # 6. Main and Allocations
    generated_mode = 3 if is_ham else 2
    generated_ham_check = 1 if is_ham else 0
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
    mode = __GENERATED_MODE__;
    ham_check = __GENERATED_HAM_CHECK__;
    set_system_params();
    max_keynum = MAX_KEYNUM_ARR - 1;
    """
    cli_parser = cli_parser.replace("__GENERATED_MODE__", str(generated_mode))
    cli_parser = cli_parser.replace("__GENERATED_HAM_CHECK__", str(generated_ham_check))
    content = content.replace('int main(void)', 'int main(int argc, char *argv[])')
    
    dynamic_allocs = cli_parser + """
    first_hinge_span = (struct hinge_span**)generated_xcalloc(max_sections + 1, sizeof(struct hinge_span*), "first hinge span table");
    current_hinge_span = (struct hinge_span**)generated_xcalloc(max_sections + 1, sizeof(struct hinge_span*), "current hinge span table");
    firstendhinge = (struct endhinge**)generated_xcalloc(max_sections + 1, sizeof(struct endhinge*), "first endhinge table");
    currentendhinge = (struct endhinge**)generated_xcalloc(max_sections + 1, sizeof(struct endhinge*), "current endhinge table");
    MC_L_Evector[0] = (double*)generated_xcalloc(max_tspans+1, sizeof(double), "left eigenvector 0");
    MC_L_Evector[1] = (double*)generated_xcalloc(max_tspans+1, sizeof(double), "left eigenvector 1");
    MC_R_Evector[0] = (double*)generated_xcalloc(max_tspans+1, sizeof(double), "right eigenvector 0");
    MC_R_Evector[1] = (double*)generated_xcalloc(max_tspans+1, sizeof(double), "right eigenvector 1");
    MC_tspans_edges = (unsigned long int**)generated_xcalloc(max_keynum+1, sizeof(unsigned long int*), "2SAP edge table");
"""
    content = re.sub(r'int\s+main\s*\([^)]*\)\s*\{', 'int main(int argc, char *argv[])\n{\n' + dynamic_allocs.replace('\\', '\\\\'), content)

    evector_alloc = """
    for(int k=1; k<=max_keynum; k++){
        if(num_outsections[k] == 0) continue;
        MC_tspans_edges[k] = unsgnlong_vecalloc(1, num_outsections[k]);
        for(int m=1; m<=num_outsections[k]; m++) MC_tspans_edges[k][m] = 0;
    }
    double calculated_dom_evalue = max_eval_LRvec(1.0) + 1.0;
    dom_evalue = calculated_dom_evalue;

    mkdir("data/MC_Evectors", 0775);
    char export_fn[128];
    sprintf(export_fn, "data/MC_Evectors/2SAP_R_Evector%s_TS_L%dM%d.txt", (ham_check ? "Ham" : ""), L, M);
    FILE *export_fp = fopen(export_fn, "w");
    if (export_fp != NULL) {
        for (int i = 1; i <= max_tspans; i++) {
            fprintf(export_fp, "%.15f\\n", MC_R_Evector[0][i]);
        }
        fclose(export_fp);
        printf("Calculated eigenvectors exported to %s\\n", export_fn);
    }
"""
    content = content.replace('double* R_Evector;', 'double* R_Evector;\n' + evector_alloc)
    content = re.sub(
        r'\n\tR_Evector = \(double\*\)malloc\(sizeof\(double\)\*\(max_tspans\+1\)\);\n'
        r'\tif\(R_Evector==NULL\)\{\n'
        r'\t\tfprintf\(stderr, "Out of memory"\);\n'
        r'\t\texit\(0\);\n'
        r'\t\}\n'
        r'.*?'
        r'\tfor\(i=1; i<= max_tspans; i\+\+\)\{\n'
        r'.*?'
        r'\t\}\n',
        '\n\tR_Evector = MC_R_Evector[0];\n',
        content,
        count=1,
        flags=re.DOTALL,
    )

    with open(dst_path, 'w') as f:
        f.write(content)

os.makedirs('generated', exist_ok=True)
process_file('deps/archive/monte_carlo/2SAP_MCsample.c', 'generated/mc_2sap.c', False)
process_file('deps/archive/monte_carlo/2SAP_MCsample_Ham.c', 'generated/mc_2sap_ham.c', True)