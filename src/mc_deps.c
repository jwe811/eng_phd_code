#include "mc_globals.h"

#if defined(CS)
#include "../deps/mc_compat/legacy_topology/LFlag_0.c"		/* include if M = 0  */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the left and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

//#include "../topology/RFlag_0.c"		/* include if M = 0  */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the right and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#include "../deps/mc_compat/legacy_utils/noncrossing.c"	/* include if M = 0 */
/* This function takes arguements (int side, int a, int b, int c, int d) */
/* and returns 1 if the walk connecting a and b does not cross */
/* the walk connecting c and d it returns 0 otherwise*/
#include "../deps/mc_compat/legacy_sections/cstatenum.c"
#else
#include "../deps/mc_compat/legacy_topology/LFlag_norder.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the left and returns 0 otherwise */
/* it also uses the global variable ordertemplate */
#include "../deps/mc_compat/legacy_topology/LFlag_endhinge_norder.c"

#define tspans_outsection t_outsection
#define tspans_nrr t_nrr
#define R_Evector R_Evector_solve
#include "../deps/mc_compat/legacy_transfer_matrix/pw_meth_ts_LRvec_fcheck.c"
#undef tspans_outsection
#undef tspans_nrr
#undef R_Evector

//#include "../topology/RFlag.c"		/* include if M <> 0 */
/* This function takes an arguement (*pointordNum)[3] */
/* it returns 1 if the two-span being considered connects to phi on the right and returns 0 otherwise */
/* it also uses the global variable ordertemplate */

#if vM*vL<7
/* The following versions are for small vM*vL values, but since vM/vL are runtime variables,
   we include all and dispatch at runtime. */
#endif

/* Standard numbering functions */
#define num_section num_section_6V
#include "../deps/mc_compat/legacy_sections/Num_section_6V.c"
#undef num_section

#define num_section num_section_8V
#include "../deps/mc_compat/legacy_sections/Num_section_8V.c"
#undef num_section

#define num_section num_section_10V
#include "../deps/mc_compat/legacy_sections/Num_section_10V.c"
#undef num_section

#define num_section num_section_12V
#include "../deps/mc_compat/legacy_sections/Num_section_12V.c"
#undef num_section

/* Non-ordered numbering functions */
#define num_section_norder num_section_norder_6V
#include "../deps/mc_compat/legacy_sections/Num_section_6V_nonordered.c"
#undef num_section_norder

#define num_section_norder num_section_norder_8V
#include "../deps/mc_compat/legacy_sections/Num_section_8V_nonordered.c"
#undef num_section_norder

#define num_section_norder num_section_norder_10V
#include "../deps/mc_compat/legacy_sections/Num_section_10V_nonordered.c"
#undef num_section_norder

#define num_section_norder num_section_norder_12V
#include "../deps/mc_compat/legacy_sections/Num_section_12V_nonordered.c"
#undef num_section_norder

/* Endhinge numbering functions */
#define num_section_endhinge_norder num_section_endhinge_norder_6V
#define num_section_endhinge_norder2 num_section_endhinge_norder2_6V
#include "../deps/mc_compat/legacy_sections/Num_section_6V_endhinge_nonordered.c"
#undef num_section_endhinge_norder
#undef num_section_endhinge_norder2

#define num_section_endhinge_norder num_section_endhinge_norder_8V
#define num_section_endhinge_norder2 num_section_endhinge_norder2_8V
#include "../deps/mc_compat/legacy_sections/Num_section_8V_endhinge_nonordered.c"
#undef num_section_endhinge_norder
#undef num_section_endhinge_norder2

#define num_section_endhinge_norder num_section_endhinge_norder_10V
#define num_section_endhinge_norder2 num_section_endhinge_norder2_10V
#include "../deps/mc_compat/legacy_sections/Num_section_10V_endhinge_nonordered.c"
#undef num_section_endhinge_norder
#undef num_section_endhinge_norder2

#define num_section_endhinge_norder num_section_endhinge_norder_12V
#define num_section_endhinge_norder2 num_section_endhinge_norder2_12V
#include "../deps/mc_compat/legacy_sections/Num_section_12V_endhinge_nonordered.c"
#undef num_section_endhinge_norder
#undef num_section_endhinge_norder2

/* Dispatcher wrappers */
unsigned long int num_section(int side) {
    int v = vM * vL;
    if (v < 7) return num_section_6V(side);
    else if (v < 9) return num_section_8V(side);
    else if (v < 11) return num_section_10V(side);
    else if (v < 13) return num_section_12V(side);
    fprintf(stderr, "Unsupported grid size for section numbering: %d\n", v);
    exit(1);
}

unsigned long int num_section_norder(int side) {
    int v = vM * vL;
    if (v < 7) return num_section_norder_6V(side);
    else if (v < 9) return num_section_norder_8V(side);
    else if (v < 11) return num_section_norder_10V(side);
    else if (v < 13) return num_section_norder_12V(side);
    fprintf(stderr, "Unsupported grid size for non-ordered section numbering: %d\n", v);
    exit(1);
}

unsigned long int num_section_endhinge_norder(void) {
    int v = vM * vL;
    if (v < 7) return num_section_endhinge_norder_6V();
    else if (v < 9) return num_section_endhinge_norder_8V();
    else if (v < 11) return num_section_endhinge_norder_10V();
    else if (v < 13) return num_section_endhinge_norder_12V();
    fprintf(stderr, "Unsupported grid size for endhinge section numbering: %d\n", v);
    exit(1);
}

unsigned long int num_section_endhinge_norder2(void) {
    int v = vM * vL;
    if (v < 7) return num_section_endhinge_norder2_6V();
    else if (v < 9) return num_section_endhinge_norder2_8V();
    else if (v < 11) return num_section_endhinge_norder2_10V();
    else if (v < 13) return num_section_endhinge_norder2_12V();
    fprintf(stderr, "Unsupported grid size for endhinge section numbering 2: %d\n", v);
    exit(1);
}

#endif



#include "../deps/mc_compat/legacy_sections/printsection.c"

#include "../deps/mc_compat/legacy_utils/int_vecalloc.c"
#include "../deps/mc_compat/legacy_utils/unsgn_vecalloc.c"

#include "../deps/utils/unsgnlong_vecalloc.c" /*unsigned long int unsgnlong_vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with unsigned long integer enteries */

#include "../deps/utils/vecalloc.c"		/*vec_ent * vecalloc(int low, int high)*/
/* used for dynamically allocating memory for vectors with enteries of type vec_ent */

#include "../deps/utils/matalloc.c"		/*mat_ent **matalloc(int rowlow, int rowhigh, int collow, int colhigh)*/
/* used for dynamically allocating memory for matrices with enteries of type mat_ent */





	
