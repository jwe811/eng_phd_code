
#ifndef MC2_NAME
#define MC2_NAME(name) name
#endif

#include "../deps/archive/topology/LFlag_norder2.c"
#include "../deps/archive/topology/LFlag_endhinge_norder2.c"

#include "../deps/archive/utils/noncrossing.c"
#define num_section MC2_NAME(num_section_0V)
#include "../deps/archive/sections/cstatenum.c"
#undef num_section

#define num_section MC2_NAME(num_section_6V)
#include "../deps/archive/sections/Num_section_6V.c"
#undef num_section

#define num_section MC2_NAME(num_section_8V)
#include "../deps/archive/sections/Num_section_8V.c"
#undef num_section

#define num_section MC2_NAME(num_section_10V)
#include "../deps/archive/sections/Num_section_10V.c"
#undef num_section

#define num_section MC2_NAME(num_section_12V)
#include "../deps/archive/sections/Num_section_12V.c"
#undef num_section


#define num_section_norder MC2_NAME(num_section_norder_6V)
#define num_section_norder2 MC2_NAME(num_section_norder2_6V)
#define num_section_endhinge_norder MC2_NAME(num_section_endhinge_norder_6V)
#define num_section_endhinge_norder2 MC2_NAME(num_section_endhinge_norder2_6V)
#define num_section_endhinge_norder3 MC2_NAME(num_section_endhinge_norder3_6V)
#define num_section_endhinge_norder4 MC2_NAME(num_section_endhinge_norder4_6V)
#include "../deps/archive/sections/Num_section_6V_nonordered.c"
#include "../deps/archive/sections/Num_section_6V_nonordered2.c"
#include "../deps/archive/sections/Num_section_6V_endhinge_nonordered.c"
#include "../deps/archive/sections/Num_section_6V_endhinge_nonordered2.c"
#undef num_section_norder
#undef num_section_norder2
#undef num_section_endhinge_norder
#undef num_section_endhinge_norder2
#undef num_section_endhinge_norder3
#undef num_section_endhinge_norder4

#define num_section_norder MC2_NAME(num_section_norder_8V)
#define num_section_norder2 MC2_NAME(num_section_norder2_8V)
#define num_section_endhinge_norder MC2_NAME(num_section_endhinge_norder_8V)
#define num_section_endhinge_norder2 MC2_NAME(num_section_endhinge_norder2_8V)
#define num_section_endhinge_norder3 MC2_NAME(num_section_endhinge_norder3_8V)
#define num_section_endhinge_norder4 MC2_NAME(num_section_endhinge_norder4_8V)
#include "../deps/archive/sections/Num_section_8V_nonordered.c"
#include "../deps/archive/sections/Num_section_8V_nonordered2.c"
#include "../deps/archive/sections/Num_section_8V_endhinge_nonordered.c"
#include "../deps/archive/sections/Num_section_8V_endhinge_nonordered2.c"
#undef num_section_norder
#undef num_section_norder2
#undef num_section_endhinge_norder
#undef num_section_endhinge_norder2
#undef num_section_endhinge_norder3
#undef num_section_endhinge_norder4

#define num_section_norder MC2_NAME(num_section_norder_10V)
#define num_section_norder2 MC2_NAME(num_section_norder2_10V)
#define num_section_endhinge_norder MC2_NAME(num_section_endhinge_norder_10V)
#define num_section_endhinge_norder2 MC2_NAME(num_section_endhinge_norder2_10V)
#define num_section_endhinge_norder3 MC2_NAME(num_section_endhinge_norder3_10V)
#define num_section_endhinge_norder4 MC2_NAME(num_section_endhinge_norder4_10V)
#include "../deps/archive/sections/Num_section_10V_nonordered.c"
#include "../deps/archive/sections/Num_section_10V_nonordered2.c"
#include "../deps/archive/sections/Num_section_10V_endhinge_nonordered.c"
#include "../deps/archive/sections/Num_section_10V_endhinge_nonordered2.c"
#undef num_section_norder
#undef num_section_norder2
#undef num_section_endhinge_norder
#undef num_section_endhinge_norder2
#undef num_section_endhinge_norder3
#undef num_section_endhinge_norder4

#define num_section_norder MC2_NAME(num_section_norder_12V)
#define num_section_norder2 MC2_NAME(num_section_norder2_12V)
#define num_section_endhinge_norder MC2_NAME(num_section_endhinge_norder_12V)
#define num_section_endhinge_norder2 MC2_NAME(num_section_endhinge_norder2_12V)
#define num_section_endhinge_norder3 MC2_NAME(num_section_endhinge_norder3_12V)
#define num_section_endhinge_norder4 MC2_NAME(num_section_endhinge_norder4_12V)
#include "../deps/archive/sections/Num_section_12V_nonordered.c"
#include "../deps/archive/sections/Num_section_12V_nonordered2.c"
#include "../deps/archive/sections/Num_section_12V_endhinge_nonordered.c"
#undef num_section_norder
#undef num_section_norder2
#undef num_section_endhinge_norder
#undef num_section_endhinge_norder2
#undef num_section_endhinge_norder3
#undef num_section_endhinge_norder4


unsigned long int MC2_NAME(num_section)(int side) {
    if (M == 0) return MC2_NAME(num_section_0V)(side);
    int v = vM * vL;
    if (v < 7) return MC2_NAME(num_section_6V)(side);
    else if (v < 9) return MC2_NAME(num_section_8V)(side);
    else if (v < 11) return MC2_NAME(num_section_10V)(side);
    else if (v < 13) return MC2_NAME(num_section_12V)(side);
    fprintf(stderr, "Unsupported grid size for section numbering: %d\n", v);
    exit(1);
}

unsigned long int MC2_NAME(num_section_norder)(int side) {
    int v = vM * vL;
    if (v < 7) return MC2_NAME(num_section_norder_6V)(side);
    else if (v < 9) return MC2_NAME(num_section_norder_8V)(side);
    else if (v < 11) return MC2_NAME(num_section_norder_10V)(side);
    else if (v < 13) return MC2_NAME(num_section_norder_12V)(side);
    fprintf(stderr, "Unsupported grid size for non-ordered section numbering: %d\n", v);
    exit(1);
}

unsigned long int MC2_NAME(num_section_endhinge_norder3_12V)(void) {
    fprintf(stderr, "Num_section_12V_endhinge_nonordered2.c is missing from archive!\n");
    exit(1);
}

unsigned long int MC2_NAME(num_section_endhinge_norder4_12V)(void) {
    fprintf(stderr, "Num_section_12V_endhinge_nonordered2.c is missing from archive!\n");
    exit(1);
}

unsigned long int MC2_NAME(num_section_norder2)(int side) {
    int v = vM * vL;
    if (v < 7) return MC2_NAME(num_section_norder2_6V)(side);
    else if (v < 9) return MC2_NAME(num_section_norder2_8V)(side);
    else if (v < 11) return MC2_NAME(num_section_norder2_10V)(side);
    else if (v < 13) return MC2_NAME(num_section_norder2_12V)(side);
    fprintf(stderr, "Unsupported grid size for non-ordered section numbering 2: %d\n", v);
    exit(1);
}

unsigned long int MC2_NAME(num_section_endhinge_norder)(void) {
    int v = vM * vL;
    if (v < 7) return MC2_NAME(num_section_endhinge_norder_6V)();
    else if (v < 9) return MC2_NAME(num_section_endhinge_norder_8V)();
    else if (v < 11) return MC2_NAME(num_section_endhinge_norder_10V)();
    else if (v < 13) return MC2_NAME(num_section_endhinge_norder_12V)();
    fprintf(stderr, "Unsupported grid size for endhinge section numbering: %d\n", v);
    exit(1);
}

unsigned long int MC2_NAME(num_section_endhinge_norder2)(void) {
    int v = vM * vL;
    if (v < 7) return MC2_NAME(num_section_endhinge_norder2_6V)();
    else if (v < 9) return MC2_NAME(num_section_endhinge_norder2_8V)();
    else if (v < 11) return MC2_NAME(num_section_endhinge_norder2_10V)();
    else if (v < 13) return MC2_NAME(num_section_endhinge_norder2_12V)();
    fprintf(stderr, "Unsupported grid size for endhinge section numbering 2: %d\n", v);
    exit(1);
}

unsigned long int MC2_NAME(num_section_endhinge_norder3)(void) {
    int v = vM * vL;
    if (v < 7) return MC2_NAME(num_section_endhinge_norder3_6V)();
    else if (v < 9) return MC2_NAME(num_section_endhinge_norder3_8V)();
    else if (v < 11) return MC2_NAME(num_section_endhinge_norder3_10V)();
    else if (v < 13) return MC2_NAME(num_section_endhinge_norder3_12V)();
    fprintf(stderr, "Unsupported grid size for endhinge section numbering 3: %d\n", v);
    exit(1);
}

unsigned long int MC2_NAME(num_section_endhinge_norder4)(void) {
    int v = vM * vL;
    if (v < 7) return MC2_NAME(num_section_endhinge_norder4_6V)();
    else if (v < 9) return MC2_NAME(num_section_endhinge_norder4_8V)();
    else if (v < 11) return MC2_NAME(num_section_endhinge_norder4_10V)();
    else if (v < 13) return MC2_NAME(num_section_endhinge_norder4_12V)();
    fprintf(stderr, "Unsupported grid size for endhinge section numbering 4: %d\n", v);
    exit(1);
}

#include "../deps/archive/sections/printsection.c"

#include "../deps/archive/utils/int_vecalloc.c"
#include "../deps/archive/utils/unsgn_vecalloc.c"

#include "../deps/archive/utils/unsgnlong_vecalloc.c"

#include "../deps/archive/utils/vecalloc.c"

#include "../deps/archive/utils/matalloc.c"
