#ifndef TM_RUNTIME_H
#define TM_RUNTIME_H

static void fatal_error(const char *message) {
    fprintf(stderr, "Fatal: %s\n", message);
    exit(EXIT_FAILURE);
}

static void *xcalloc(size_t count, size_t size, const char *label) {
    void *ptr = calloc(count, size);
    if (!ptr) {
        fprintf(stderr, "Fatal: unable to allocate %s (%zu x %zu bytes)\n", label, count, size);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

static void *xmalloc(size_t size, const char *label) {
    void *ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "Fatal: unable to allocate %s (%zu bytes)\n", label, size);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

static void checked_snprintf(char *buffer, size_t size, const char *fmt, ...) {
    va_list args;
    int written;

    va_start(args, fmt);
    written = vsnprintf(buffer, size, fmt, args);
    va_end(args);

    if (written < 0 || (size_t)written >= size) {
        fatal_error("formatted path exceeded its destination buffer");
    }
}

static void ensure_directory(const char *path) {
    if (mkdir(path, 0775) != 0 && errno != EEXIST) {
        fprintf(stderr, "Fatal: could not create directory '%s': %s\n", path, strerror(errno));
        exit(EXIT_FAILURE);
    }
}

static FILE *xfopen(const char *path, const char *open_mode) {
    FILE *fp = fopen(path, open_mode);
    if (!fp) {
        fprintf(stderr, "Fatal: could not open '%s': %s\n", path, strerror(errno));
        exit(EXIT_FAILURE);
    }
    return fp;
}

static void checked_fwrite(const void *ptr, size_t size, size_t count, FILE *fp, const char *label) {
    if (fwrite(ptr, size, count, fp) != count) {
        fprintf(stderr, "Fatal: failed writing %s\n", label);
        exit(EXIT_FAILURE);
    }
}

#define SECTION_HASH_SIZE 131072
struct section_hash_node {
    unsigned long int hash;
    unsigned long int section_num;
    struct section_hash_node *next;
} *section_hash_table[SECTION_HASH_SIZE];

static unsigned long int next_section_num = 1;

static unsigned long int get_or_add_section(unsigned long int section_hash) {
    unsigned long int bucket = section_hash % SECTION_HASH_SIZE;
    struct section_hash_node *node = section_hash_table[bucket];

    while (node) {
        if (node->hash == section_hash) return node->section_num;
        node = node->next;
    }

    if (next_section_num > actual_max_sections) {
        fprintf(stderr, "Fatal: actual_max_sections overflow (%lu)\n", actual_max_sections);
        exit(EXIT_FAILURE);
    }

    node = (struct section_hash_node *)xmalloc(sizeof(struct section_hash_node), "section hash node");
    node->hash = section_hash;
    node->section_num = next_section_num;
    node->next = section_hash_table[bucket];
    section_hash_table[bucket] = node;

    sectionkey[next_section_num] = section_hash;
    return next_section_num++;
}

#endif
