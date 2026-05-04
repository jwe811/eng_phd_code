#ifndef MC_RUNTIME_H
#define MC_RUNTIME_H

static inline void *mc_xcalloc(size_t count, size_t size, const char *label)
{
	void *ptr = calloc(count, size);
	if (!ptr) {
		fprintf(stderr, "Fatal: unable to allocate %s (%zu x %zu bytes)\n", label, count, size);
		exit(EXIT_FAILURE);
	}
	return ptr;
}

static inline void *mc_xmalloc(size_t size, const char *label)
{
	void *ptr = malloc(size);
	if (!ptr) {
		fprintf(stderr, "Fatal: unable to allocate %s (%zu bytes)\n", label, size);
		exit(EXIT_FAILURE);
	}
	return ptr;
}

#endif
