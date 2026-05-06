#ifndef MC_RUNTIME_H
#define MC_RUNTIME_H

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

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

static inline void mc_checked_snprintf(char *buffer, size_t size, const char *format, ...)
{
	int needed;
	va_list args;

	va_start(args, format);
	needed = vsnprintf(buffer, size, format, args);
	va_end(args);

	if (needed < 0 || (size_t)needed >= size) {
		fprintf(stderr, "Fatal: formatted path exceeded buffer size\n");
		exit(EXIT_FAILURE);
	}
}

#endif
