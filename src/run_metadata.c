#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "run_metadata.h"

#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif

static char metadata_command[2048] = "";

void run_metadata_set_command(int argc, char **argv)
{
	size_t used = 0;
	int i;

	metadata_command[0] = '\0';
	for (i = 0; i < argc; i++) {
		int written = snprintf(
			metadata_command + used,
			sizeof(metadata_command) - used,
			"%s%s",
			i == 0 ? "" : " ",
			argv[i]);
		if (written < 0) {
			metadata_command[0] = '\0';
			return;
		}
		if ((size_t)written >= sizeof(metadata_command) - used) {
			used = sizeof(metadata_command) - 1;
			metadata_command[used] = '\0';
			return;
		}
		used += (size_t)written;
	}
}

void run_metadata_write(
	const char *output_path,
	const char *producer,
	const char *artifact,
	int mode,
	int lattice_l,
	int lattice_m,
	int span,
	unsigned int seed,
	int run,
	double eigenvalue)
{
	char metadata_path[4096];
	time_t now;
	FILE *fp;
#ifdef _OPENMP
	int threads = omp_get_max_threads();
#else
	int threads = 1;
#endif

	if (!output_path || !*output_path) {
		return;
	}
	if (snprintf(metadata_path, sizeof(metadata_path), "%s.meta", output_path) >= (int)sizeof(metadata_path)) {
		fprintf(stderr, "Warning: metadata path too long for %s\n", output_path);
		return;
	}

	fp = fopen(metadata_path, "w");
	if (!fp) {
		fprintf(stderr, "Warning: could not write metadata file %s\n", metadata_path);
		return;
	}

	now = time(NULL);
	fprintf(fp, "producer=%s\n", producer ? producer : "unknown");
	fprintf(fp, "artifact=%s\n", artifact ? artifact : "unknown");
	fprintf(fp, "output=%s\n", output_path);
	fprintf(fp, "command=%s\n", metadata_command[0] ? metadata_command : "unknown");
	fprintf(fp, "git_commit=%s\n", GIT_COMMIT);
	fprintf(fp, "timestamp_unix=%ld\n", (long)now);
	fprintf(fp, "mode=%d\n", mode);
	fprintf(fp, "L=%d\n", lattice_l);
	fprintf(fp, "M=%d\n", lattice_m);
	fprintf(fp, "span=%d\n", span);
	fprintf(fp, "seed=%u\n", seed);
	fprintf(fp, "run=%d\n", run);
	fprintf(fp, "eigenvalue=%.17g\n", eigenvalue);
	fprintf(fp, "openmp_threads=%d\n", threads);
	fclose(fp);
}
