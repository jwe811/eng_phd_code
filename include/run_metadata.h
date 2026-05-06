#ifndef RUN_METADATA_H
#define RUN_METADATA_H

void run_metadata_set_command(int argc, char **argv);
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
	double eigenvalue);

#endif
