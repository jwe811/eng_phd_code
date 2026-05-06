#ifndef MC_PATHS_H
#define MC_PATHS_H

typedef struct McModePaths {
	const char *sample_dir;
	const char *sample_prefix;
	const char *creator_dir;
	const char *creator_prefix;
	const char *label;
} McModePaths;

const McModePaths *mc_mode_paths(int mode);

#endif
