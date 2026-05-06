#include <stddef.h>

#include "mc_paths.h"

static const McModePaths mode_paths[] = {
	{"data/MonteCarlo/SAPs", "MCpolys", "All_SAPs", "AllSAPs", "SAPs"},
	{"data/MonteCarlo/HamSAPs", "MCpolysHam", "All_HamSAPs", "AllHamSAPs", "HamSAPs"},
	{"data/MonteCarlo/2SAPs", "MC2SAPs", "All_2SAPs", "All2SAPs", "2SAPs"},
	{"data/MonteCarlo/Ham2SAPs", "MC2SAPsHam", "All_Ham2SAPs", "AllHam2SAPs", "Ham2SAPs"},
};

const McModePaths *mc_mode_paths(int mode)
{
	if (mode < 0 || mode >= (int)(sizeof(mode_paths) / sizeof(mode_paths[0]))) {
		return NULL;
	}
	return &mode_paths[mode];
}
