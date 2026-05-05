CC = gcc
INCDIR = include
SRCDIR = src
BUILDDIR = build
BINDIR = bin

COMMON_WARNFLAGS = -Wall -Wextra -Wformat=2
STRICT_WARNFLAGS = $(COMMON_WARNFLAGS) -Wsign-compare -Wshadow
LEGACY_WARNFLAGS = $(COMMON_WARNFLAGS) -Wno-sign-compare -Wno-shadow

# Base CFLAGS
CFLAGS_BASE = -O3 -I$(INCDIR) -lm -fopenmp
STRICT_CFLAGS = $(CFLAGS_BASE) $(STRICT_WARNFLAGS)
MC_CFLAGS = $(CFLAGS_BASE) $(LEGACY_WARNFLAGS)
GEN_CFLAGS = $(CFLAGS_BASE) $(LEGACY_WARNFLAGS) -Wno-unused-variable -Wno-unused-parameter -Wno-unused-result -Wno-misleading-indentation -I$(SRCDIR)

# Target executables
TM_OUT = $(BINDIR)/tm_master
MC_OUT = $(BINDIR)/mc_master
CREATOR_OUT = $(BINDIR)/creator_all

# Source files
TM_SRC = $(SRCDIR)/MASTER_TMcalc.c $(SRCDIR)/tm_spectral.c
MC_SRC = $(SRCDIR)/MASTER_MCsample.c $(SRCDIR)/mc_sysparams.c $(SRCDIR)/mc_globals.c \
         $(SRCDIR)/mc_builder.c $(SRCDIR)/mc_utils.c \
         $(SRCDIR)/mc_memory.c $(SRCDIR)/mc_validation.c $(SRCDIR)/mc_deps.c \
         $(SRCDIR)/mc_spectral.c $(SRCDIR)/mc_sampler_weights.c \
         $(SRCDIR)/mc_2sap_integrated.c $(SRCDIR)/mc_2sap_ham_integrated.c
CREATOR_SRC = $(SRCDIR)/MASTER_CreatorAll.c $(SRCDIR)/mc_sysparams.c $(SRCDIR)/mc_globals.c \
              $(SRCDIR)/mc_builder.c $(SRCDIR)/mc_utils.c \
              $(SRCDIR)/mc_memory.c $(SRCDIR)/mc_validation.c $(SRCDIR)/mc_deps.c \
              $(SRCDIR)/mc_spectral.c $(SRCDIR)/mc_2sap_integrated.c $(SRCDIR)/mc_2sap_ham_integrated.c

# Object files
TM_OBJS = $(patsubst $(SRCDIR)/%.c, $(BUILDDIR)/%.tm.o, $(TM_SRC))
MC_OBJS = $(patsubst $(SRCDIR)/%.c, $(BUILDDIR)/%.o, $(MC_SRC))
CREATOR_OBJS = $(patsubst $(SRCDIR)/%.c, $(BUILDDIR)/%.creator.o, $(CREATOR_SRC))

.PHONY: all clean tm sampler creator test verify parity-audit directories

all: directories $(TM_OUT) $(MC_OUT) $(CREATOR_OUT)

directories:
	@mkdir -p $(BUILDDIR) $(BINDIR)

tm: directories $(TM_OUT)

sampler: directories $(MC_OUT)

creator: directories $(CREATOR_OUT)

# TM Master Target
$(TM_OUT): $(TM_OBJS)
	$(CC) $(TM_OBJS) -o $@ $(STRICT_CFLAGS)

$(BUILDDIR)/%.tm.o: $(SRCDIR)/%.c $(INCDIR)/tm_runtime.h $(INCDIR)/tm_spectral.h
	$(CC) -c $< -o $@ $(STRICT_CFLAGS)

# MC Master Target
$(MC_OUT): $(MC_OBJS)
	$(CC) $(MC_OBJS) -o $@ $(MC_CFLAGS)

$(CREATOR_OUT): $(CREATOR_OBJS)
	$(CC) $(CREATOR_OBJS) -o $@ $(MC_CFLAGS)

# Compile MC objects
$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(INCDIR)/mc_globals.h $(INCDIR)/mc_runtime.h $(INCDIR)/mc_spectral.h $(INCDIR)/mc_sampler_weights.h
	$(CC) -c $< -o $@ $(MC_CFLAGS)

$(BUILDDIR)/%.creator.o: $(SRCDIR)/%.c $(INCDIR)/mc_globals.h $(INCDIR)/mc_runtime.h $(INCDIR)/mc_spectral.h
	$(CC) -c $< -o $@ $(MC_CFLAGS)

clean:
	rm -rf $(BUILDDIR) $(BINDIR)
	rm -rf data/TMresults data/MC_Evectors data/CreatorAll
	rm -f data/*.txt data/*.bin

test: directories $(TM_OUT)
	@mkdir -p data/TMresults
	python3 scripts/audit_engine.py -L 1 -M 1 -m 0

verify: directories $(TM_OUT)
	@mkdir -p data/TMresults
	bash scripts/verify_all.sh

parity-audit: all
	python3 scripts/parity_audit.py --no-build
