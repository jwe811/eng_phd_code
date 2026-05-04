CC = gcc
INCDIR = include
SRCDIR = src
GENDIR = generated
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
MC2SAP_OUT = $(BINDIR)/mc_2sap
MC2SAP_HAM_OUT = $(BINDIR)/mc_2sap_ham

# Source files
TM_SRC = $(SRCDIR)/MASTER_TMcalc.c
MC_SRC = $(SRCDIR)/MASTER_MCsample.c $(SRCDIR)/mc_sysparams.c $(SRCDIR)/mc_globals.c \
         $(SRCDIR)/mc_legacy.c $(SRCDIR)/mc_builder.c $(SRCDIR)/mc_utils.c \
         $(SRCDIR)/mc_memory.c $(SRCDIR)/mc_validation.c $(SRCDIR)/mc_deps.c

# Object files
MC_OBJS = $(patsubst $(SRCDIR)/%.c, $(BUILDDIR)/%.o, $(MC_SRC))

.PHONY: all clean tm sampler test verify directories

all: directories $(TM_OUT) $(MC_OUT) $(MC2SAP_OUT) $(MC2SAP_HAM_OUT)

directories:
	@mkdir -p $(BUILDDIR) $(BINDIR) $(GENDIR)

tm: directories $(TM_OUT)

sampler: directories $(MC_OUT) $(MC2SAP_OUT) $(MC2SAP_HAM_OUT)

# TM Master Target
$(TM_OUT): $(TM_SRC) $(INCDIR)/tm_runtime.h
	$(CC) $(TM_SRC) -o $@ $(STRICT_CFLAGS)

# MC Master Target
$(MC_OUT): $(MC_OBJS)
	$(CC) $(MC_OBJS) -o $@ $(MC_CFLAGS)

# Compile MC objects
$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(INCDIR)/mc_globals.h $(INCDIR)/mc_runtime.h
	$(CC) -c $< -o $@ $(MC_CFLAGS)

# Generated targets
$(GENDIR)/mc_2sap.c $(GENDIR)/mc_2sap_ham.c: deps/archive/monte_carlo/2SAP_MCsample.c deps/archive/monte_carlo/2SAP_MCsample_Ham.c scripts/build_2sap_generic.py
	python3 scripts/build_2sap_generic.py

$(MC2SAP_OUT): $(GENDIR)/mc_2sap.c $(SRCDIR)/mc_sysparams.c
	$(CC) $(GENDIR)/mc_2sap.c $(SRCDIR)/mc_sysparams.c -o $@ -lm $(GEN_CFLAGS)

$(MC2SAP_HAM_OUT): $(GENDIR)/mc_2sap_ham.c $(SRCDIR)/mc_sysparams.c
	$(CC) $(GENDIR)/mc_2sap_ham.c $(SRCDIR)/mc_sysparams.c -o $@ -lm $(GEN_CFLAGS)

clean:
	rm -rf $(BUILDDIR) $(BINDIR)
	rm -rf data/TMresults
	rm -f data/*.txt data/*.bin

test: directories $(TM_OUT)
	@mkdir -p data/TMresults
	python3 scripts/audit_engine.py -L 1 -M 1 -m 0

verify: directories $(TM_OUT)
	@mkdir -p data/TMresults
	bash scripts/verify_all.sh
