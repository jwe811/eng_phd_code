CC = gcc
CFLAGS = -O3 -lm -fopenmp
SRC = src/MASTER_TMcalc.c
OUT = tm_master

MC_SRC = src/mc_master.c src/mc_sysparams.c src/mc_globals.c src/mc_legacy.c src/mc_builder.c src/mc_utils.c src/mc_memory.c src/mc_validation.c src/mc_deps.c
MC_OUT = mc_master

.PHONY: all clean tm sampler test verify

all: $(OUT) $(MC_OUT) mc_2sap mc_2sap_ham

tm: $(OUT)

sampler: $(MC_OUT) mc_2sap mc_2sap_ham

$(OUT): $(SRC)
	$(CC) $(SRC) -o $(OUT) $(CFLAGS)

$(MC_OUT): $(MC_SRC)
	$(CC) $(MC_SRC) -o $(MC_OUT) $(CFLAGS)

src/mc_2sap.c src/mc_2sap_ham.c: src/archive_deps/monte_carlo/2SAP_MCsample.c src/archive_deps/monte_carlo/2SAP_MCsample_Ham.c scripts/build_2sap_generic.py
	python3 scripts/build_2sap_generic.py

mc_2sap: src/mc_2sap.c src/mc_sysparams.c
	$(CC) src/mc_2sap.c src/mc_sysparams.c -o mc_2sap -lm -w -O3

mc_2sap_ham: src/mc_2sap_ham.c src/mc_sysparams.c
	$(CC) src/mc_2sap_ham.c src/mc_sysparams.c -o mc_2sap_ham -lm -w -O3

clean:
	rm -f $(OUT) $(MC_OUT) mc_2sap mc_2sap_ham data/*.txt data/*.bin
	rm -f $(OUT) $(MC_OUT) data/*.txt data/*.bin

test: $(OUT)
	python3 audit_engine.py -L 1 -M 1 -m 0

verify: $(OUT)
	bash verify_all.sh
