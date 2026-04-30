CC = gcc
CFLAGS = -O3 -lm -fopenmp
SRC = src/MASTER_TMcalc.c
OUT = tm_master

MC_SRC = src/MASTER_MCsample.c
MC_OUT = mc_master

.PHONY: all clean tm sampler test verify

all: $(OUT) $(MC_OUT)

tm: $(OUT)

sampler: $(MC_OUT)

$(OUT): $(SRC)
	$(CC) $(SRC) -o $(OUT) $(CFLAGS)

$(MC_OUT): $(MC_SRC)
	$(CC) $(MC_SRC) -o $(MC_OUT) $(CFLAGS)

clean:
	rm -f $(OUT) $(MC_OUT) data/*.txt data/*.bin

test: $(OUT)
	python3 audit_engine.py -L 1 -M 1 -m 0

verify: $(OUT)
	bash verify_all.sh
