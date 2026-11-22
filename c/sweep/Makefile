CC      ?= gcc
ARCH    ?= -march=native
CFLAGS  ?= -O3 $(ARCH) -fopenmp -std=gnu11 -Wall -Wextra
LDFLAGS ?= -fopenmp -lz

SRC = src/main.c src/tables.c src/sweep.c src/run.c src/brute.c src/patterns.c src/greedy.c src/search.c src/bigbottom.c src/bigsplit.c
HDR = src/scnp.h src/bigbottom.h

all: scnp

scnp: $(SRC) $(HDR)
	$(CC) $(CFLAGS) -o $@ $(SRC) $(LDFLAGS)

# Same program without AVX-512 (exercises the generic vector code path).
scnp-generic: $(SRC) $(HDR)
	$(CC) -O3 -march=x86-64-v2 -fopenmp -std=gnu11 -Wall -Wextra -o $@ $(SRC) $(LDFLAGS)

test: scnp scnp-generic
	./tests/run_tests.sh

clean:
	rm -f scnp scnp-generic

.PHONY: all test clean
