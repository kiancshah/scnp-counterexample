CC      ?= gcc
CFLAGS  ?= -O3 -march=native -fopenmp -std=gnu11 -Wall -Wextra -Wno-maybe-uninitialized

all: scnp_fast scnp_fast_generic scnp_fast12

scnp_fast: scnp_fast.c
	$(CC) $(CFLAGS) -o $@ $<

# Same program without AVX-512 (exercises the portable code path).
scnp_fast_generic: scnp_fast.c
	$(CC) -O3 -march=x86-64-v2 -fopenmp -std=gnu11 -Wall -Wextra -o $@ $<

clean:
	rm -f scnp_fast scnp_fast_generic scnp_fast12

.PHONY: all clean

# n = 12 build (128-byte H vectors; steps generated on the fly)
scnp_fast12: scnp_fast.c
	$(CC) $(CFLAGS) -DMAXN=12 -o $@ $<
