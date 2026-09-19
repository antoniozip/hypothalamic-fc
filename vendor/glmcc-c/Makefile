# One source, two binary names. `glmcc_fixed` is what scripts/regenerate_glmcc.py,
# scripts/glmcc_positive_control.py and tests/ invoke; `glmcc` is what the older
# control scripts (run_jitter_control.py, run_gapped_ongoing*.py,
# run_glmcc_c_missing.sh) invoke. They were separate sources once, and the stale
# copy kept producing sign-inverted weights after the live one was fixed.

CC      = gcc
CFLAGS  = -O3 -fopenmp -Wall -Wextra
LDFLAGS = -lm
SRC     = glmcc_fixed.c

.PHONY: all clean check

all: glmcc_fixed glmcc

glmcc_fixed: $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

# Same program under the name the older scripts call.
glmcc: $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

# Parity against Kobayashi's Python reference, all four modes.
check: all
	cd ../.. && python3 tests/test_glmcc_c_parity.py

clean:
	rm -f glmcc_fixed glmcc
