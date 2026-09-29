# Host build + regression test of the JDPS C module against the Python reference.
#   make test        # export vectors with Python, build, run
CC      ?= gcc
CFLAGS  ?= -std=c99 -O2 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion
PYTHON  ?= python3

all: cs_jdps_test

cs_jdps_test: cs_jdps.c cs_jdps.h test/cs_jdps_test.c
	$(CC) $(CFLAGS) -D_DEFAULT_SOURCE -o $@ cs_jdps.c test/cs_jdps_test.c -lm

test/vectors.txt: ../cs_agc/algorithms.py ../cs_agc/export_vectors.py
	cd .. && $(PYTHON) -m cs_agc.export_vectors c/test/vectors.txt

test: cs_jdps_test test/vectors.txt
	./cs_jdps_test test/vectors.txt

clean:
	rm -f cs_jdps_test test/vectors.txt

.PHONY: all test clean
