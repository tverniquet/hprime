#
# Kind of a brute force makefile that compiles everything together
#
TARGET = hprime
TARGET_DEBUG = hprime-debug

DIRS= src/common \
      src/calc_blocks \
      src/prime

SRC=$(foreach mydir, $(DIRS), $(shell find $(mydir)/ -name '*.c' | grep -v 'unused'))

HEADERS=$(foreach mydir, $(DIRS), $(shell find $(mydir)/ -name '*.h'))
HEADERDIRS = $(shell find src/ -type d)

RELEASE_PROG_MAIN = src/prog/main.c
DEBUG_PROG_MAIN = src/prog/test.c

LDFLAGS = -lm -lrt -lpthread
CC  = gcc-5
CFLAGS = -D_ISOC11_SOURCE -Wall -Wextra $(foreach var, $(HEADERDIRS), -I $(var))
OPTIMISE =  -O3 -march=native
OPTIMISE_DEBUG =  -march=native

.PHONY: release debug all clean dummy

all: release

release: bin/$(TARGET)

debug:   bin/$(TARGET_DEBUG)

bin/$(TARGET_DEBUG): $(SRC) $(DEBUG_PROG_MAIN) $(HEADERS) Makefile
	@mkdir -p $(dir $@)
	$(CC) -g -o $@ $(CFLAGS) $(OPTIMISE_DEBUG) $(SRC) $(DEBUG_PROG_MAIN) $(LDFLAGS)

bin/$(TARGET): $(SRC) $(RELEASE_PROG_MAIN) $(HEADERS) Makefile
	@mkdir -p $(dir $@)
	$(CC) -g -o $@ $(CFLAGS) $(OPTIMISE) $(SRC) $(RELEASE_PROG_MAIN) $(LDFLAGS)

# Quick way to run some benchmarks Eg make release run
run:
	bin/$(TARGET) 0 1000000000
	bin/$(TARGET) 0 10000000000

# Test the newest plan against the previous (working) plan
test: bin/$(TARGET_DEBUG)
	bin/$(TARGET_DEBUG) 0 1000000000 1 0

clean:
	-rm -rf bin/$(TARGET) bin/$(TARGET_DEBUG)

