# ============================================================
# MinEl_PPG — build file
# ============================================================

CC       ?= gcc
CFLAGS   ?= -O3 -std=c11 -D_DEFAULT_SOURCE -DNDEBUG \
            -ffp-contract=off -fno-math-errno -pipe \
            -Wall -Wextra -Wshadow
CFLAGS   += -MMD -MP
LDFLAGS  ?=
LDLIBS   ?= -lm -pthread

# Optional: `make MARCH=native` for local tuning only
MARCH    ?=
ifneq ($(MARCH),)
CFLAGS   += -march=$(MARCH)
endif

SRCS  = numerics.c minel_ppg.c msptd_pb.c
OBJ   = $(SRCS:.c=.o)

BINS  = run_detector bench throughput
DEPS  = $(OBJ:.o=.d) run_detector.d bench.d throughput.d

.PHONY: all clean

all: $(BINS)

run_detector: run_detector.o $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bench: bench.o $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

throughput: throughput.o $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.d $(BINS) *.gcda *.gcno

-include $(DEPS)