CC = gcc
CFLAGS = -Wall -O2 -std=c99
LDFLAGS = -lm

# Main programs
PROGRAMS = ising1d ising2d ising3d ising2d_fss ising3d_fss
TEST_PROGRAM = test_suite

# Object files
OBJS =

all: $(PROGRAMS) $(TEST_PROGRAM)

ising1d: ising1d.c pcg_random.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

ising2d: ising2d.c pcg_random.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

ising3d: ising3d.c pcg_random.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

ising2d_fss: ising2d_fss.c pcg_random.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

ising3d_fss: ising3d_fss.c pcg_random.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

test_suite: test_suite.c pcg_random.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

test: test_suite
	./test_suite

clean:
	rm -f $(PROGRAMS) $(TEST_PROGRAM) *.o

# Quick test runs
test-1d:
	@echo "Running 1D Ising model (N=20) for first 10 temperature points:"
	@./ising1d 20 | head -10

test-2d:
	@echo "Running 2D Ising model (N=5) for first 5 temperature points:"
	@timeout 30s ./ising2d 5 | head -5 || echo "Test completed (may have timed out)"

test-3d:
	@echo "Running 3D Ising model (N=5) for first 5 temperature points:"
	@timeout 30s ./ising3d 5 | head -5 || echo "Test completed (may have timed out)"

benchmark:
	@echo "=== Performance Benchmark ==="
	@echo "1D Ising (N=100):"
	@time -p ./ising1d 100 > /dev/null
	@echo "2D Ising (N=10):"
	@time -p ./ising2d 10 | head -1 > /dev/null

help:
	@echo "Available targets:"
	@echo "  all       - Build all programs"
	@echo "  test      - Run test suite"
	@echo "  test-1d   - Quick test of 1D model"
	@echo "  test-2d   - Quick test of 2D model"
	@echo "  test-3d   - Quick test of 3D model"
	@echo "  benchmark - Performance benchmark"
	@echo "  clean     - Remove executables"
	@echo "  help      - Show this help"

.PHONY: all test clean test-1d test-2d test-3d benchmark help