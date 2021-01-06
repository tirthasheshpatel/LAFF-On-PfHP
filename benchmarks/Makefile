.PHONY: all dirs benchmark_level_1_BLAS_copy benchmark_level_1_BLAS_scal benchmark_level_1_BLAS_axpy benchmark_level_1_BLAS_dot benchmark_level_1_BLAS_nrm2 benchmark_level_1_BLAS_asum benchmark_level_2_BLAS_gemv benchmark_level_2_BLAS_symv benchmark_level_2_BLAS_trsv benchmark_level_2_BLAS_ger benchmark_level_2_BLAS_syr benchmark_level_2_BLAS_syr2 benchmark_level_1_BLAS benchmark_level_2_BLAS clean

# Make sure you have BLIS installed in your home directory (or some other BLAS library).
BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created.
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O3 -I$(BLAS_INC) -m64 -mavx2 -std=c99 -march=native -fopenmp -D_POSIX_C_SOURCE=200809L
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed.
NREPEATS   := 3#         number of times each experiment is repeated.  The best time is reported.
NFIRST     := 48#        smallest size to be timed.
NLAST_SMALL:= 100#       largest size to be timed for slow implementations.
NLAST      := 2000#      largest size to be timed for fast implementations.
NINC       := 50#        increment between sizes.

LDFLAGS    := -lpthread -m64 -lm -fopenmp

UTIL_OBJS_LEVEL_1  := FLA_Clock.o RandomMatrix.o
UTIL_OBJS_LEVEL_2  := FLA_Clock.o RandomMatrix.o RandomTril.o RandomSym.o

# benchmark all the BLAS operations.
all: dirs benchmark_level_1_BLAS benchmark_level_2_BLAS

# create directories to store the generated data
dirs:
	@mkdir -p out plots

# compile the level-1 driver.
driver_level_1_BLAS.x: level_1_BLAS.o $(UTIL_OBJS_LEVEL_1)
	@echo "Compiling the driver for Level-1 BLAS operations."
	@$(LINKER) level_1_BLAS.o $(UTIL_OBJS_LEVEL_1) $(BLAS_LIB) -o driver_level_1_BLAS.x $(LDFLAGS)

# compile the level-2 driver.
driver_level_2_BLAS.x: level_2_BLAS.o $(UTIL_OBJS_LEVEL_2)
	@echo "Compiling the driver for Level-2 BLAS operations."
	@$(LINKER) level_2_BLAS.o $(UTIL_OBJS_LEVEL_2) $(BLAS_LIB) -o driver_level_2_BLAS.x $(LDFLAGS)

# Benchmark BLAS 'copy' operation.
benchmark_level_1_BLAS_copy: driver_level_1_BLAS.x
	@echo "Benchmarking Level-1 'copy' operation."
	@./driver_level_1_BLAS.x copy $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_copy_output.csv
	@python3 plotter.py "level_1_BLAS_copy"

# Benchmark BLAS 'scal' operation.
benchmark_level_1_BLAS_scal: driver_level_1_BLAS.x
	@echo "Benchmarking Level-1 'scal' operation."
	@./driver_level_1_BLAS.x scal $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_scal_output.csv
	@python3 plotter.py "level_1_BLAS_scal"

# Benchmark BLAS 'axpy' operation.
benchmark_level_1_BLAS_axpy: driver_level_1_BLAS.x
	@echo "Benchmarking Level-1 'axpy' operation."
	@./driver_level_1_BLAS.x axpy $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_axpy_output.csv
	@python3 plotter.py "level_1_BLAS_axpy"

# Benchmark BLAS 'dot' operation.
benchmark_level_1_BLAS_dot: driver_level_1_BLAS.x
	@echo "Benchmarking Level-1 'dot' operation."
	@./driver_level_1_BLAS.x dot $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_dot_output.csv
	@python3 plotter.py "level_1_BLAS_dot"

# Benchmark BLAS 'nrm2' operation.
benchmark_level_1_BLAS_nrm2: driver_level_1_BLAS.x
	@echo "Benchmarking Level-1 'nrm2' operation."
	@./driver_level_1_BLAS.x nrm2 $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_nrm2_output.csv
	@python3 plotter.py "level_1_BLAS_nrm2"

# Benchmark BLAS 'asum' operation.
benchmark_level_1_BLAS_asum: driver_level_1_BLAS.x
	@echo "Benchmarking Level-1 'asum' operation."
	@./driver_level_1_BLAS.x asum $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_asum_output.csv
	@python3 plotter.py "level_1_BLAS_asum"

# benchmark all the Level-1 BLAS operations.
benchmark_level_1_BLAS: benchmark_level_1_BLAS_copy benchmark_level_1_BLAS_scal benchmark_level_1_BLAS_axpy benchmark_level_1_BLAS_dot benchmark_level_1_BLAS_nrm2 benchmark_level_1_BLAS_asum

# Benchmark BLAS 'gemv' operation.
benchmark_level_2_BLAS_gemv: driver_level_2_BLAS.x
	@echo "Benchmarking Level-2 'gemv' operation."
	@./driver_level_2_BLAS.x gemv $(NREPEATS) $(NFIRST) $(NLAST) $(NINC) > out/level_2_BLAS_gemv_output.csv
	@python3 plotter.py "level_2_BLAS_gemv"

# Benchmark BLAS 'symv' operation.
benchmark_level_2_BLAS_symv: driver_level_2_BLAS.x
	@echo "Benchmarking Level-2 'symv' operation."
	@./driver_level_2_BLAS.x symv $(NREPEATS) $(NFIRST) $(NLAST) $(NINC) > out/level_2_BLAS_symv_output.csv
	@python3 plotter.py "level_2_BLAS_symv"

# Benchmark BLAS 'trsv' operation.
benchmark_level_2_BLAS_trsv: driver_level_2_BLAS.x
	@echo "Benchmarking Level-2 'trsv' operation."
	@./driver_level_2_BLAS.x trsv $(NREPEATS) $(NFIRST) $(NLAST) $(NINC) > out/level_2_BLAS_trsv_output.csv
	@python3 plotter.py "level_2_BLAS_trsv"

# Benchmark BLAS 'ger' operation.
benchmark_level_2_BLAS_ger: driver_level_2_BLAS.x
	@echo "Benchmarking Level-2 'ger' operation."
	@./driver_level_2_BLAS.x ger $(NREPEATS) $(NFIRST) $(NLAST) $(NINC) > out/level_2_BLAS_ger_output.csv
	@python3 plotter.py "level_2_BLAS_ger"

# Benchmark BLAS 'syr' operation.
benchmark_level_2_BLAS_syr: driver_level_2_BLAS.x
	@echo "Benchmarking Level-2 'syr' operation."
	@./driver_level_2_BLAS.x syr $(NREPEATS) $(NFIRST) $(NLAST) $(NINC) > out/level_2_BLAS_syr_output.csv
	@python3 plotter.py "level_2_BLAS_syr"

# Benchmark BLAS 'syr2' operation.
benchmark_level_2_BLAS_syr2: driver_level_2_BLAS.x
	@echo "Benchmarking Level-2 'syr2' operation."
	@./driver_level_2_BLAS.x syr2 $(NREPEATS) $(NFIRST) $(NLAST) $(NINC) > out/level_2_BLAS_syr2_output.csv
	@python3 plotter.py "level_2_BLAS_syr2"

# benchmark all the Level-2 BLAS operations.
benchmark_level_2_BLAS: benchmark_level_2_BLAS_gemv benchmark_level_2_BLAS_symv benchmark_level_2_BLAS_trsv benchmark_level_2_BLAS_ger benchmark_level_2_BLAS_syr benchmark_level_2_BLAS_syr2


# clean everything.
clean:
	rm -f *.o *~ core *.x *.pdf
