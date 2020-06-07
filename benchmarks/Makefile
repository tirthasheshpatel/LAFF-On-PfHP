# HOME      := /Users/rvdg

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
NLAST_SMALL:= 500#       largest size to be timed for slow implementations.
NLAST      := 1500#      largest size to be timed for fast implementations.
NINC       := 48#        increment between sizes.

LDFLAGS    := -lpthread -m64 -lm -fopenmp

UTIL_OBJS  := FLA_Clock.o RandomMatrix.o

# compile the driver.
driver_level_1_BLAS.x: level_1_BLAS.o $(UTIL_OBJS)
	$(LINKER) level_1_BLAS.o $(UTIL_OBJS) $(BLAS_LIB) -o driver_level_1_BLAS.x $(LDFLAGS)

# Benchmark BLAS 'copy' operation.
benchmark_level_1_BLAS_copy: driver_level_1_BLAS.x
	./driver_level_1_BLAS.x copy $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_copy_output.csv
	python3 plotter.py "level_1_BLAS_copy"
	make clean

# Benchmark BLAS 'scal' operation.
benchmark_level_1_BLAS_scal: driver_level_1_BLAS.x
	./driver_level_1_BLAS.x scal $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_scal_output.csv
	python3 plotter.py "level_1_BLAS_scal"
	make clean

# Benchmark BLAS 'axpy' operation.
benchmark_level_1_BLAS_axpy: driver_level_1_BLAS.x
	./driver_level_1_BLAS.x axpy $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_axpy_output.csv
	python3 plotter.py "level_1_BLAS_axpy"
	make clean

# Benchmark BLAS 'dot' operation.
benchmark_level_1_BLAS_dot: driver_level_1_BLAS.x
	./driver_level_1_BLAS.x dot $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_dot_output.csv
	python3 plotter.py "level_1_BLAS_dot"
	make clean

# Benchmark BLAS 'nrm2' operation.
benchmark_level_1_BLAS_nrm2: driver_level_1_BLAS.x
	./driver_level_1_BLAS.x nrm2 $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_nrm2_output.csv
	python3 plotter.py "level_1_BLAS_nrm2"
	make clean

# Benchmark BLAS 'asum' operation.
benchmark_level_1_BLAS_asum: driver_level_1_BLAS.x
	./driver_level_1_BLAS.x asum $(NREPEATS) 1000 5000 50 > out/level_1_BLAS_asum_output.csv
	python3 plotter.py "level_1_BLAS_asum"
	make clean

# benchmark all the Level-1 BLAS operations.
benchmark_level_1_BLAS:
	make benchmark_level_1_BLAS_copy
	make benchmark_level_1_BLAS_scal
	make benchmark_level_1_BLAS_axpy
	make benchmark_level_1_BLAS_dot
	make benchmark_level_1_BLAS_nrm2
	make benchmark_level_1_BLAS_asum

# clean everything.
clean:
	rm -f *.o *~ core *.x *.pdf
