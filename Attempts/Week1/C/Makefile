BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

CC          := gcc
LINKER      := $(CC)
CFLAGS      := -O3 -I$(BLAS_INC) -m64 -mavx2 -std=c99 -march=native -fopenmp -D_POSIX_C_SOURCE=200809L
FFLAGS      := $(CFLAGS)

NREPEATS    := 3
NFIRST      := 48
NLAST_SMALL := 500
NLAST       := 1500
NINC        := 48

LDFLAGS     := -lpthread -m64 -lm -fopenmp

UTIL_OBJS   := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o


############################### Gemm -> IJP ####################################

OBJS_IJP := driver.o Gemm_IJP.o

driver_IJP.x: $(OBJS_IJP) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IJP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJP.x $(LDFLAGS)

IJP: driver_IJP.x
	./driver_IJP.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_IJP.txt
	cat data/output_IJP.txt

###############################     END     ####################################
############################### Gemm -> IPJ ####################################

OBJS_IPJ := driver.o Gemm_IPJ.o

driver_IPJ.x: $(OBJS_IPJ) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IPJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IPJ.x $(LDFLAGS)

IPJ: driver_IPJ.x
	./driver_IPJ.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_IPJ.txt
	cat data/output_IPJ.txt

###############################     END     ####################################
############################### Gemm -> JIP ####################################

OBJS_JIP := driver.o Gemm_JIP.o

driver_JIP.x: $(OBJS_JIP) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JIP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JIP.x $(LDFLAGS)

JIP: driver_JIP.x
	./driver_JIP.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_JIP.txt
	cat data/output_JIP.txt

###############################     END     ####################################
############################### Gemm -> JPI ####################################

OBJS_JPI := driver.o Gemm_JPI.o

driver_JPI.x: $(OBJS_JPI) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JPI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JPI.x $(LDFLAGS)

JPI: driver_JPI.x
	./driver_JPI.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_JPI.txt
	cat data/output_JPI.txt

###############################     END     ####################################
############################### Gemm -> PIJ ####################################

OBJS_PIJ := driver.o Gemm_PIJ.o

driver_PIJ.x: $(OBJS_PIJ) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PIJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PIJ.x $(LDFLAGS)

PIJ: driver_PIJ.x
	./driver_PIJ.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_PIJ.txt
	cat data/output_PIJ.txt

###############################     END     ####################################
############################### Gemm -> PJI ####################################

OBJS_PJI := driver.o Gemm_PJI.o

driver_PJI.x: $(OBJS_PJI) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PJI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJI.x $(LDFLAGS)

PJI: driver_PJI.x
	./driver_PJI.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_PJI.txt
	cat data/output_PJI.txt

###############################     END     ####################################
############################### All Orderings ##################################

ALL_ORDERINGS:
	make IJP
	make IPJ
	make JIP
	make JPI
	make PIJ
	make PJI
	make clean

###############################     END     ####################################

# clean everything
clean:
	rm -f *.o *~ core *.x *.pdf
