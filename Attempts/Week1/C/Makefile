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
	./driver_IJP.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_IJP.log
	cat data/output_IJP.log

###############################     END     ####################################
############################### Gemm -> IPJ ####################################

OBJS_IPJ := driver.o Gemm_IPJ.o

driver_IPJ.x: $(OBJS_IPJ) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IPJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IPJ.x $(LDFLAGS)

IPJ: driver_IPJ.x
	./driver_IPJ.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_IPJ.log
	cat data/output_IPJ.log

###############################     END     ####################################
############################### Gemm -> JIP ####################################

OBJS_JIP := driver.o Gemm_JIP.o

driver_JIP.x: $(OBJS_JIP) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JIP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JIP.x $(LDFLAGS)

JIP: driver_JIP.x
	./driver_JIP.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_JIP.log
	cat data/output_JIP.log

###############################     END     ####################################
############################### Gemm -> JPI ####################################

OBJS_JPI := driver.o Gemm_JPI.o

driver_JPI.x: $(OBJS_JPI) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JPI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JPI.x $(LDFLAGS)

JPI: driver_JPI.x
	./driver_JPI.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_JPI.log
	cat data/output_JPI.log

###############################     END     ####################################
############################### Gemm -> PIJ ####################################

OBJS_PIJ := driver.o Gemm_PIJ.o

driver_PIJ.x: $(OBJS_PIJ) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PIJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PIJ.x $(LDFLAGS)

PIJ: driver_PIJ.x
	./driver_PIJ.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_PIJ.log
	cat data/output_PIJ.log

###############################     END     ####################################
############################### Gemm -> PJI ####################################

OBJS_PJI := driver.o Gemm_PJI.o

driver_PJI.x: $(OBJS_PJI) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PJI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJI.x $(LDFLAGS)

PJI: driver_PJI.x
	./driver_PJI.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_PJI.log
	cat data/output_PJI.log

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
############################### Gemv -> I_Dots #################################

OBJS_I_Dots := driver_gemv.o Dots.o Gemv_I_Dots.o

driver_I_Dots.x: $(OBJS_I_Dots) $(UTIL_OBJS)
	$(LINKER) $(OBJS_I_Dots) $(UTIL_OBJS) $(BLAS_LIB) -o driver_I_Dots.x $(LDFLAGS)

I_Dots: driver_I_Dots.x
	./driver_I_Dots.x $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_I_Dots.log
	cat data/output_I_Dots.log

###############################     END     ####################################
############################### Gemv -> J_Axpy #################################

OBJS_J_Axpy := driver_gemv.o Axpy.o Gemv_J_Axpy.o

driver_J_Axpy.x: $(OBJS_J_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_J_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_Axpy.x $(LDFLAGS)

J_Axpy: driver_J_Axpy.x
	./driver_J_Axpy.x $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_J_Axpy.log
	cat data/output_J_Axpy.log

###############################     END     ####################################
############################### J_Gemv_J_Axpy ##################################

OBJS_J_Gemv_J_Axpy := driver.o Gemm_J_Gemv.o Gemv_J_Axpy.o Axpy.o

driver_J_Gemv_J_Axpy.x: $(OBJS_J_Gemv_J_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_J_Gemv_J_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_Gemv_J_Axpy.x $(LDFLAGS)

J_Gemv_J_Axpy: driver_J_Gemv_J_Axpy.x
	./driver_J_Gemv_J_Axpy.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_J_Gemv_J_Axpy.log
	cat data/output_J_Gemv_J_Axpy.log

###############################     END     ####################################
############################### J_Gemv_I_Dots ##################################

OBJS_J_Gemv_I_Dots := driver.o Gemm_J_Gemv.o Gemv_I_Dots.o Dots.o

driver_J_Gemv_I_Dots.x: $(OBJS_J_Gemv_I_Dots) $(UTIL_OBJS)
	$(LINKER) $(OBJS_J_Gemv_I_Dots) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_Gemv_I_Dots.x $(LDFLAGS)

J_Gemv_I_Dots: driver_J_Gemv_I_Dots.x
	./driver_J_Gemv_I_Dots.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_J_Gemv_I_Dots.log
	cat data/output_J_Gemv_I_Dots.log

###############################     END     ####################################
###############################   Outer_J   ####################################

OUTER_J:
	make J_Gemv_I_Dots
	make J_Gemv_J_Axpy
	make JIP
	make JPI
	make clean

###############################     END     ####################################
###############################  Ger_J_Axpy ####################################

OBJS_Ger_J_Axpy := driver_ger.o Ger_J_Axpy.o Axpy.o

driver_Ger_J_Axpy.x: $(OBJS_Ger_J_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_Ger_J_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Ger_J_Axpy.x $(LDFLAGS)

Ger_J_Axpy: driver_Ger_J_Axpy.x
	./driver_Ger_J_Axpy.x $(NFIRST) $(NLAST) $(NINC) > data/output_Ger_J_Axpy.log
	cat data/output_Ger_J_Axpy.log

###############################     END     ####################################
###############################  Ger_I_Axpy ####################################

OBJS_Ger_I_Axpy := driver_ger.o Ger_I_Axpy.o Axpy.o

driver_Ger_I_Axpy.x: $(OBJS_Ger_I_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_Ger_I_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Ger_I_Axpy.x $(LDFLAGS)

Ger_I_Axpy: driver_Ger_I_Axpy.x
	./driver_Ger_I_Axpy.x $(NFIRST) $(NLAST) $(NINC) > data/output_Ger_I_Axpy.log
	cat data/output_Ger_I_Axpy.log

###############################     END     ####################################
############################### P_Ger_I_Axpy ###################################

OBJS_P_Ger_I_Axpy := driver.o Gemm_P_Ger.o Ger_I_Axpy.o Axpy.o

driver_P_Ger_I_Axpy.x: $(OBJS_P_Ger_I_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_P_Ger_I_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_Ger_I_Axpy.x $(LDFLAGS)

P_Ger_I_Axpy: driver_P_Ger_I_Axpy.x
	./driver_P_Ger_I_Axpy.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_P_Ger_I_Axpy.log
	cat data/output_P_Ger_I_Axpy.log

###############################     END     ####################################
############################### P_Ger_J_Axpy ###################################

OBJS_P_Ger_J_Axpy := driver.o Gemm_P_Ger.o Ger_J_Axpy.o Axpy.o

driver_P_Ger_J_Axpy.x: $(OBJS_P_Ger_J_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_P_Ger_J_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_Ger_J_Axpy.x $(LDFLAGS)

P_Ger_J_Axpy: driver_P_Ger_J_Axpy.x
	./driver_P_Ger_J_Axpy.x $(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC) > data/output_P_Ger_J_Axpy.log
	cat data/output_P_Ger_J_Axpy.log

###############################     END     ####################################
###############################   Outer_P   ####################################

OUTER_P:
	make P_Ger_I_Axpy
	make P_Ger_J_Axpy
	make PIJ
	make PJI
	make clean

###############################     END     ####################################


all:
	make ALL_ORDERINGS
	make OUTER_J
	make OUTER_P
	make I_Dots
	make J_Axpy
	make Ger_J_Axpy
	make Ger_I_Axpy
	make clean
	echo -e "\n\nall tests pass\n\n"

# clean everything
clean:
	rm -f *.o *~ core *.x *.pdf
