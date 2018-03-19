CC = gcc
CC_DBG=gcc -O0 -g2
MCC = mcc --ompss
SMPCC = smpcc --ompss
SMPCC_DBG = smpcc --ompss -g2 -O0 -k --debug
#INSTRUMENTATION = --instrument
BSC_CC = $(SMPCC)


MKL_LINK = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread
MKL_INC = -I${MKLROOT}/include

LAPACK_LINK = /apps/CBLAS/lib/cblas_LINUX.a -L/apps/LAPACK/3.5.0/GCC/lib -lblas
LAPACK_INC = -I/apps/CBLAS/include


### TODO INTEL_MKL | LAPACK
MATH = INTEL_MKL
ifeq ($(MATH), INTEL_MKL)
	MATH_INC = $(MKL_INC)
	MATH_LINK = $(MKL_LINK)
endif

ifeq ($(MATH), LAPACK)
	MATH_INC = $(LAPACK_INC)
	MATH_LINK = $(LAPACK_LINK)
endif

CFLAGS = -c -Wall -D$(MATH) $(MATH_INC) $(INSTRUMENTATION)
LDFLAGS = $(INSTRUMENTATION) $(MATH_LINK) -lm

ALL = dcg

ALLSRC = hb_io.c csparse.c cg_setup.c cg_main.c cg_config.c cg_aux.c

ALLOBJ = $(ALLSRC:.c=.o)

$(ALL) : $(ALLOBJ)
	$(BSC_CC) -o $@ $^ $(LDFLAGS)

%.o : %.c
	$(BSC_CC) $(CFLAGS) -o $@ $<

#For the compatibility
.c.o:
	$(BSC_CC) $(CFLAGS) -o $@ $<

.PHONY: clean

clean:
	rm -f $(ALLOBJ) $(ALL) smpcc_*.c mcc_*.c
