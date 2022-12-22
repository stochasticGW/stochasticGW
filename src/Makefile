### Compiler options
FCMPI = mpifort
MPIFLG = -DMPI -O3
PREPROCFLG = -Mpreprocess

### Preprocessor options
USEGPU = yes
USELIBXC = yes

### Paths to libraries

# FFTW (required)
FFT_PATH = -L/usr/bin/fftw-3.3.9/lib -lfftw3

# LIBXC (optional, required to use functionals other than PW-LDA)
LIBXC_PATH = -I/home/filip/Documents/LBL/libxc/nv-5.2.2/include -L/home/filip/Documents/LBL/libxc/nv-5.2.2/lib -lxcf03 -lxc
#LIBXC_LIB = -I/usr/bin/LIBXC/include -L/usr/bin/LIBXC/lib -lxcf03 -lxc

# GPU libraries (required for GPU execution)
CUFLG = -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/cuda/lib64 -lcudart
CULIB = -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/math_libs/lib64 -lcufft -lcublas

LIBS = $(FFT_PATH)

# Process flags
LIBS = $(FFT_PATH)

USEGPU := $(strip ${USEGPU})
GPU_FILES=
ifeq (${USEGPU},yes)
   LIBS += $(CUFLG) $(CULIB)
   FCMPI += -cuda -acc=gpu -fast -Minfo=accel -DGPU_ENABLED=1
   GPU_FILES += gpu/*.f90
else
   FCMPI += -DGPU_ENABLED=0
endif

USELIBXC := $(strip ${USELIBXC})
XC_FILES=
ifeq (${USELIBXC},yes)
    LIBS += $(LIBXC_PATH)
    FCMPI += -DLIBXC_ENABLED=1
    XC_FILES += XCI/*.f90
else
    FCMPI += -DLIBXC_ENABLED=0
endif

FCMPI += ${PREPROCFLG}

# compile
all :  clean sgw clnomod 

sgw: simmpi gwm atoms kbtop ppm kbmod
	$(FCMPI) $(MPIFLG) -o sgw.x kb/*f90 libgw/*f libgw/*f90 $(GPU_FILES) *f90 $(XC_FILES) $(LIBS)

simmpi:
	$(FCMPI) $(MPIFLG) -c ./libgw/0_library_mpi_module.f90 ./libgw/0_simple_mpi.f90

gwm:
	$(FCMPI) $(MPIFLG) -c 1_gwm.f90

atoms:
	$(FCMPI) $(MPIFLG) -c 2_atoms.f90

kbtop:
	$(FCMPI) $(MPIFLG) -c 2_kb_top_module.f90

ppm:
	$(FCMPI) $(MPIFLG) -c 3_ppm.f90

kbmod:
	$(FCMPI) $(MPIFLG) -c kb/1_kb_mod.f90

clnomod:
	rm -rf *o *mod kb/*o libgw/*o 

clean:
	rm -rf *o *mod kb/*o libgw/*o sgw.x 

