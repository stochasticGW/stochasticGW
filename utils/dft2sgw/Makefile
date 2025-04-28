FC=gfortran
MPIFLG=
FFT_PATH = -L/usr/bin/fftw-3.3.9/lib -lfftw3
LIBS = $(FFT_PATH)

all : clean dft2sgw

dft2sgw:
	$(FC) $(MPIFLG) -o dft2sgw.x utils.f90 procoptions.f90 fftinterpolate.f90 readinfile.f90 readdftout.f90 readcube.f90 writesgw.f90 dft2sgw.f90 $(LIBS)

clean:
	rm -rf *.o *.mod dft2sgw.x
