FC=gfortran
MPIFLG=

all : clean dft2sgw

dft2sgw:
	${FC} ${MPIFLG} -o dft2sgw.x utils.f90 procoptions.f90 readinfile.f90 readdftout.f90 readcube.f90 writesgw.f90 dft2sgw.f90

clean:
	rm -rf *.o *.mod dft2sgw.x
