FC=gfortran
MPIFLG=
#FC=mpifort 
#MPIFLG= -DMPI -O3

all : clean qe2sgw qe2sgwstepwise

qe2sgw:
	${FC} ${MPIFLG} -o qe2sgw.x qe2sgw.f90

qe2sgwstepwise:
	${FC} ${MPIFLG} -o qe2sgw_stepwise.x qe2sgw_stepwise.f90

clean:
	rm -rf *.o *.mod qe2sgw.x qe2sgw_stepwise.x
