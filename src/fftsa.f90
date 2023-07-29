!***********************************************************************
!
!                      SUBROUTINE FFTSA
!
!     THIS IS A SIMPLE FFT ROUTINE FOR COMPLEX VECTORS OF 
!     LENGTH N=2**M.
!
!     ARGUEMENTS:  
!        N         LENGTH OF THE VECTOR
!        A         COMPLEX VECTOR OF LENGTH N. ON INPUT, A(k)=f[r(k)]
!                  ON OUTPUT, A(k)=sum(j=1,N) f[r(j)] * exp[-i*r(j)*p(k)] 
!        M         N=2**M
!
!***********************************************************************
      SUBROUTINE FFTSA(N,A,M)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 A(*),U,W,T,ZI
      PIPI = DACOS(-1.D0)
      N11=2**M
      IF(N.NE.N11) THEN
         WRITE(6,*) ' *** FATAL ERROR IN ROUTINE FFTSA '
         WRITE(6,*) ' *** N NOT EQUAL 2**M: (N,M) - ', N, M
         STOP
      ENDIF
      NM1=N11-1
      NV2=N11/2
      J=1
      DO 7 I=1,NM1
      IF (I.GE.J) GOTO 5
      T=A(J)
      A(J) = A(I)
      A(I) = T
5     K=NV2
6     IF (K.GE.J) GOTO 7
      J =J- K
      K=K/2
      GOTO 6
7     J=J+K
      DO 20 L=1, M
      LE = 2** L
      LE1 = LE/2
      U=(1.D0,0.D0)
      ANG=PIPI / LE1
      ZI=(0.D0,1.D0)
      W= DCOS(ANG)+ZI* DSIN(ANG)
      DO 20 J=1, LE1
      DO 10 I=J,N11, LE
      IP=I + LE1
      T= A(IP)*CONJG(U)
      A(IP) = A(I) - T
10    A(I)  = A(I) + T
20    U=U*W
      RETURN
      END


