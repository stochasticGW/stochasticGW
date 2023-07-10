!
! From Numerical Recipes (www.nr.com) with tiny changes
!  
!  Numerical Recipes routine for preparing spline coeff.  Given X, Y arrays, N
!  and two numbers, YP1 and YP2 (which should be set at 10^32 exactly
!  for natural splines) 
!  the routine reutrns Y2, used for spline interpolations.
!
!  The routine was modified to use allocation, so NMAX is not used for U
!  any more (f77-->f90)


      SUBROUTINE SPLINE_dn(X,Y,N,Y2)
      implicit none
      integer  i, N,k
      real*8   yp1, ypn
      real*8  X(N),Y(N),Y2(N)
      real*8, allocatable ::  u(:)
      real*8  sig, p, qn, un

      allocate(u(n),stat=i)

      if(i/=0) then
         write(6,*)' problem allocating u in spline; code = ',i
         stop
      endif

      yp1 = 10.d0**32
      ypn = 10.d0**32
      
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.d0
        U(1)=0.d0
      ELSE
        Y2(1)=-0.5d0
        U(1)=(3.d0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.d0
        Y2(I)=(SIG-1.d0)/P
        U(I)=(6.d0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.d0
        UN=0.d0
      ELSE
        QN=0.5d0
        UN=(3.d0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.d0)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      deallocate(u)
      RETURN
      END


