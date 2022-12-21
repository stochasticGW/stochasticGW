!                                                   
! From Numerical Recipes (www.nr.com) with tiny chnges.
!
!  Given XA, YA, Y2A (output of spline) N and X -- this routine returns Y
!  modified for real*8, implicit none, etc.  See numerical recipes , chap. 3
!

      SUBROUTINE SPLINT_dn(XA,YA,Y2A,N,X,Y)
      implicit none
      integer n, klo, khi, k
      real*8    x, y, h, a, b
      real*8    XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.d0) then
         write(6,*) 'Bad XA input.'
         stop
      endif
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.d0
      RETURN
      END
