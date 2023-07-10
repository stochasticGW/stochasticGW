
      SUBROUTINE c_det_only(Ain, N, c_det)
      implicit none
      integer info, job,lda, N, ipvt(20000)
      COMPLEX*16 A(N,N), Ain(N,N),WORK(20000), DET(2) , c_det

      A = Ain
      lda = N

      JOB = 10                ! only
      IF(N.GT.20000) THEN
          WRITE(6,*)' N = ',N,' DECREASE IT OR INCREASE DIM. '
          STOP
      ENDIF

      CALL ZGEFA(A,LDA,N,IPVT,INFO)
      if (info /= 0) then
         if(n <= 10) then
            write(86,*)' A ', A
            write(86,*)' A ',A
         endif
         write(6,*)' problem in c_mtrx_with_det ',info;stop
      endif
      CALL ZGEDI(A,LDA,N,IPVT,DET,WORK,JOB)

      if(abs(det(2)-nint(real(det(2))))>1.d-13)stop   ! checks that det(2) int
      if(abs(det(2))>70) then
         write(6,*)' possible problem. c_det(2) : ',det
         stop
      endif
      c_det = det(1)* 10.d0**(nint(real(det(2))))
      END


      SUBROUTINE c_inv_mtrx_with_det(A, A_inv, N, c_det)
      implicit none
      integer info, job,lda, N, ipvt(20000)
      COMPLEX*16 A(N,N), A_inv(N, N), WORK(20000), DET(2) , c_det

      lda = N
      A_inv = A

      JOB = 11                ! with det
      IF(N.GT.20000) THEN
          WRITE(6,*)' N = ',N,' DECREASE IT OR INCREASE DIM. '
          STOP
      ENDIF

      CALL ZGEFA(A_inv,LDA,N,IPVT,INFO)
      if (info /= 0) then
         if(n <= 10) then
            write(86,*)' A ', A
            write(86,*)' A_inv ',A_inv
         endif
         write(6,*)' problem in c_inv_mtrx_with_det ',info;stop
      endif
      CALL ZGEDI(A_inv,LDA,N,IPVT,DET,WORK,JOB)

      if(abs(det(2)-nint(real(det(2))))>1.d-13)stop   ! checks that det(2) int
      c_det = det(1)* 10.d0**(nint(real(det(2))))
      END


      SUBROUTINE c_inv_mtrx(A, A_inv, N)
      implicit none
      integer info, job,lda, N, ipvt(20000)
      COMPLEX*16 A(N,N), A_inv(N, N), WORK(20000), DET(2) 

      lda = N
      A_inv = A

      JOB = 1                ! ONLY INVERESE
      IF(N.GT.20000) THEN
          WRITE(6,*)' N = ',N,' DECREASE IT OR INCREASE DIM. '
          STOP
      ENDIF

      CALL ZGEFA(A_inv,LDA,N,IPVT,INFO)
      if (info /= 0) then
         if(n <= 10) then
            write(86,*)' A '
            write(86,*)' A_inv ',A_inv
         endif

         write(6,*)' problem in c_inv_mtrx ',info;stop
      endif
      CALL ZGEDI(A_inv,LDA,N,IPVT,DET,WORK,JOB)

      END

     

      SUBROUTINE ZGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      COMPLEX*16 A(LDA,N),DET(2),WORK(1)
C
C     ZGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY ZGECO OR ZGEFA.
C
C     ON ENTRY
C
C        A       COMPLEX*16(LDA, N)
C                THE OUTPUT FROM ZGECO OR ZGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM ZGECO OR ZGEFA.
C
C        WORK    COMPLEX*16(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     COMPLEX*16(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. ABS(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF ZGECO HAS SET RCOND .GT. 0.0 OR ZGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS ZAXPY,ZSCAL,ZSWAP
C     FORTRAN DABS,DCMPLX,MOD
C
C     INTERNAL VARIABLES
C
      COMPLEX*16 T
      REAL*8 TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
      COMPLEX*16 ZDUM , ztemp(n)
      COMPLEX*16 ZDUMR,ZDUMI

C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = (1.0D0,0.0D0)
         DET(2) = (0.0D0,0.0D0)
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (ABS(DET(1)) .EQ. 0.0D0) GO TO 60
   10       IF (ABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = CMPLX(TEN,0.0D0)*DET(1)
               DET(2) = DET(2) - (1.0D0,0.0D0)
            GO TO 10
   20       CONTINUE
   30       IF (ABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/CMPLX(TEN,0.0D0)
               DET(2) = DET(2) + (1.0D0,0.0D0)
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = (1.0D0,0.0D0)/A(K,K)
            T = -A(K,K)
            CALL ZSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = (0.0D0,0.0D0)
!!               CALL ZAXPY(K,T,A(1,K),1,A(1,J),1)
               a(1:k,J) =a(1:k,j)+ a(1:k,k)*t 
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = (0.0D0,0.0D0)
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
!!               CALL ZAXPY(N,T,A(1,J),1,A(1,K),1)
               a(1:N,k) =a(1:N,k)+ a(1:N,j)*t 
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) then
!!!!!!!               cALL ZSWAP(N,A(1,K),1,A(1,L),1)   ! note, erased
               ztemp(1:n) = a(1:n,K)
               a(1:n,K) = a(1:n,L)
               a(1:n,L) = ztemp(1:n)
            endif
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END


      SUBROUTINE ZGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(1),INFO
      COMPLEX*16 A(LDA,1)
C
C     ZGEFA FACTORS A COMPLEX*16 MATRIX BY GAUSSIAN ELIMINATION.
C
C     ZGEFA IS USUALLY CALLED BY ZGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR ZGECO) = (1 + 9/N)*(TIME FOR ZGEFA) .
C
C     ON ENTRY
C
C        A       COMPLEX*16(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT ZGESL OR ZGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN ZGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS ZAXPY,ZSCAL,IZAMAX
C     FORTRAN DABS
C
C     INTERNAL VARIABLES
C
      COMPLEX*16 T
      INTEGER IZAMAX,J,K,KP1,L,NM1
C
      COMPLEX*16 ZDUM
      REAL*8 ABS
      COMPLEX*16 ZDUMR,ZDUMI
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IZAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (ABS(A(L,K)) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -(1.0D0,0.0D0)/A(K,K)
            CALL ZSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE

! new,avoids compiler problems
!!               CALL ZAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
                 a(k+1:n,J) =a(k+1:n,j)+ a(k+1:n,k)*t 

   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (ABS(A(N,N)) .EQ. 0.0D0) INFO = N
      RETURN
      END


      SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     JACK DONGARRA, 3/11/78.
C
      COMPLEX*16 ZX(1),ZY(1),ZA
      complex*16 zztop(n)
      IF(N.LE.0)RETURN
      IF (ABS(ZA) .EQ. 0.0D0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        zztop(i) = ZY(IY) + ZA*ZX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      iy=1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      do i=1,n
         iy=iy+incy
         zy(iy) = zztop(i)
      enddo
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        Zztop(I) = ZY(I) + ZA*ZX(I)
   30 CONTINUE
      iy=1
      do i=1,n
         zy(iy) = zztop(i)
         iy=iy+incy
      enddo
      RETURN
      END

      INTEGER FUNCTION IZAMAX(N,ZX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, 1/15/85.
C
      COMPLEX*16 ZX(1)
      REAL*8 SMAX
      INTEGER I,INCX,IX,N
C
      IZAMAX = 0
      IF(N.LT.1)RETURN
      IZAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      SMAX = ABS(ZX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(ABS(ZX(IX)).LE.SMAX) GO TO 5
         IZAMAX = I
         SMAX = ABS(ZX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 SMAX = ABS(ZX(1))
      DO 30 I = 2,N
         IF(ABS(ZX(I)).LE.SMAX) GO TO 30
         IZAMAX = I
         SMAX = ABS(ZX(I))
   30 CONTINUE
      RETURN
      END

