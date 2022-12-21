      SUBROUTINE r_inv_mtrx_lrg(A,A_inv,N, LDA)
      IMPLICIT none
      INTEGER IPVT(20000), N, Job, LDA, INFO
      REAL*8 A(LDA,LDA), A_inv(LDA, LDA), WORK(20000), DET(2) 

      A_inv = A

      JOB = 1
C ! ONLY INVERESE
      IF(N.GT.20000) THEN
          WRITE(6,*)' N = ',N,' DECREASE IT OR INCREASE DIM. '
          STOP
      ENDIF

      CALL RZGEFA(A_inv,LDA,N,IPVT,INFO)
      if (info /= 0) then
         write(6,*)' problem in r_inv_mtrx ';stop
      endif
      CALL RRZGEDI(A_inv,LDA,N,IPVT,DET,WORK,JOB)

      END

      SUBROUTINE r_inv_mtrx(A,A_inv,N)
      IMPLICIT none
      INTEGER IPVT(20000), N, Job, LDA, INFO
      REAL*8 A(N,N), A_inv(N, N), WORK(20000), DET(2) 

      LDA = N
      A_inv = A

      JOB = 1
C ! ONLY INVERESE
      IF(N.GT.20000) THEN
          WRITE(6,*)' N = ',N,' DECREASE IT OR INCREASE DIM. '
          STOP
      ENDIF

      CALL RZGEFA(A_inv,LDA,N,IPVT,INFO)
      if (info /= 0) then
         write(6,*)' problem in r_inv_mtrx ';stop
      endif
      CALL RRZGEDI(A_inv,LDA,N,IPVT,DET,WORK,JOB)

      END


      SUBROUTINE r_inv_det_mtrx_good(A,A_inv,DET,N)
      IMPLICIT none
      INTEGER N
      REAL*8 A(N,N), A_inv(N, N),  DETA(2) ,DET
      call r_inv_det_mtrx_tricky(A,A_inv,DETA,N)
      DET = DETA(1) * 10.0**DETA(2)
      end

      SUBROUTINE r_inv_det_mtrx_tricky(A,A_inv,DET,N)
      IMPLICIT none
      INTEGER IPVT(20000), N, Job, LDA, INFO
      REAL*8 A(N,N), A_inv(N, N), WORK(20000), DET(2) 


      LDA = N
      A_inv = A

      JOB = 11
      IF(N.GT.20000) THEN
          WRITE(6,*)' N = ',N,' DECREASE IT OR INCREASE DIM. '
          STOP
      ENDIF

      CALL RZGEFA(A_inv,LDA,N,IPVT,INFO)
      if (info /= 0) then
         write(6,*)' problem in r_inv_mtrx ';det=0;return
      endif
      CALL RRZGEDI(A_inv,LDA,N,IPVT,DET,WORK,JOB)

      END

     

      SUBROUTINE RRZGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER LDA,N,IPVT(*),JOB
      REAL*8 A(LDA,*),DET(2),WORK(*)
C
C     RRZGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY ZGECO OR RZGEFA.
C
C     ON ENTRY
C
C        A       REAL*8(LDA, N)
C                THE OUTPUT FROM ZGECO OR RZGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM ZGECO OR RZGEFA.
C
C        WORK    REAL*8(N)
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
C        DET     REAL*8(2)
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
C        AND IF ZGECO HAS SET RCOND .GT. 0.0 OR RZGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS RZAXPY,RZSCAL,RZSWAP
C     FORTRAN DABS, CMPLX,MOD
C
C     INTERNAL VARIABLES
C
      REAL*8 T
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
      REAL*8 ZDUM
      DOUBLE PRECISION ABS1
      REAL*8 ZDUMR,ZDUMI
      real*8 temp(n)
CC      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
      ABS1(ZDUM) = DABS(ZDUM) 
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
            IF (ABS1(DET(1)) .EQ. 0.0D0) GO TO 60
   10       IF (ABS1(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) =  CMPLX(TEN,0.0D0)*DET(1)
               DET(2) = DET(2) - (1.0D0,0.0D0)
            GO TO 10
   20       CONTINUE
   30       IF (ABS1(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/ CMPLX(TEN,0.0D0)
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
            CALL RZSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = (0.0D0,0.0D0)
               !!!! CALL RZAXPY(K,T,A(1,K),1,A(1,J),1)
               a(1:k,j) = a(1:k,j)+ t*a(1:k,k)
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
               !!! CALL RZAXPY(N,T,A(1,J),1,A(1,K),1)
               a(1:n,k) = a(1:n,k)+ t*a(1:n,j)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) then
                !!!! CALL RZSWAP(N,A(1,K),1,A(1,L),1)
               temp = a(:,k)
               a(:,k) = a(:,l)
               a(:,l) = temp
            endif
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END


      SUBROUTINE RZGEFA(A,LDA,N,IPVT,INFO)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER LDA,N,IPVT(*),INFO
      REAL*8 A(LDA,*)
      !real*8 dum(lda,*)
C
C     RZGEFA FACTORS A REAL*8 MATRIX BY GAUSSIAN ELIMINATION.
C
C     RZGEFA IS USUALLY CALLED BY ZGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR ZGECO) = (1 + 9/N)*(TIME FOR RZGEFA) .
C
C     ON ENTRY
C
C        A       REAL*8(LDA, N)
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
C                     INDICATE THAT ZGESL OR RRZGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN ZGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS RZAXPY,RZSCAL,JZAMAX
C     FORTRAN DABS
C
C     INTERNAL VARIABLES
C
      REAL*8 T
      INTEGER JZAMAX,J,K,KP1,L,NM1
C
      REAL*8 ZDUM
      DOUBLE PRECISION ABS1
      REAL*8 ZDUMR,ZDUMI
      real*8 term

CC      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
      ABS1(ZDUM) = DABS(ZDUM) 
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
         L = JZAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (ABS1(A(L,K)) .EQ. 0.0D0) GO TO 40
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
            CALL RZSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
!!!!!         CALL RZAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
               a(k+1:n,j) = a(k+1:n,j)+ t*a(k+1:n,k)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (ABS1(A(N,N)) .EQ. 0.0D0) INFO = N
      RETURN
      END


      SUBROUTINE RZAXPY(N,ZA,ZX,INCX,ZY,INCY)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     JACK DONGARRA, 3/11/78.
C
      REAL*8 ZX(*),ZY(*),ZA
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
        ZY(IY) = ZY(IY) + ZA*ZX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
      
   20 DO 30 I = 1,N
       ZY(I) = ZY(I) + ZA*ZX(I)
   30 CONTINUE
      RETURN
      END

      INTEGER FUNCTION JZAMAX(N,ZX,INCX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, 1/15/85.
C
      REAL*8 ZX(*)
      DOUBLE PRECISION SMAX
      INTEGER I,INCX,IX,N
C
      JZAMAX = 0
      IF(N.LT.1)RETURN
      JZAMAX = 1
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
         JZAMAX = I
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
         JZAMAX = I
         SMAX = ABS(ZX(I))
   30 CONTINUE
      RETURN
      END

      SUBROUTINE  RZSCAL(N,ZA,ZX,INCX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C    SCALES A VECTOR BY A CONSTANT.
C    JACK DONGARRA, 3/11/78.
C
      REAL*8 ZA,ZX(*)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1
C
      IX = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        ZX(IX) = ZA*ZX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        ZX(I) = ZA*ZX(I)
   30 CONTINUE
      RETURN
      END

      SUBROUTINE  RZSWAP (N,ZX,INCX,ZY,INCY)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     INTERCHANGES TWO VECTORS.
C     JACK DONGARRA, 3/11/78.
C
      REAL*8 ZX(*),ZY(*)
      REAL*8 ZTEMP
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        ZTEMP = ZX(IX)
        ZX(IX) = ZY(IY)
        ZY(IY) = ZTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
   20 DO 30 I = 1,N
        ZTEMP = ZX(I)
        ZX(I) = ZY(I)
        ZY(I) = ZTEMP
   30 CONTINUE
      RETURN
      END

