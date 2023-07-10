      subroutine driver_diag()  ! dummy subroutine;checks rdiag
      implicit none
      integer M; parameter(M=4)
      integer i, j
      real*8 A(M,M), Vc(M,M), w(M), U(M,M), A_dum(M, M)

      do i=1, M
         do j=1, M
            A(i, j) = i*2+j*2-i**2*10-j**2*10+(i-j)**4
         enddo
      enddo

      write(6,*)' A '
      write(6,88)A

      A_dum = A
      A_dum(:,2) = A(:, 3)
      A_dum(:,3) = A(:, 2)

      A = A_dum
      A_dum(2,:) = A(3, :)
      A_dum(3,:) = A(2, :)
      A = A_dum

      write(6,*)' A '
      write(6,88)A

      call r_diag_normlz(A, Vc, w, M)

      write(6,*)' A '
      write(6,88) A
      write(6,*)' w '
      write(6,88) w
      write(6,*)' Vc '
      do i=1, M
      write(6,88) Vc(i,:) 
      enddo

 88   format(1x,4f15.6)
      
      end
      
c-----
c generic subroutine to diagonalize A(M,M) with e.values onlu
c          A is here symmetric and real. If not, apply rg or cg.
c----
      subroutine r_diag_evl_only(A, w, M)
      implicit none
      
      integer M, st, matz, ierr, i, h, j
      real*8 A(M, M), w(M)
      real*8, allocatable :: vc_dum(:,:)
      
      real*8, dimension(:), allocatable :: fv1, fv2
c----
c  first check symmetry
c---
      do i=1,M ; do j=1, M
         if(abs(A(i,j)-A(j,i))>1.d-9) then
            write(6,*)' problem -- A is not symmetrical '
            write(6,*)' i, j, A ',i, j, A(i,j), A(j,i)
            stop
         endif
      enddo;     enddo
c-----
c  diagonalize
c----

      allocate (fv1(M), fv2(M), stat= st); if(st.ne.0) stop
      matz = 0
      call rs(M, M, A, w, matz, Vc_dum, fv1, fv2, ierr)
      if(ierr /=0) then; write(6,*)' ie',ierr;stop;endif
      deallocate(fv1, fv2)
!      do i=1, M
!         write(6,*)i,w(i)
!      enddo
!      call flush(6)
!      write(6,*)' minval ',minval(w)
!      write(6,*)' maxval ',maxval(w)
!      call flush(6)
!      return
      end
c-----
c generic subroutine to diagonalize A(M,M); results fulfil A*vc = eps*A. 
c                       also: each vector in vc is orthonormalized.
c                       and: eigenvalues and eigenvectors are ordered.
c          A is here symmetric and real. If not, apply rg or cg.
c----
      subroutine r_diag_normlz(A, Vc, w, M)
      implicit none
      
      integer, parameter :: flag_check=-1
      integer M, st, matz, ierr, i, h, j
      real*8 A(M, M), w(M), Vc(M,M), sumv, ratio
      real*8 delta; external delta
      real*8, dimension(:), allocatable :: fv1, fv2
      real*8, dimension(:, :), allocatable :: mat_dum
      
      integer, dimension(:), allocatable :: i_old_ofnew
      real*8, dimension(:),   allocatable :: w_dum
      real*8, dimension(:,:), allocatable :: Vc_dum
      
c----
c  first check symmetry
c---
      if(flag_check==1) then
      do i=1,M ; do j=1, M
         if(abs(A(i,j)-A(j,i))>1.d-6) then
!        if(abs(A(i,j)-A(j,i))>1.d-9) then
            write(6,*)' problem: matrix you called is not symmetric '
            write(6,*)' i, j, A ',i, j, A(i,j), A(j,i)
            stop
         endif
      enddo;     enddo
      endif
c-----
c  diagonalize
c----

      allocate (fv1(M), fv2(M), stat= st); if(st.ne.0) stop
      matz = 1
      call rs(M, M, A, w, matz, Vc, fv1, fv2, ierr)
      if(ierr /=0) then; write(6,*)' ie',ierr;stop;endif

c----
c Now normalize
c---
      do i=1, M
         Vc(:,i) = Vc(:, i)/(dot_product(Vc(:,i), Vc(:, i)))**0.5d0
      enddo
!
! return if not checking
!
      if(flag_check/=1) return

c-------
c check orthogonality
c------- 
      do i=1, M
         do j=1, M
            sumv = dot_product(Vc(:, i), Vc(:, j))
            if(abs(sumv-delta(i,j)).gt.1.d-8) then
               write(6,*) i, j, ' sumv ', sumv
               stop
            endif
         enddo
      enddo

c-----
c now order eigenvalues. First call ordering routine.
c----
      allocate (i_old_ofnew(M), w_dum(M), 
     x                         Vc_dum(M,M),stat=st); if(st.ne.0)stop

      call order_r_index(w,i_old_ofnew,M)
      w_dum = w
      Vc_dum = Vc
      do i=1, M
            w(i) = w_dum(   i_old_ofnew(i))
         Vc(:,i) = Vc_dum(:,i_old_ofnew(i))
      enddo

      deallocate(i_old_ofnew,Vc_dum,w_dum)

c----
c check that w is indeed the right value - can be changed
c----

      allocate (mat_dum(M, M),  stat= st); if(st.ne.0) stop      

      do i=1, M
         mat_dum(:, i) = Vc(:, i) * w(i)
      enddo


      ratio = sum( (mat_dum - matmul(A, Vc) )**2 )/ sum(w**2)
      if(ratio.gt.1.d-6) then
         write(6,*)' ratio ', ratio
         stop
      endif


      deallocate (fv1, fv2, mat_dum)
      
      end




