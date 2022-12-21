!                     
!                     
!    The routine(s) in this file are a part of the  
!                     StochasticGW                  
!    suite, developed 2012-2018, and copyrighted    
!    to the authors of StochasticGW , see:          
!                                                   
!                 www.stochasticgw.com              
!                                                   
!   If you use or modify any part of this routine   
!   the header should be kept, unmodified.          
!                                                   
!                                                   
!                                                   
      subroutine cdriver_diag()  ! dummy subroutine; checks cdiag
      implicit none
      integer M; parameter(M=3)
      integer i, j, flag_c_diag_method
      complex*16 A(M,M), Vc(M,M), w(M), U(M,M), A_dum(M, M)
      complex*16, parameter :: ci=(0.d0,1.d0)

      do i=1, M
         do j=1, M
            A(i, j) = i*2+j*2-i**2*10-j**2*10+(i-j)**4 
     x              +  ci*(i+j)**3+ci*(i-j)**4
         enddo
      enddo


      write(6,*)' A '
      write(6,88)A

      write(6,*)' A '
      write(6,88)A

      write(6,*)' method '
      read( 5,*)  flag_c_diag_method

      call cs_diag_normlz(A, Vc, w, M, 1, flag_c_diag_method)

      write(6,*)' A '
      write(6,88) A
      write(6,*)' w '
      write(6,88) w
      write(6,*)' Vc '
      do i=1, M
      write(6,88) Vc(i,:) 
      enddo

 88   format(1x,6f12.4)
      end

      subroutine cs_easy_diag_normlz(A, Vc, w, M)
      implicit none
      
      integer M
      integer fl_c_diag_method
      integer matz

      complex*16 A(M,M), Vc(M, M), w(M)
      
      fl_c_diag_method = 1
      matz = 1
      
      call cs_diag_normlz(A, Vc, w, M, matz, fl_c_diag_method)

      end

c-----
c generic subroutine to diagonalize A(M,M); results fulfil A*vc = eps*A. 
c                       also: each vector in vc is orthonormalized.
c                       and: eigenvalues and eigenvectors are ordered.
c     A is complex_symmetric
c----
      subroutine cs_diag_normlz(A, Vc, w, M, matz, flag_c_diag_method)
      implicit none
      
      integer M, matz,iflg_order, iflg_check,
     x                             fl_check_diag,flag_c_diag_method, st

      complex*16 A(M, M), Vc(M, M), w(M)
      real*8, allocatable:: Ar(:,:), Vr(:,:), wr(:)

      fl_check_diag = 1
      iflg_check = fl_check_diag 


      
      if(sum(abs(aimag(A)))>1.d-10*sum(abs(A))) then
         select case (flag_c_diag_method)
         case (1)
            call cs_usual_diag_normlz(A, Vc, w, M, matz, iflg_check)
         case (2)
            !!         call cs_new_smpl(A,Vc,w,M,matz,iflg_check)
            write(6,*)' problem in cs_new ';stop
         case (3)
            call cs_diag_normlz_modified(A, Vc, w, M, matz)
         case default
            write(6,*)' flag_c_diag_method '
            write(6,*)  flag_c_diag_method
            stop
         end select
      
      else
         allocate(wr(M), Ar(M,M), Vr(M,M), stat= st)
         call check(st,0,' stwr         ')
         Ar = A
         call r_diag_normlz(Ar, Vr, wr, M)
         Vc = Vr
         w  = wr
         deallocate(wr, Ar, Vr)
      endif

      end
      subroutine cs_usual_diag_normlz(A, Vc, w, M, matz, iflg_check)
      implicit none

      complex*16, parameter :: ci=(0.d0,1.d0)
      
      integer M, st, matz, ierr, i, h, j, iflg_order, iflg_check
      complex*16 A(M, M), w(M), Vc(M,M)


      real*8 sumv, ratio, maxa
      real*8 delta; external delta
      real*8, dimension(:), allocatable :: fv1, fv2, fv3, wr, wi
      real*8, dimension(:, :), allocatable :: ar,ai,zr,zi
      
      integer, dimension(:),      allocatable :: i_old_ofnew
      complex*16, dimension(:),   allocatable :: w_dum
      complex*16, dimension(:,:), allocatable :: Vc_dum, mat_dum

      iflg_order=1

!!    write(6,*)' matz ',matz
      if(matz/=0.and.matz/=1) stop
c----
c  first check symmetry
c---
      maxa = maxval(abs(A))
      do i=1,M ; do j=1, M
         if(abs(A(i,j)-A(j,i))>(1.d-6)*maxA) then
            write(6,*)' i, j, A ',i, j, A(i,j), A(j,i)
            stop
         endif
      enddo;     enddo
      A= (A+transpose(A))/2.d0

c-----
c  diagonalize
c----
      allocate (fv1(M), fv2(M), fv3(M))
      allocate(ar(M,M), stat=st);if(st/=0)stop
      allocate(ai(M,M), stat=st);if(st/=0)stop
      allocate(wr(M), stat=st);if(st/=0)stop
      allocate(wi(M), stat=st);if(st/=0)stop
      
      ar(:, :) = A(:, :);      ai(:, :) = aimag(A(:, :))
      if(matz/=0)allocate(zr(M,M),zi(M,M),stat=st);if(st/=0)stop
      
      call cg(M,M,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)

      if(ierr /=0) then; write(6,*)' ierr',ierr;stop;endif
      deallocate(ar,ai)
      deallocate(fv1, fv2, fv3)
      
      w = wr + ci*wi

      if(matz/=0)then; Vc = zr + ci*zi; deallocate(zr, zi);endif
c-----
c now order eigenvalues. First call ordering routine.
c----
      if(iflg_order == 1) then

         allocate (i_old_ofnew(M), w_dum(M))
         w_dum = w

         if(matz/=0)then
            allocate(  Vc_dum(M,M),stat=st); if(st.ne.0)stop
            Vc_dum = Vc
         endif

         call order_r_index(wr,i_old_ofnew,M)
         
         do i=1, M; w(i) = w_dum(   i_old_ofnew(i)); enddo
         if(matz/=0)then
            do i=1,M;Vc(:,i)=Vc_dum(:,i_old_ofnew(i));enddo
            deallocate(vc_dum)
         endif
         deallocate(i_old_ofnew,w_dum)
      endif

      deallocate(wr, wi)

      if(matz==0)return
c----
c Now orthogonalize
c---
      call c_grahm_s(Vc, M, M, 1.d0)
c-------
c check orthogonality. Careful when w's degenerate (ortho first).
c------- 

      if(iflg_check == 1) then
         do i=1, M
            do j=1, M
               sumv = sum(Vc(:, i)*Vc(:, j))
               if(abs(sumv-delta(i,j)).gt.1.d-5) then
                  write(6,*) i, j, ' sumv ', sumv
                  stop
               endif
            enddo
         enddo
c---- 
c     check that w is indeed the right value - can be changed
c---- 
         allocate (mat_dum(M, M),  stat= st); if(st.ne.0) stop      
         
         do i=1, M
            mat_dum(:, i) = Vc(:, i) * w(i)
         enddo
         
         ratio = sum( abs(mat_dum - matmul(A, Vc) ))/ sum(abs(w))
         if(ratio.gt.1.d-6) then
            write(6,*)' ratio ', ratio
            stop
         endif
         
         deallocate (mat_dum)
      endif
      
      end

c-----
c generic subroutine to diagonalize A(M,M); results fulfil A*vc = eps*A. 
c                       and: eigenvalues and eigenvectors are ordered.
c----put matz = 1 for eigenvectors; 0 otherwise

      subroutine cg_diag(A, Vc, w, M, matz, fl_check_diag)
      implicit none

      complex*16, parameter :: ci=(0.d0,1.d0)
      
      integer M, st, matz, ierr, i, h, j, iflg_order, iflg_check
      integer fl_check_diag
      complex*16 A(M, M), w(M), Vc(M,M)
      real*8 sumv, ratio
      real*8 delta; external delta
      real*8, dimension(:), allocatable :: fv1, fv2, fv3, wr, wi
      real*8, dimension(:, :), allocatable :: ar,ai,zr,zi
      integer, dimension(:),      allocatable :: i_old_ofnew
      complex*16, dimension(:),   allocatable :: w_dum
      complex*16, dimension(:,:), allocatable :: Vc_dum, mat_dum

      iflg_order=1
      iflg_check = fl_check_diag


!!    write(6,*)' matz ',matz
      if(matz/=0.and.matz/=1) stop
c-----
c  diagonalize
c----
      allocate (fv1(M), fv2(M), fv3(M))
      allocate(ar(M,M), stat=st);if(st/=0)stop
      allocate(ai(M,M), stat=st);if(st/=0)stop
      allocate(wr(M), stat=st);if(st/=0)stop
      allocate(wi(M), stat=st);if(st/=0)stop
      
      ar(:, :) = A(:, :);      ai(:, :) = aimag(A(:, :))
      if(matz/=0)allocate(zr(M,M),zi(M,M),stat=st);if(st/=0)stop
      
      call cg(M,M,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)

      if(ierr /=0) then; write(6,*)' ierr',ierr;stop;endif
      deallocate(ar,ai)
      deallocate(fv1, fv2, fv3)
      
      w = wr + ci*wi
      
      if(matz/=0)then; Vc = zr + ci*zi; deallocate(zr, zi);endif
c-----
c     now order eigenvalues. First call ordering routine.
c---- 
      if(iflg_order == 1) then
         
         allocate (i_old_ofnew(M), w_dum(M))
         w_dum = w
         
         if(matz/=0)then
            allocate(  Vc_dum(M,M),stat=st); if(st.ne.0)stop
            Vc_dum = Vc
         endif
         
         call order_r_index(wr,i_old_ofnew,M)
         
         do i=1, M; w(i) = w_dum(   i_old_ofnew(i)); enddo
         if(matz/=0)then
            do i=1,M;Vc(:,i)=Vc_dum(:,i_old_ofnew(i));enddo
            deallocate(vc_dum)
         endif
         deallocate(i_old_ofnew,w_dum)
         deallocate(wr, wi)
         
         if(matz==0)return
      endif  !check
c---- 
c     check that w is indeed the right value - can be changed
c---- 
      if(iflg_check == 1) then
         allocate (mat_dum(M, M),  stat= st); if(st.ne.0) stop      
         
         do i=1, M
            mat_dum(:, i) = Vc(:, i) * w(i)
         enddo
         
         ratio = sum( abs(mat_dum - matmul(A, Vc) ))/ sum(abs(w))
         if(ratio.gt.1.d-6) then
            write(6,*)' ratio ', ratio
            stop
         endif
         
         deallocate (mat_dum)
      endif
      
      end

c------------------------------------------------------------------------------------------
c     modified cs_diag_normlz routine for a generic (non-symmetric) complex matrix
c     Ken Lopata May 21, 2008
c------------------------------------------------------------------------------------------
      subroutine cg_diag_new(A, Vc, w, M, matz)
      implicit none
      
      complex*16, parameter :: ci=(0.d0,1.d0)
      
      integer M, st, matz, ierr, i, h, j, iflg_order
      complex*16 A(M, M), w(M), Vc(M,M)
      
      
      real*8 sumv, ratio
      real*8 delta; external delta
      real*8, dimension(:), allocatable :: fv1, fv2, fv3, wr, wi
      real*8, dimension(:, :), allocatable :: ar,ai,zr,zi
      
      integer, dimension(:),      allocatable :: i_old_ofnew
      complex*16, dimension(:),   allocatable :: w_dum
      complex*16, dimension(:,:), allocatable :: Vc_dum, mat_dum
      
      iflg_order=1
      
!     !    write(6,*)' matz ',matz
      if (matz/=0.and.matz/=1) then
         stop
      endif
      
c-----
c     diagonalize
c---- 
      allocate (fv1(M), fv2(M), fv3(M))
      allocate(ar(M,M), stat=st);if(st/=0)stop
      allocate(ai(M,M), stat=st);if(st/=0)stop
      allocate(wr(M), stat=st);if(st/=0)stop
      allocate(wi(M), stat=st);if(st/=0)stop
      
      ar(:, :) = A(:, :);      ai(:, :) = aimag(A(:, :))
      if(matz/=0) then 
         allocate(zr(M,M),zi(M,M),stat=st)
      endif
      if(st/=0) then 
         stop
      endif
      
      call cg(M,M,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
      
      if(ierr /=0) then
         write(6,*)' ierr',ierr;stop
      endif

      deallocate(ar,ai)
      deallocate(fv1, fv2, fv3)
      
      w = wr + ci*wi
      
      if(matz/=0) then
         Vc = zr + ci*zi
         deallocate(zr, zi)
      endif
c-----
c     now order eigenvalues. First call ordering routine.
c---- 
      if(iflg_order == 1) then
         
         allocate (i_old_ofnew(M), w_dum(M))
         w_dum = w
         
         if(matz/=0)then
            allocate(  Vc_dum(M,M),stat=st); if(st.ne.0)stop
            Vc_dum = Vc
         endif
         
         call order_r_index(wr,i_old_ofnew,M)
         
         do i=1, M; w(i) = w_dum(   i_old_ofnew(i))
         end do
         
         if(matz/=0)then
            do i=1,M
               Vc(:,i)=Vc_dum(:,i_old_ofnew(i))
            enddo
            deallocate(vc_dum)
         endif
         deallocate(i_old_ofnew,w_dum)
      endif
         
      deallocate(wr, wi)

      if(matz==0) then 
         return
      endif
      
c---- 
c     check that w is indeed the right value - can be changed
c---- 
      allocate (mat_dum(M, M),  stat= st)
      if(st.ne.0) then 
         stop      
      endif
            
      do i=1, M
         mat_dum(:, i) = Vc(:, i) * w(i)
      enddo
            
      ratio=sum(abs(mat_dum-matmul(A, Vc) ))/ sum(abs(w))
      
      if(ratio.gt.1.d-6) then
         write(6,*)' ratio ', ratio
         stop
      endif
            
      deallocate (mat_dum)
      end subroutine

