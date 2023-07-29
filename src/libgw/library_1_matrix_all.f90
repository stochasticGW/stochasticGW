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

!
! Matrix programs.
! Here are the two matrix programs that you'll need:
!  mat_inv   and mat_diag     as well as plot_mat
!
! To use these routines, place the statement      
!                   use mat_module
! in your program, immideately after the 'program' or 'subroutine' statement
! (before the implicit none).
! 
! Calling:
!
!
!  call mat_inv(A,A_inv)
! 
!        where A (input) is a general real*8 matrix, and so is A_inv (output).
!     OR:
!       A and A_inv can be both complex*16 matrices.
!
!       A and A_inv should be of the same size (N*N, where N is the # of
!       rows of A)
!
!
!
! call mat_diag(A,A_vc, A_evl)
!
!        where A (input) is either a SYMMETRIC real*8 matrix
!               or ANY complex*16 matrix
!              A_vc (ouput) is a (non-symmetric, but orthogonal) real*8 matrix
!                     (i.e., the eigenvector matrix)
!         and 
!              A_evl  (output) is a list of eigenvalues .
!
!        Note that A_evl is NOT the diagonal eigenvalue matrix, D,
!        defined as 
!                   D=0; do i=1,N; D(i,i)=A_evl(i); enddo
!
!        The eigenvalue equation is matmul(A,A_vc) = matmul(A_vc, D)
!
!        That equation is checked.
!
!
module mat_module
  interface mat_inv
     module procedure mat_inv_r8
     module procedure mat_inv_c16
  end interface
  
  interface mat_diag
     module procedure mat_diag_r8_symm
     module procedure mat_diag_c16
     module procedure mat_h_diag
  end interface

  interface plot_mat
     module procedure plot_mat_r8
     module procedure plot_mat_c16
  end interface

  interface mat_trace
     module procedure mat_trace_r8
     module procedure mat_trace_c16
  end interface
contains
  real*8 function mat_trace_r8(A)
    implicit none
    
    real*8 A(:,:)
    integer N, M, j

    N = Ubound(A, 1)- Lbound(A,1) + 1
    M = Ubound(A, 2)- Lbound(A,2) + 1
    
    if(M/=N) then
       write(6,*)' problem with sizes in A ',N,M
       write(6,*)' size(A) ',size(A)
       stop
    endif

    mat_trace_r8 = 0.d0
    do j=0,N-1
       mat_trace_r8 = mat_trace_r8 + A(Lbound(A,1)+j,Lbound(A,2)+j)
    enddo
       
  end function mat_trace_r8

  complex*16 function mat_trace_c16(A)
    implicit none
    
    complex*16 A(:,:)
    integer N, M, j

    N = Ubound(A, 1)- Lbound(A,1) + 1
    M = Ubound(A, 2)- Lbound(A,2) + 1
    
    if(M/=N) then
       write(6,*)' problem with sizes in A ',N,M
       write(6,*)' size(A) ',size(A)
       stop
    endif

    mat_trace_c16 = 0.d0
    do j=0,N-1
       mat_trace_c16 = mat_trace_c16 + A(Lbound(A,1)+j,Lbound(A,2)+j)
    enddo
       
  end function mat_trace_c16
    
  subroutine mat_inv_r8(A, A_inv)
    implicit none
    real*8 A(:,:), A_inv(:,:)
    real*8, allocatable ::  Unity(:,:)
    integer N, M, st, i
    integer, parameter :: flag_check = 1   
    real*8  diff

    N = Ubound(A, 1)- Lbound(A,1) + 1
    M = Ubound(A, 2)- Lbound(A,2) + 1
    
    if(M/=N.or.size(A_inv,1)/=N.or.size(A_inv,2)/=N) then
       write(6,*)' problem with sizes in mat_inv_r8_symm '
       write(6,*)' size(A) ',size(A)
       write(6,*)' size(A_inv) ',size(A_inv)
       stop
    endif

    call r_inv_mtrx(A, A_inv, N)
    if(flag_check /= 1) return

!
!   Pre checking: prepare unity
!
    allocate(unity(N,N), stat=st)
    if(st/=0) then
       write(6,*)' problem allocating diag_evl, unity in mat_inv '
       stop
    endif

    Unity    = 0.d0
    do i=1, N
       Unity(   i,i) = 1.d0
    enddo
!
!   Now check
!
    diff = maxval(abs(matmul(A, A_inv)-unity))/N

!    write(6,*)' deviation from A*Ainv = 1 ',diff
    if(diff>1.d-4) then
       write(6,*)' deviation too big -- stopping '
!       stop
    endif
!    write(6,*)
!
!    space clearing 
!    

    deallocate(Unity)
  end subroutine mat_inv_r8

  subroutine mat_inv_c16(A, A_inv)
    implicit none
    complex*16 A(:,:), A_inv(:,:)
    real*8, allocatable ::  Unity(:,:)
    integer N, M, st, i
    integer, parameter :: flag_check = 1   
    real*8  diff

    N = Ubound(A, 1)- Lbound(A,1) + 1
    M = Ubound(A, 2)- Lbound(A,2) + 1
    
    if(M/=N.or.size(A_inv,1)/=N.or.size(A_inv,2)/=N) then
       write(6,*)' problem with sizes in mat_inv_r8_symm '
       write(6,*)' size(A) ',size(A)
       write(6,*)' size(A_inv) ',size(A_inv)
       stop
    endif

    call c_inv_mtrx(A, A_inv, N)
    if(flag_check /= 1) return

!
!   Pre checking: prepare unity
!
    allocate(unity(N,N), stat=st)
    if(st/=0) then
       write(6,*)' problem allocating diag_evl, unity in mat_inv '
       stop
    endif

    Unity    = 0.d0
    do i=1, N
       Unity(   i,i) = 1.d0
    enddo
!
!   Now check
!
    diff = sum(abs(matmul(A, A_inv)-unity))/sum(abs(unity))

!    write(6,*)' deviation from A*Ainv = 1 ',diff
    if(diff>1.d-4) then
       write(6,*)' deviation c too big -- stopping '
!       stop
    endif
!    write(6,*)
!
!    space clearing 
!    

    deallocate(Unity)
  end subroutine mat_inv_c16

  subroutine mat_diag_r8_symm(A, A_vc, A_evl)
    implicit none
    real*8 A(:,:), A_vc(:,:), A_evl(:)
    real*8, allocatable :: Diag_evl(:,:), unity(:,:)
    real*8  diff
    integer N, M,st,i
    integer, parameter :: flag_check = -1   

!
!   sizes 
!
    N = Ubound(A, 1)- Lbound(A,1) + 1
    M = Ubound(A, 2)- Lbound(A,2) + 1

    if(M/=N.or.size(A_vc,1)/=N.or.size(A_vc,2)/=N.or.size(A_evl)/=N) then
       write(6,*)' problem with sizes in mat_diag_r8_symm '
       write(6,*)' size(A) ',size(A)
       write(6,*)' size(A_vc) ',size(A_vc)
       write(6,*)' size(A_evl) ',size(A_evl)
       stop
    endif
    
!
!   call a ready made program
!
    call r_diag_normlz(A, A_vc, A_evl, N)

    if(flag_check /= 1) return
!
!   pre-checking: prepare Diag_evl, Unity
!

    allocate(Diag_evl(N,N), unity(N,N), stat=st)
    if(st/=0) then
       write(6,*)' problem allocating diag_evl, unity in mat_diag '
       stop
    endif

    Unity    = 0.d0
    Diag_evl = 0.d0
    do i=1, N
       Diag_evl(i,i) = A_evl(i)
       Unity(   i,i) = 1.d0
    enddo
!
!   first checking the A_vc and A_evl are indeed e.vectors and e.values.
!

    diff = sum(abs(matmul(A, A_vc)-matmul(A_vc, Diag_evl)))/  &
           sum(abs(matmul(A, A_vc)))

!    write(6,*)' deviation from A*Avc=Avc*Diag ',diff
    if(diff>1.d-4) then
       write(6,*)' deviation a too big -- stopping '
!       stop
    endif
!    write(6,*)
!
!   now checking that A_vc is orthogonal (A_vc*A_vc_tranpose=I)
!

    diff = sum(abs(matmul(A_vc, transpose(A_vc))-unity))/sum(abs(unity))

!    write(6,*)' deviation from A_vc being orthogonal, i.e., '
!    write(6,*)' A_vc*A_vc_transpoe=I   is ',diff
        if(diff>1.d-4) then
       write(6,*)' deviation b too big -- stopping ',diff
!       stop
    endif
!    write(6,*)
!
!    space clearing 
!    

    deallocate(Diag_evl, Unity)
  end subroutine mat_diag_r8_symm

  subroutine mat_diag_c16(A, A_vc, A_evl)
    implicit none
    complex*16 A(:,:), A_vc(:,:), A_evl(:)
    real*8, parameter :: tol_h = 1.d-4
    real*8  diff
    integer N, M,st,i
    integer, parameter :: flag_check = 1   
    N = size(A,1)
    M = size(A,2)
!
!   sizes 
!
    N = Ubound(A, 1)- Lbound(A,1) + 1
    M = Ubound(A, 2)- Lbound(A,2) + 1

    if(M/=N.or.size(A_vc,1)/=N.or.size(A_vc,2)/=N.or.size(A_evl)/=N) then
       write(6,*)' problem with sizes in mat_diag_r8_symm '
       write(6,*)' size(A) ',size(A)
       write(6,*)' size(A_vc) ',size(A_vc)
       write(6,*)' size(A_evl) ',size(A_evl)
       stop
    endif
    
    diff = sum(abs(A-transpose(conjg(A))))/sum(abs(A)) 

!   write(6,*) ' diff from in mat_diag_c16 from hermit. ', diff

    if(diff<tol_h) then
      ! write(6,*)' calling hermitian diag '
       call ch_diag(A, A_vc, A_evl, M)
    else
       !write(6,*)' calling general (non-hermitian) diag '
       call cg_diag(A, A_Vc, A_evl, M, 1, flag_check)
    endif
  end subroutine mat_diag_c16

  subroutine mat_h_diag(A, A_vc, r_evl)
    implicit none
    real*8     r_evl(:)
    complex*16 A(:,:), A_vc(:,:)
    complex*16, allocatable :: a_evl(:)
    real*8, parameter :: tol_h = 1.d-8
    real*8  diff
    integer N, M,st,i
    N = size(A,1); M=size(A,2)
    if(M/=N.or.size(A_vc,1)/=N.or.size(A_vc,2)/=N.or.size(r_evl)/=N) then
       write(6,*)' problem mat_h_diag;  sizes of A, A_vc, r_evl ',size(A), size(A_vc), size(r_evl);    stop
    endif
    allocate(a_evl(N), stat=st); call check(st,0,' a_evl ')
    
    diff = sum(abs(A-transpose(conjg(A))))/sum(abs(A)) !;   write(6,*) ' diff in mat_h from hermit. ', diff
    if(diff>tol_h) stop

    !write(6,*)' pre calling of ch_diag '; call flush(6)
    call ch_diag(A, A_vc, A_evl, N);
    !write(6,*)' called ch_diag  in mat_h_diag '; call flush(6)
    r_evl = A_evl
    deallocate(a_evl)
    !write(6,*)' finished mat_h_diag ';call flush(6)
  end subroutine mat_h_diag

  subroutine ch_diag(A, A_vc, A_evl, m)

    implicit none
    integer m, ierr,i
    integer, parameter :: icheck=1
    real*8,  parameter :: tol_h_diag = 1.e-8 ! effective only if icheck=1
    real*8             :: diff
    real*8, allocatable :: ar(:,:), ai(:,:), w(:), &
         zr(:,:), zi(:,:), fv1(:), fv2(:),fm1(:,:), Unity(:,:)

    complex*16, parameter :: ci = (0.d0,1.d0)
    complex*16 A(m,m), A_vc(m,m), A_evl(m)

    allocate(ar(m,m), ai(m,m), w(m), zr(m,m),&
         zi(m,m), fv1(m),fv2(m),fm1(2,m), Unity(m,m), stat=ierr)
    call check(ierr,0,' stch ')

    Unity = 0.d0
    do i=1, m
       Unity(i,i) = 1.d0
    enddo

    ar = dble( A)
    ai = aimag(A)

    call ch(m,m,ar,ai,w,1,zr,zi,fv1,fv2,fm1,ierr)

    A_evl = w
    A_vc  = zr + ci*zi

    if(icheck==1) then
       diff = 0.d0
       do i=1, m
          diff = diff + sum(abs(matmul(A,A_vc(:,i))-A_evl(i)*A_vc(:,i)))
       enddo
       diff = diff/sum(abs(A))

 !       write(6,*)' V V^H- 1 ',sum(abs(matmul(A_vc,transpose(conjg(A_vc)))-Unity))
!       write(7788,*)diff; 
!       call flush(7788)

       if(diff>tol_h_diag) then
          write(6,*)' dev. from ch*v - v ', diff
          write(6,*)' V V^H- 1 ',sum(abs(matmul(A_vc,transpose(conjg(A_vc)))-Unity))
          write(6,*)' minval(Avl),minval(abs(Avl)),maxval(Avl) ',&
               minval(w),minval(abs(w)),maxval(w)
          stop
       endif
    end if

    deallocate(ar, ai, w, zr, zi, fv1,fv2,fm1, Unity)
  end subroutine ch_diag

  subroutine plot_mat_c16(A, iplot)
    implicit none
    complex*16 A(:,:)
    integer N,M, i, iplot
    N = size(A,1)
    M = size(A,2)
    if(N<1.or.M<1.or.N*M>1000) then
       write(6,*)' not set to plot matrix with dimensions ',N,M
       stop
    endif
    
    do i=lbound(A,1), ubound(A,1)
       write(iplot,77)A(:,i)
       if(i>3)write(iplot,*)
77     format(' ',3(f11.6,'+i*',f11.6,' '))
    enddo
  end subroutine plot_mat_c16

  subroutine plot_mat_r8(A, iplot)
    implicit none
    real*8 A(:,:)
    integer N,M, i, iplot
    N = size(A,1)
    M = size(A,2)
    if(N<1.or.M<1.or.N*M>1000) then
       write(6,*)' not set to plot matrix with dimensions ',N,M
       stop
    endif
    
    do i=lbound(A,1), ubound(A,1)
       write(iplot,77)A(:,i)
       if(i>6)write(iplot,*)
77     format(' ',6f12.6)
    enddo
  end subroutine plot_mat_r8
end module

       








