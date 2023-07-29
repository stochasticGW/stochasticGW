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
subroutine c_mtrx_exp(c_A, c_exp_A, N) 
  implicit none
  integer N, j, matz, jp, i
  integer, parameter :: fl_check_diag=1
  integer, parameter ::  fl_c_diag_method = 1

  real*8  diff, delta ; external delta

  complex*16 c_A(N,N), c_exp_A(N, N)

  complex*16, dimension(:),   allocatable :: c_VL
  complex*16, dimension(:,:), allocatable :: c_Vc, c_A_try, c_Vc_inv
  real*8    , dimension(:,:), allocatable :: unity

  integer isymm
  integer ifirst ; save ifirst; data ifirst /1 /

  if(maxval(abs(C_A))<1.d-10) then
     C_exp_A = 0.d0
     do i=1, N
        C_exp_A(i, i) = 1.d0
     enddo
     return
  endif

  ifirst = ifirst + 1
  matz = 1

  allocate(c_vc(N, N), c_vc_inv(N, N), c_VL(N), stat = j); if(j/=0) stop
  allocate(unity(N, N), c_A_try(N, N),          stat = j); if(j/=0) stop

  unity = 0.d0
  do i=1, N
     unity(i, i) = 1.d0
  enddo

  if(sum(abs(c_A - transpose(C_A)))<1.d-8) then
     isymm = 1
     call cs_diag_normlz(c_A, c_Vc, c_VL, N, matz, fl_c_diag_method)
     c_Vc_inv = transpose(c_Vc)
  else
     isymm = 0
     call cg_diag(C_A, C_Vc, C_VL, N, matz, 0)
     call c_inv_mtrx(c_Vc, c_Vc_inv, N)   
  endif
  check_d: if(fl_check_diag == 1) then
     
     diff =0
     do j=1, N
        diff = diff + sum(abs(matmul(c_A,c_Vc(:,j))-c_Vc(:,j)*c_vl(j)))
     enddo
     if(diff>1.d-8) then
        write(6,*)' diff (cA*Cv-Cv*cVl) ', diff
        stop
     endif
     
     diff = sum(abs(matmul(c_vc,c_vc_inv)-unity))
     if(diff > 1.d-8) then 
        write(6,*)' overalp diff ',diff
        stop
     endif
     
     !
     ! check more  can erase
     ! 
     do j=1, N
        c_A_try(:, j) = c_vC(:, j)* c_Vl(j)
     enddo

     c_A_try = matmul( c_A_try, c_Vc_inv)
     diff = sum(abs(c_A_try - c_A))
     
     if(diff > 1.d-8) then
        write(6,*)' diff (c_A_try - C_A))',diff
        stop
     endif
  else
     call check(fl_check_diag,-1,'fl_chk_d')
  endif check_d

  do j=1, N
     c_exp_A(:, j) = c_vC(:, j)* exp(c_Vl(j))
  enddo
  c_exp_A = matmul( c_exp_A, c_Vc_inv)

  deallocate(c_vc, c_vl, c_vc_inv, unity, c_A_try) 

end subroutine c_mtrx_exp

subroutine c_mtrx_exp_times_vec(c_A, c_exp_A, N, ket, ket_new, orb)
  implicit none
  integer N, j, matz, orb

  integer, parameter :: fl_check_diag = 1
  integer, parameter :: fl_c_diag_method = 1

  real*8  diff_A  
  complex*16 c_A(N,N), c_exp_A(N, N)
  complex*16 ket(N, Orb), ket_new(N, Orb)

  complex*16, dimension(:),   allocatable :: c_VL
  complex*16, dimension(:,:), allocatable :: c_Vc


  allocate(c_vc(N, N), c_VL(N), stat = j); if(j/=0) stop

  matz = 1
  call cs_diag_normlz(c_A, c_Vc, c_VL, N, matz, fl_c_diag_method)

  check_d: if(fl_check_diag == 1) then
     !
     ! check  can erase
     ! 
     do j=1, N
        c_exp_A(:, j) = c_vC(:, j)* c_Vl(j)
     enddo
     
     c_exp_A = matmul( c_exp_A, transpose( c_Vc))
     diff_A = sum(abs(c_exp_A - c_A))
     
     if(diff_A > 1.d-8) then
        write(6,*)' diff_A ',diff_A
        stop
     endif

  else
     call check(fl_check_diag,-1,'fl_chk_d')
  endif check_d

  ket_new = matmul(transpose(c_Vc),ket)

  do j=1, N
     ket_new(j, :) = exp(c_Vl(j))*ket_new(j,:)
  enddo

  ket_new =  matmul(c_Vc,ket_new)

  deallocate(c_vc, c_vl) 

end subroutine c_mtrx_exp_times_vec

subroutine cs_diag_normlz_modified(c_A, c_Vc, c_VL, N, matz)
  implicit none
  integer N, matz, j, ierr
  complex*16 c_A(N, N), c_Vc(N, N), c_Vl(N)
  complex*16, parameter :: eye=(0.d0,1.d0)
  

  real*8, dimension(:),   allocatable :: wr, wi, fv1, fv2, fv3
  real*8, dimension(:,:), allocatable :: Zr, Zi, ar, ai
  integer, parameter :: fl_check_diag = 1

  allocate(wr(N), wi(N), zr(N, N), zi(N, N), stat=j);  if(j/=0) stop
  allocate(              ar(N, N), ai(N, N), stat=j);  if(j/=0) stop
  allocate(fv1(N), fv2(N), fv3(N),           stat=j);  if(j/=0) stop

  matz = 1
  ar = c_A
  ai = aimag(C_A)
  write(6,*)' N ',N
  write(6,*)' Ar ',Ar
  write(6,*)' AI ',AI
  call cg(N,N,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
  call check(ierr, 0, ' ierr, zero ')

  c_vl =  wr+eye*wi   !cmplx(wr, wi)
  c_vc =  zr+eye*zi   !cmplx(zr, zi)

  
  select case (fl_check_diag)
  case (1)
     call c_grahm_s_fast(c_vc, N, N, 1.d0)
  case (-1)
     do j=1, N
        c_vc(:,j)= c_vc(:,j)/  sqrt(sum(c_vc(:,j)**2))
     enddo
  case default
     write(6,*)' problem, fl ',fl_check_diag; stop
  end select

  deallocate(ar, ai, wr, wi, zr, zi, fv1, fv2, fv3)
end subroutine cs_diag_normlz_modified

subroutine c_grahm_s_fast(CA, N, Nvec, Wght)  ! no checks, one loop only.
  implicit none                               ! when doubting, use c_grahm_s(
  integer N, Nvec, i, j
  complex*16 CA(N, Nvec)
  real*8 Wght
  

  do i=1, N
     do j=1, i-1
        CA(:, i) = CA(:, i) - sum(CA(:,i)*CA(:,j))*CA(:,j)
     enddo
     CA(:, i) = CA(:, i)/   sqrt(sum(CA(:,i)**2))
  enddo
end subroutine c_grahm_s_fast
  




