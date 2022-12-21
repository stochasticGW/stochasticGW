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
subroutine c_orthog_ket_wrt_bra(vec, ket_f, n, A2, w_int, c_det)
  implicit none
  integer n, A2, i 
  real*8  w_int, diff
  complex*16 c_det, vec(n, A2), ket_f(n, A2)
  complex*16, dimension(:,:), allocatable :: s_ov, s_inv, unity

  integer fl_check_diag;   common/flag_check_com/fl_check_diag  

  call check(abs(fl_check_diag),1,'abs_fl_check_diag')
  allocate(s_ov(A2,A2), s_inv(A2,A2), unity(A2,A2),stat=i);if(i/=0) stop
  unity = 0.d0; do i=1,A2; unity(i,i) = 1.d0; enddo
  s_ov  = matmul(transpose(vec),ket_f)*w_int

  call c_inv_mtrx_with_det(s_ov, s_inv, A2, c_det)

  ket_f = matmul(ket_f, s_inv)

!!  check_diag : if(fl_check_diag == 1) then
     diff = sum(abs(matmul(s_ov,s_inv) - unity))
     if(diff > 1.d-5) then
        write(6,*)' diff_orthog_a '
        write(6,*)' product ',matmul(s_ov, s_inv)
        write(6,*)' diff ',diff
        stop
     endif

     diff = sum(abs(matmul(transpose(vec),ket_f)*w_int - unity))
     if(diff > 1.d-5) then
        write(6,*)' diff_orthog_b ',diff
        write(6,*)' matmul(transpose(vec),ket_f)*w_int '
        write(6,*)  matmul(transpose(vec),ket_f)*w_int
        stop
     endif

!!  endif check_diag

  deallocate(s_ov, s_inv, unity)
end subroutine c_orthog_ket_wrt_bra

subroutine c_norm_bos_orth_fer_core(ket,N,A2,Orb,Orb_core,w_int,fl_b,c_det)
  implicit none
  integer N, A2, fl_b, fl_grd, orb, orb_core
  complex*16 ket(N, A2), c_det
  real*8 w_int


  if(fl_b == 1) then
     ! normalize for bosons
     c_det = c_det *  sqrt(sum(ket**2)*w_int/A2)**A2
     ket = ket /      sqrt(sum(ket**2)*w_int/A2)  ! exact nrmlz frml not impo
  else if(fl_b == -1) then
     call c_grahm_s_core_det(ket, N, A2, Orb, Orb_core, w_int, c_det)
  else
     write(6,*)' oops. fl_b = ',fl_b; stop
  endif

end subroutine c_norm_bos_orth_fer_core

subroutine c_norm_bos_or_orthonorm_fer(ket, N, A2, w_int, fl_b, c_det)
  implicit none
  integer N, A2, fl_b, fl_grd
  complex*16 ket(N, A2), c_det
  real*8 w_int

  if(fl_b == 1) then
     ! normalize for bosons
     c_det = c_det *  sqrt(sum(ket**2)*w_int/A2)  
     ket = ket /      sqrt(sum(ket**2)*w_int/A2)  ! exact nrmlz frml not impo
  else if(fl_b == -1) then
     call c_grahm_s_det(ket, N, A2, w_int, c_det)
  else
     write(6,*)' oops. fl_b = ',fl_b; stop
  endif

end subroutine c_norm_bos_or_orthonorm_fer

subroutine  c_grahm_s(ket, N, A2, w_int)
  implicit none
  integer N, A2, irepeat, j, i
  real*8 w_int, delta
  complex*16 ket(N, A2), c_term

  do irepeat=1, 4   ! it's important to repeat twice.

     do i=1, A2
        do j=1, i-1
           ket(:,i) = ket(:,i) - sum(ket(:,i)*ket(:,j))*w_int*ket(:,j)
        enddo

        ket(:,i) = ket(:,i) /   sqrt(sum(ket(:,i)*ket(:,i))*w_int)
     enddo     
  enddo

  do i=1, A2
     do j=1, A2
        c_term = sum(ket(:,i)*ket(:,j))*w_int
        if(abs(c_term - delta(i,j))>1.d-10) then
           write(6,*) ' in c_grahm_s , i, j , c_term ',i, j, c_term
           stop
        endif
     enddo
  enddo
  
end subroutine c_grahm_s

subroutine  c_grahm_s_det(ket, N, A2, w_int, c_det)
  implicit none
  integer N, A2, irepeat, j, i
  real*8 w_int, delta
  complex*16 ket(N, A2), c_term, c_det

  c_det = 1.d0

  do irepeat=1,4   ! it's important to repeat twice.

     do i=1, A2
        do j=1, i-1
           ket(:,i) = ket(:,i) - sum(ket(:,i)*ket(:,j))*w_int*ket(:,j)
        enddo

        c_term =   sqrt(sum(ket(:,i)*ket(:,i))*w_int)
        c_det = c_det       * c_term
        ket(:,i) = ket(:,i) / c_term
     enddo     
  enddo

  do i=1, A2
     do j=1, A2
        c_term = sum(ket(:,i)*ket(:,j))*w_int
        if(abs(c_term - delta(i,j))>1.d-10) then
           write(6,*) ' in c_grahm_s , i, j , c_term ',i, j, c_term
           stop
        endif
     enddo
  enddo
  
end subroutine c_grahm_s_det

subroutine  c_grahm_s_core_det(ket, N, A2, Orb, Orb_core, w_int, c_det)
  implicit none
  integer N, A2, irepeat, j, i, Orb, Orb_core
  real*8 w_int, delta
  complex*16 ket(N, A2), c_term, c_det

  c_det = 1.d0

  if(Orb>N.or.A2>Orb.or.Orb_core>A2) then
     write(6,*)' Orb_core, A2, Orb,  A2 '
     write(6,*)  Orb_core, A2, Orb,  A2
     stop
  endif

  do irepeat=1,4   ! it's important to repeat twice.

     do i=1, Orb
        do j=1, min(i-1, Orb_core)
           ket(:,i) = ket(:,i) - sum(ket(:,i)*ket(:,j))*w_int*ket(:,j)
        enddo

        if(i <= Orb_core) then
           c_term =   sqrt(sum(ket(:,i)*ket(:,i))*w_int)
           c_det = c_det       * c_term
           ket(:,i) = ket(:,i) / c_term
        endif
     enddo     
  enddo

  do i=1, Orb
     do j=1, Orb_core
        c_term = sum(ket(:,i)*ket(:,j))*w_int
        if(abs(c_term - delta(i,j))>1.d-10) then
           write(6,*) ' in c_grahm_s , i, j , c_term ',i, j, c_term
           stop
        endif
     enddo
  enddo
  
end subroutine c_grahm_s_core_det

subroutine r_norm_bos_or_orthonorm_fer(ket, N, A2, w_int, fl_b)
  implicit none
  integer N, A2, fl_b, fl_grd
  real*8 ket(N, A2)
  real*8 w_int
  
  if(fl_b == 1) then

     ! normalize for bosons
     
     ket = ket / sqrt(sum(ket**2)*w_int/A2)  ! exact nrmlz frml not impo
     
  else if(fl_b == -1) then
  
     call r_grahm_s(ket, N, A2, w_int)

  else
     write(6,*)' oops. fl_b = ',fl_b; stop
  endif

end subroutine r_norm_bos_or_orthonorm_fer

subroutine  r_grahm_s(ket, N, A2, w_int)
  implicit none
  integer N, A2, irepeat, j, i
  real*8 w_int, delta
  real*8 ket(N, A2), r_term
  
  do irepeat=1, 2   ! it's important to repeat twice.

     do i=1, A2
        do j=1, i-1
           ket(:,i) = ket(:,i) - sum(ket(:,i)*ket(:,j))*w_int*ket(:,j)
        enddo
        ket(:,i) = ket(:,i) / sqrt(sum(ket(:,i)*ket(:,i))*w_int)
     enddo
     
  enddo

  do i=1, A2
     do j=1, A2
        r_term = sum(ket(:,i)*ket(:,j))*w_int
        if(abs(r_term - delta(i,j))>1.d-10) then
           write(6,*) ' in r_grahm_s , i, j , r_term ',i, j, r_term
           stop
        endif
     enddo
  enddo
  
end subroutine r_grahm_s

