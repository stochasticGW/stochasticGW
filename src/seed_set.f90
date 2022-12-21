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
subroutine set_seed_pt(ib)
  use simple_mpi,    only : color_size, color_rank, rank
  use gwm,           only : icolor, ncolors, ns, ns_blk, imc, is_map, nsp
  use gwm,           only : icounter
  use seed_gw_modu,  only : line=>line_seed
  use seed_gw_modu,  only : nseed
  use seed_gw_modu,  only : line_seed_pt
  use seed_gw_modu,  only : seed_array
  use seed_gw_modu,  only : seed_read
  use seed_gw_modu,  only : check_seed
  implicit none
  integer         :: ib
  integer         :: is_m, imc_m, i
  integer(kind=8) :: j1, j2, j3, j4, j5, j6, j7

  is_m  = is_map(ib, color_rank)
  call check_le(1,is_m,' one is_m ')
  if(nsp==2) is_m = (is_m-1) /nsp + 1  
  imc_m = icolor+(imc-1)*ncolors
  
  line= line_seed_pt 
  call check_seed

  j1 = is_m
  j2 = 42626
  j3 = imc_m
  j4 = 88935
  j6 = icounter
  j7 = 3827163
  do i=1,nseed
     j5 = i
     seed_array(i,line) = seed_read(i,line) + j1*j2 + j3*j4*j5 + j6*j7
  enddo
  call ran_ps_putseed(seed_array(:,line))
end subroutine set_seed_pt

subroutine set_seed_gam(igam)
  use simple_mpi,    only : color_size, color_rank
  use gwm,           only : icolor, ncolors, imc, igam_map
  use gwm,           only : icounter
  use seed_gw_modu,  only : line=>line_seed
  use seed_gw_modu,  only : nseed
  use seed_gw_modu,  only : line_seed_gam
  use seed_gw_modu,  only : seed_array
  use seed_gw_modu,  only : seed_read
  use seed_gw_modu,  only : check_seed
  implicit none
  integer igam_m, imc_m, igam, i
  integer(kind=8) :: j1, j2, j3, j4, j5, j6, j7
  igam_m  = igam_map(igam, color_rank) 
  imc_m   = icolor + (imc-1)*   ncolors
  
  line= line_seed_gam
  call check_seed

  j1 = igam_m
  j2 = 933273
  j3 = imc_m
  j4 = 866941
  j6 = icounter
  j7 = 3827163
  do i=1,nseed
     j5 = i
     seed_array(i,line) = seed_read(i,line) + j1*j2 + j3*j4*j5 + j6*j7
  enddo
  call ran_ps_putseed(seed_array(:,line))
end subroutine set_seed_gam

subroutine set_seed_zeta()
  use simple_mpi,    only : color_size 
  use gwm,           only : icolor, ncolors, imc
  use gwm,           only : icounter
  use seed_gw_modu,  only : line=>line_seed
  use seed_gw_modu,  only : nseed
  use seed_gw_modu,  only : line_seed_zeta
  use seed_gw_modu,  only : seed_array
  use seed_gw_modu,  only : seed_read
  use seed_gw_modu,  only : check_seed
  implicit none
  integer         :: imc_m, i
  integer(kind=8) :: j1, j2, j3, j6, j7
  
  call check_seed

  imc_m = icolor + (imc-1)* ncolors
  
  line= line_seed_zeta

  j1 = imc_m
  j2 = 54623
  j6 = icounter
  j7 = 3827163
  do i=1,nseed
     j3 = i
     seed_array(i,line) = seed_read(i,line) + j1*j2*j3 + j6*j7
  enddo
  call ran_ps_putseed(seed_array(:,line))
end subroutine set_seed_zeta

subroutine set_seed_g()
  use simple_mpi,    only : color_size
  use gwm,           only : icolor, ncolors, imc
  use gwm,           only : icounter
  use seed_gw_modu,  only : line=>line_seed
  use seed_gw_modu,  only : nseed
  use seed_gw_modu,  only : line_seed_g
  use seed_gw_modu,  only : seed_array
  use seed_gw_modu,  only : seed_read
  use seed_gw_modu,  only : check_seed
  implicit none
  integer         :: imc_m, i
  integer(kind=8) :: j1, j2, j3, j6, j7

  call check_seed

  imc_m = icolor + (imc-1)* ncolors
  line= line_seed_g
  j1 = imc_m
  j2 = 76367
  j6 = icounter
  j7 = 3827163
  do i=1,nseed
     j3 = i
     seed_array(i,line) = seed_read(i,line) + j1*j2*j3 + j6*j7
  enddo
  call ran_ps_putseed(seed_array(:,line))
end subroutine set_seed_g

subroutine set_seed_general()
  use gwm,           only : icounter
  use simple_mpi,    only : rank
  use seed_gw_modu,  only : line=>line_seed
  use seed_gw_modu,  only : nseed
  use seed_gw_modu,  only : line_seed_general
  use seed_gw_modu,  only : seed_array
  use seed_gw_modu,  only : seed_read
  use seed_gw_modu,  only : check_seed
  implicit none
  integer         :: i
  integer(kind=8) :: j1, j2, j3, j6, j7
  
  call check_seed

  line= line_seed_general

  j1 = rank
  j2 = 285901
  j6 = icounter
  j7 = 3827163
  do i=1,nseed
     j3 = i
     seed_array(i,line) = seed_read(i,line) + j1*j2*j3 + j6*j7
  enddo
  call ran_ps_putseed(seed_array(:,line))
end subroutine set_seed_general

subroutine set_seed_exchange(is)
  use simple_mpi,    only : color_size, color_rank
  use gwm,           only : icolor, ncolors, imc, ns, is_map
  use gwm,           only : icounter
  use seed_gw_modu,  only : line=>line_seed
  use seed_gw_modu,  only : nseed
  use seed_gw_modu,  only : line_seed_exchange
  use seed_gw_modu,  only : seed_array
  use seed_gw_modu,  only : seed_read
  use seed_gw_modu,  only : check_seed
  implicit none
  integer         :: imc_m, i, is, is_m 
  integer(kind=8) :: j1, j2, j3, j4, j5, j6, j7

  is_m  = is_map(is, color_rank)
  imc_m = icolor + (imc-1)*      ncolors
  
  line= line_seed_exchange
  call check_seed

  j1 = imc_m 
  j2 = 68571
  j3 = imc_m
  j4 = 87375
  j6 = icounter
  j7 = 3827163
  do i=1,nseed
     j5 = i
     seed_array(i,line) = seed_read(i,line) + j1*j2+ j3*j4*j5 + j6*j7
  enddo
  call ran_ps_putseed(seed_array(:,line))
end subroutine set_seed_exchange












