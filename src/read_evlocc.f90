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
subroutine read_evlocc(ip)
  use gwm
  use gwm, only : occall => occ_det_full
  use simple_mpi, only : rank, bcast_scalar_r8, bcast_scalar_i, bcast_r8

  implicit none
  integer  ip
  integer  st
  character*9 ch
  
  rnk0: if(rank==0) then
     if(binwf) then;    read(ip  )ch, ntop
     else;              read(ip,*)ch, ntop
     endif
     if(ch/='nstates') then; write(6,*)" ERROR reading wf.txt. Expected nstates,got ",ch; stop
     end if
     
     if(allocated(occall)) stop ' occall error '
     allocate(occall(ntop),stat=st); call check0(st,' occall ')
     !
     ! checking
     !
     if(ntop<rnel * nsp/2d0) then
        write(6,*)' not enough states available; ntstates(input) vs. total charge ',ntop,rnel
     end if
     
     !
     ! evl: alloc, determine
     !
     if(allocated(evl)) stop ' ERROR -- evl shouldnt be occupied yet '
     allocate(evl(ntop), stat=st); call check0(st,' evl ')
     
     if(binwf) then
        read(ip  )ch
        read(ip  )evl(:)
     else
        read(ip,*)ch
        read(ip,*)evl(:)
     endif
     
     if(ch/='evls') then
        write(6,*)" ERROR in reading wf.txt. Expected evls , got ",ch
        stop
     end if
     
     if(orb_indx>0) then
        call check_le(orb_indx,ntop,' orb_indx ntop ')
        eorb=evl(orb_indx)
        write(6,*)
        write(6,*)" eorb = eigenvalue( state ",orb_indx,") = ",real(eorb) 
        write(6,*)
     endif
     
     call setocc(evl,occall,rnel,ntop)
  end if rnk0
  
  !
  ! bcast
  !
  call bcast_scalar_r8(eorb)
  call bcast_scalar_i(ntop)
  
  if(.not.allocated(occall)) then 
     allocate(occall(ntop),stat=st); call check0(st,' occall ')
  endif
  call bcast_r8(occall, size(occall), 0)
  
  if(.not.allocated(evl)) then 
     allocate(evl(ntop),stat=st); call check0(st,' evl ')
  endif
  call bcast_r8(evl, size(evl), 0)

  call set_nrm_trace 
end subroutine read_evlocc

subroutine distribute_occ_det
  use simple_mpi, only : rank, bcast_scalar_r8, crnk => color_rank
  use gwm,        only : occall => occ_det_full
  use gwm
  implicit none
  integer is, isf, st
  real*8  oc1
  
  if(ns >0) then
     allocate(occ_det(ns),stat=st);   if(st/=0) stop "Occ_det allocation failed!"
  end if
  
  !call debug_distribute
  
  is = 0
  do isf=1,ntop  ! note: change later ntop 
     if(rank==0) oc1 = occall(isf)
     call bcast_scalar_r8(oc1)
     
     if(is_start(crnk).le.isf.and.isf.le.is_end(crnk)) then      
        is = is+ 1
        occ_det(is) = oc1
     end if
  end do
end subroutine distribute_occ_det

