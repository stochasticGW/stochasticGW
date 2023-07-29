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
subroutine read_orb_1
  use gwm,        only : orb=>orb_rd, orb_kind, sp0, nsp, eorb
  use gwm,        only : nx, ny, nz, dx, dy, dz, n, dv
  use simple_mpi, only : bcast_scalar_i
  use simple_mpi, only : bcast_scalar_r8
  use simple_mpi, only : rank,bcast_r8
  implicit none
  logical                 :: fex
  character*30            :: a
  character*9             :: ch
  integer                 :: st, ip, msp
  complex*16, allocatable :: uc(:)
  
  ip = 111

  rnk0write : if(rank==0) then
     write(6,*)
     write(6,*) " ############# READING DFT ORBITAL FROM orb.txt #############"
     inquire(file='orb.txt',exist=fex);     if(.not.fex) stop 'orb.txt missing'
     open(ip,file='orb.txt',status='old');   rewind(ip)
     a=' from file: orb.txt '
  end if rnk0write

  call read_xyz(ip)

  if(.not.allocated(orb)) then; allocate(orb(n),stat=st); if(st/=0) stop ' orb problems '
  else; call check(size(orb),n,' size_orb, n ')
  endif

  rnk0rd : if(rank==0) then
     select case(nsp)
     case(1)
        sp0 = 1
        write(6,*)' Reading orbital energy (overriding eorb from INPUT if given there ) '
        read(ip,*)ch, eorb
     case(2)
        write(6,*)' Reading orbital energy and spin (overriding eorb from INPUT if given there ) '
        read(ip,*)ch, eorb, sp0
     case default
        stop ' ERROR: nsp should be 1 or 2 '
     end select

     if(ch/='orb')then; write(6,*)' ERROR,expect orb , got ',ch,a; stop 
     endif
     call check_lele(1,sp0,nsp,' one-sp0-nsp ')

     select case(orb_kind)
     case(1);  read(ip,*) orb
     case(2)
        allocate(uc(n),stat=st);call check0(st,' uc ') 
        read(ip,*) uc
        orb = dble(uc)
        deallocate(uc)
     case default; write(6,*)' ERROR: orb_kind not 1(=real*8) nor 2(=complex*16) but ',orb_kind
        call flush(6); stop
     end select

     write(17,*) 
     write(17,*) "-> sum(psi^2)",sum(orb**2)*dv
     if(abs(sum(orb**2)*dv -1.d0)>1d-7) then
        write(6,*) " WARNING: the orb.txt orbital is not correct -- as it isnt normalized !"
        write(6,*) " Specifcially, its integral =",sum(orb**2)*dv
        write(6,*) " The calculation is continued but... Good luck!" 
        call flush(6)
     end if
     call flush(17)

     close(ip)
  end if rnk0rd
  call bcast_scalar_i(   sp0   )
  call bcast_scalar_r8(  eorb  )
  call bcast_r8(orb,size(orb),0)

end subroutine read_orb_1

