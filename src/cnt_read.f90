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
 subroutine read_cnt
  use atoms,     only : an => atom_name, un=> atom_upper_name
  use atoms,     only : cnt, ch_a, valch_a
  use gwm,       only : na, fu=> flg_units_nuc, angstrom_in_au, use_unified_input
  use simple_mpi,only : rank, bcast_i, bcast_r8, bcast_scalar_i

  implicit none
  integer     :: i,j,ip,st
  real*8      :: x,y,z
  character*2 :: ccha
  logical     :: found

  ! Skip if unified read is used
  if (use_unified_input) return

  ip=19 ! file number
  i=0
  if(rank==0) then
     if(fu==1)write(6,*)' Positions in bohr, read from: cnt.ini '
     if(fu==2)write(6,*)' Positions in Angstrom, read from: cnt.ini '

     if(fu==1)write(17,*)' Positions in bohr, read from: cnt.ini '
     if(fu==2)write(17,*)' Positions in Angstrom read from: cnt.ini '

     call flush(6)
     open(ip,file='cnt.ini',status='old')
     rewind(ip)
     do 
        read(ip,*,end=11) ccha
        i= i+1
     enddo
11   rewind(44)
     na = i
     write(6,*) ' Number of atoms :   ',na
     write(17,*)' Number of atoms :   ',na
   end if
  call bcast_scalar_i(na)
  
!  allocate(cnt(3,na),ch_a(na),valch_a(na),stat=st); call check0(st,' cnt_allocate ') !!! PTMOD valch_a removed below
  allocate(cnt(3,na),ch_a(na),stat=st); call check0(st,' cnt_allocate ')

  rnk0 : if(rank==0) then
     rewind(ip)
     i=0
     atom_read : do
        read(ip,*,end=99) ccha,x,y,z
        i=i+1
        cnt(:,i) = (/x,y,z/) 
        found = .false.
        atom_type : do j=1,size(an)
           if(ccha==an(j).or.ccha==un(j)) then
              if(found)stop ' read-atom-label fits more than 1 atom-types '
              found=.true.
              ch_a(i) = j
           endif
        end do atom_type
        if(.not.found) then
           write(6,*)' read-atom-label = ',ccha,' does not fit any element '
           stop
        end if
     enddo atom_read
     close(ip)

     select case(fu)
     case(1)
        cnt = cnt
     case(2) 
        cnt = cnt* angstrom_in_au
     case default
        stop ' Error: flg_units_nuc wrong'
     end select
  end if rnk0
99 continue

  call bcast_i( ch_a, size(ch_a), 0)
  call bcast_r8(cnt,size(cnt),0)

  if(rank==0) then
     write(17,*)' Atomic charges in molecule: ',ch_a
     call flush(6)
     call flush(17) 
  end if
end subroutine read_cnt
