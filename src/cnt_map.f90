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
subroutine map_cnt
  use atoms,     only : ch_a
  use gwm,       only : na
  implicit none
  call map_cnt_1st  
contains
  subroutine map_cnt_1st
    use atoms,      only : an => atom_name, un=> atom_upper_name 
    use atoms,      only : matop,mapat,mapai
    use simple_mpi, only : rank
    implicit none
    integer, allocatable :: used_atom(:)
    integer              :: ma, ic, ia,st,rep

    if(allocated(mapai)) stop 'ERROR: mapai already allocated '
    allocate(mapai(na),stat=st); if(st/=0)stop ' mapai '

    allocate(used_atom(size(an)), stat=st); if(st/=0)stop ' used_atom  '
    used_atom = 0

    do ia=1,na
        call check_lele(1,ch_a(ia),size(an),' 1,charge(atom), n118 ')
       used_atom(ch_a(ia))= 1
    end do
    
    do rep=1,2
       ma = 0
       do ic=1,size(used_atom)
          if(used_atom(ic)==1) then
             ma = ma+1
             if(rep==2) mapat(ma) = ic
          endif
       enddo
       matop = ma
       if(rep==1) allocate(mapat(matop),stat=st); if(st/=0) stop ' mapat '
    end do

   do ia=1,na
      do ma=1,matop
         if(ch_a(ia)==mapat(ma)) mapai(ia) = ma
      end do
   end do
   
   if(rank==0) then
      write(6,*) ' Number of elements: ',matop
      write(6,*)
      write(17,*)' number of elements: ',matop
      write(17,*)' mapai ',mapai(1:na); call flush(9000)
      do ma=1,matop
         write(17,*)' number of atoms of charge ',&
              mapat(ma), ' is ',count(mapai==ma)
      enddo
      call flush(6)
      call flush(17)
   end if
   deallocate(used_atom)
 end subroutine map_cnt_1st

end subroutine map_cnt


     
