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
subroutine vk_set_exch
  use gwm, only : vk, vk_scale_cnst
  use gwm, only : vk_exch, vk_old
  use gwm, only : dim_periodic
  use gwm, only : method_exch
  use simple_mpi, only : rank
  implicit none
  integer, save :: i1=1

  if(i1/=1) return
  i1=-1

  call set_seed_general

  call alloc

  vk_old = vk
  if(vk_scale_cnst.le.0d0) stop ' ERROR: vk_scale_cnst not defined '

  if_periodic : if(dim_periodic.gt.0) then
  
     select case(method_exch)
     case ('no_fuzzy_no_0')
        
        vk_exch = vk
        vk_exch(1)=0d0     
        
     case('fuzzy','v00')
        
        select case(dim_periodic)
        case(2)
           call vk_exch_mc_2d_fuzzy
        case(3)
           call vk_exch_mc_3d_fuzzy
        case default
           stop ' vk_exch: fuzzy, but dim_periodic/=2,3 '
        end select

        vk_exch = vk_exch / size(vk_exch) ! note, DN
        if(method_exch=='v00') vk_exch(2:)=vk(2:) 
        
     case default
        stop ' method_exch error '
     end select

     if(rank==0) then
!        open(22,file=trim(adjustl(method_exch))//'.txt',status='replace')
!        write(22,*)size(vk_exch)
!        write(22,*)vk_exch
!        call flush(22)
!        close(22)
     end if

  else ! non-periodic
     vk_exch = vk
  end if if_periodic

!  if(rank==0) then
!     write(6,*)' vk(1:3), ',real(vk(1:3))
!     write(6,*)' vkx(1:3) ',real(vk_exch(1:3))
!  endif

contains
  subroutine alloc
    implicit none
    integer st
    if(.not.allocated(vk)) stop ' problem: vk not allocated '
    allocate(vk_exch(size(vk)),stat=st); call check0(st,' vk_exch ')
    allocate(vk_old( size(vk)),stat=st); call check0(st,' vk_old  ')
  end subroutine alloc
end subroutine vk_set_exch

