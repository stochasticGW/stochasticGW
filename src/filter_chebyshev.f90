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
subroutine filter_chebyshev(isp, vks, n, ms, p, normz, dh, havg, th, nc, nchb, dv, rank)
  use gwm, only : scratch_dir
  implicit none
  integer n, ms, nc, nchb, rank, ik, k, i
  integer isp
  integer, external :: is_it_power2
  real*8  vks(n), normz(ms), dh, havg, th(0:nchb), dv
  complex*16 p(n*ms)
  complex*16, allocatable, dimension(:) :: po, pe, tmp

  if(rank==0)then;write(17,*)' entering filter; nc = ',nc;call flush(17);
  endif
  
  allocate(po(n*ms), pe(n*ms), tmp(n*ms),stat=i); if(i/=0) stop ' i  '   
     
  pe = p
  po = pe
  p = 0d0
  if(rank==0) then
    write(17,*) "sum p,po,pe - entering=",sum(p),sum(po),sum(pe)
    call flush(17)
  end if

  do k=0,nc-1
     if(mod(k,2)==0) then
        p=p+pe*th(k)
        call hc_scld_comp(isp, pe, po, tmp, vks, n, ms, havg, dh)  ! 1,0 --> 1,2
        if(k==0) po = (po+pe)/2d0
     else
        p=p+po*th(k)
        call hc_scld_comp(isp, po, pe, tmp, vks, n, ms, havg, dh)  
     end if
     if(is_it_power2(k)==1.and.rank==0) then; 
        write(17,*)' k in filter, |pfiltered| ',k,real(sum(abs(pe)))
        call flush(17) 
     end if
  enddo
  call normlz_p(p)
  deallocate(po,pe,tmp)
contains
  subroutine normlz_p(p)
    implicit none
    integer  is
    complex*16 p(n,ms)
    if(rank==0) then
       write(17,*) "shape,normpt=",shape(normz); 
       write(17,*) "sum p=",sum(p)
       call flush(17)
    endif

    do is=1,ms
       normz(is) = sqrt(sum(abs(p(:,is))**2)*dv)
       p(:,is)=p(:,is)/normz(is)
    enddo
    
    if(rank==0) then
       write(17,*) "normlz=",normz; call flush(17) ; 
    end if
  end subroutine normlz_p
end subroutine filter_chebyshev
