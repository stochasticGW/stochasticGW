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
subroutine setocc(evl,occ,rnel,ntop)
  use gwm, only : tpin => tp_dns, nsp, toll=>tollocc
  implicit none
  integer i, in, ntop
  real*8 evl(ntop)
  real*8 occ(ntop)
  real*8 rnel
  real*8 mu
  real*8 tp
  real*8 etop
  real*8 ebot
  real*8 de, ddh
  
  integer, parameter :: indxmx = 2000
  real*8,  parameter :: tpmin = 1d-5

  tp = tpin

  if(tp<tpmin) then
     write(6,*)' Rescaling the density temperature (in Hartree) to be ',tpmin
     tp = tpmin
  endif

  if(rnel> 2d0/nsp*ntop ) stop ' not enough dft states '
  

  etop = maxval(evl)
  ebot = minval(evl)
  ddh   = etop - ebot
  etop = etop + 0.05*ddh 
  ebot = ebot - 0.05*ddh


  in = 0
  muloop : do
     in = in+1
     mu=(etop+ebot)/2d0
     
     occ = 0d0
     do i=1,ntop
        de = (evl(i)-mu)/tp
        if(de<-100) then
           occ(i) = 2d0/dble(nsp)
        else if(de<100) then
           occ(i) = 2d0/dble(nsp) / (1d0+exp(de))
        endif
     enddo
     
     if(abs(sum(occ)-rnel)<toll) exit muloop
     if(in> indxmx) stop ' too much time searching for mu '

     if(sum(occ)>rnel) then 
        etop=mu
     else
        ebot=mu
     end if
  end do muloop
  
  write(17,*)' Chemical potential for deterministic density = ',mu
  write(17,*)' Summed occupations          ',sum(occ)
  write(17,*)' summed occupations - charge ',sum(occ)-rnel

end subroutine setocc

  
