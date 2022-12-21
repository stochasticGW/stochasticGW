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
!
! called only fromn rank=0
!
subroutine read_wf1(ur, n, ip, isf, isp)
  use gwm, only : ntop, binwf, nx, ny, nz, dx, dy, dz, nsp
  implicit none
  integer        :: n, ip, isf, isp, i, m, mtop, msp
  integer        :: mx, my, mz
  real*8         :: ddx,ddy,ddz
  real*8         :: ur(n)
  character*9    :: ch

  isf1: if(isf==1) then

     if(binwf) then;  
        read(ip  )ch, mx;  read(ip  ) ch, my;  read(ip  ) ch, mz; call checkm
        read(ip  )ch, ddx; read(ip  ) ch, ddy; read(ip  ) ch, ddz; call checkd
        read(ip  )ch, msp; call check(nsp,msp,' msp, nsp ')
        read(ip  )ch, mtop; if(ch/='nstates') then; write(6,*)" ERROR wf nstates ",ch;stop; endif
        call check(ntop,mtop,' nmtop ')
        read(ip  )ch;  read(ip  ); if(ch/='evls')then;write(6,*)" ERROR wf evls ,",ch;stop; endif
        read(ip  )ch;    if(ch/='orbitals') then; write(6,*)" ERROR wf orbitals ",ch; stop; endif
     else;              
        read(ip,*)ch, mx;  read(ip,*) ch, my;  read(ip,*) ch, mz; call checkm
        read(ip,*)ch, ddx; read(ip,*) ch, ddy; read(ip,*) ch, ddz; call checkd
        read(ip,*)ch, msp; call check(nsp,msp,' msp, nsp ')
        read(ip,*)ch, mtop; 
        if(ch/='nstates') then; 
           write(6,*)" ERROR wf.txt. nstates ",ch; stop; 
        endif
        call check(ntop,mtop,' nmtop ')
        read(ip,*)ch;read(ip,*);if(ch/='evls')then; write(6,*)" ERROR wf evls  ,",ch;  stop; endif
        read(ip,*)ch;if(ch/='orbitals')     then; write(6,*)" ERROR wf orbitals  ",ch; stop; endif
     end if
  end if isf1

  if(binwf)then;read(ip  )i,m;call check(i,isf,'isf');call check(m,isp,'msp');read(ip  ) ur(:)
  else;         read(ip,*)i,m;call check(i,isf,'isf');call check(m,isp,'msp');read(ip,*) ur(:)
  endif

contains 
  subroutine checkm
    implicit none
    call check(nx,mx,' nx, mx ')
    call check(ny,my,' ny, my ')
    call check(nz,mz,' nz, mz ')
  end subroutine checkm
  
  subroutine checkd
    implicit none
    call check_r(ddx,dx,' ddx, dx ')
    call check_r(ddy,dy,' ddy, dy ')
    call check_r(ddz,dz,' ddz, dz ')
  end subroutine checkd

end subroutine read_wf1


