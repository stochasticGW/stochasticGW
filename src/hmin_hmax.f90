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
subroutine h_minmax
  use gwm, only : n, dv, dhscl, dh, havg, hmin, hmax
  use gwm, only : vks, imc, tp, nchbmx
  use gwm, only : powr_maxit, powr_ftol
!  use gwm, only : functional, de_bnl
  use simple_mpi, only : rank, bcast_scalar_r8
  implicit none
  integer, save :: i1=1
  integer       :: ii!, isp
  real*8        :: ee,ee2
  real*8        :: oldrayqt,conver
  complex*16    :: c21,c12
  complex*16,allocatable :: pa(:),pb(:),p1(:),tmp(:)
  if(i1/=1)   return
  i1=-1

!  isp=1
  ! need to be done over all runks due to bcast in random.inp reading
  allocate(pa(n),pb(n),p1(n),tmp(n))
  call set_seed_general()
  call rand_c(pa,n,dv)
  call rand_c(pb,n,dv)
  
  p1=0d0
  tmp=0d0
  
  pa=pa/sqrt(sum(abs(pa)**2)*dv)
  pb=pb/sqrt(sum(abs(pb)**2)*dv)
  
  p1=0d0; call hc(pa, p1, vks, n, 1)
  c21 = sum(conjg(pb)*p1)*dv
  
  p1=0d0; call hc(pb, p1, vks, n, 1)
  c12 = sum(conjg(pa)*p1)*dv
  
  if(rank==0) then
     write(17,*) "--- checking Hermicity of the Hamiltonian applied here:"
     write(17,*) "<psi1|h*psi2> and <psi2|h*psi1>",c12,c21
     write(17,*)
     call flush(17)
  end if

  if(abs(c21-conjg(c12))>1d-6) then
     write(6,*)' ERROR in hermiticity; c12, c21 = ',c12,c21
     call flush(6)
     stop
  endif

  !Now - power iter.
  oldrayqt=0.d0
  pa=0.d0; pb=0.d0;tmp=0.d0
  call rand_c(pa,n,dv)
  pa=pa/sqrt(sum(abs(pa)**2)*dv)
  
  do ii=0,powr_maxit
     p1=0d0; call hc(pa, p1, vks, n, 1)
     ee = sum(conjg(pa)*p1*dv)/sum(conjg(pa)*pa*dv)

     if (mod(ii,10).eq.0) then
        conver=abs((oldrayqt-ee)/ee)
        oldrayqt=ee
        if(rank==0) then
          write(17,'(X,A,X,I5,X,A,X,f15.8,X,A,X,ES10.3,X)') &
          'itn:',ii,'ee_extreme',ee,'fdel',conver; call flush(17)
        endif
     endif

     pa=p1/sqrt(sum(abs(p1)**2)*dv)
     if (conver.le.powr_ftol) exit
  end do

  if(rank==0) then
     write(17,*)' ee_extreme ',ee; call flush(17)
  endif
  
  ! now the other extreme eigenvalue
  oldrayqt=0.d0
  pa=0d0
  call rand_c(pa,n,dv)
  pa=pa/sqrt(sum(abs(pa)**2*dv))
     
  do ii=0,powr_maxit
     p1=0d0; call hc(pa, p1, vks, n, 1)
     p1 = p1- ee*pa
     ee2 = sum(conjg(pa)*p1*dv)/sum(conjg(pa)*pa*dv)

     if (mod(ii,10).eq.0) then
        conver=abs((oldrayqt-ee2)/ee2)
        oldrayqt=ee2
        if(rank==0) then
           write(17,'(X,A,X,I5,X,A,X,f15.8,X,A,X,ES10.3,X)') &
           'itn:',ii,'ee_other',ee2,'fdel',conver; call flush(17)
        endif
     endif

     pa=p1/sqrt(sum(abs(p1)**2)*dv)
     if (conver.le.powr_ftol) exit
  enddo

  hmin = min(ee,ee+ee2)
  hmax = max(ee,ee+ee2)

!  if(functional=='bnl') hmin = hmin + de_bnl 

  if(rank==0) then
     write(17,*)' ee_other ',ee2; call flush(17)
!     write( 6,*)
!     write( 6,*)" ############ HAMILTONIAN INFO: ############ "
!     write( 6,*)
!     write(17,*)" ############ HAMILTONIAN INFO: ############ "
!     write(17,*)' extracted minmax  ',ee,ee+ee2
!     write(17,*)' Extracted hmin, hmax        ',real(hmin),real(hmax)

!     havg = (hmax+hmin)/2d0
!     dh   = (hmax-hmin)/2d0
!     dh = dh*dhscl
!
!     write(6,'(X,A,F16.8)')' hmin       = ',real(hmin)
!     write(6,'(X,A,F16.8)')' hmax       = ',real(hmax)
!     write(6,'(X,A,F16.8)')' dh scaling = ',real(dhscl)
!     write(6,'(X,A,F16.8)')' dh         = ',real(dh) 
!     write(6,'(X,A,F16.8)')' havg       = ',real(havg)
!     write(17,*)' havg=',real(havg),'; dh(scaled by ',real(dhscl),') is ',real(dh)
!     call flush(6)  
!     call flush(17) 
  endif
  deallocate(pa,pb,p1,tmp)

!  call bcast_scalar_r8(havg)
!  call bcast_scalar_r8(dh)

!  if (nchbmx==-1) nchbmx = dh/tp*30d0  !note new formula
!  if(rank==0) then
!     write(17,*)' nchbmx = ',nchbmx
!     call flush(17)
!     write(6,'(X,A,I16)')' nchbmx     = ',nchbmx
!     call flush(6)
!  endif
end subroutine h_minmax

subroutine proc_hminmax_info

  use gwm, only : dhscl, dh, havg, hmin, hmax
  use gwm, only : tp, nchbmx
  use simple_mpi, only : rank, bcast_scalar_r8
  implicit none

  if(rank==0) then
     write( 6,*)
     write( 6,*)" ############ HAMILTONIAN INFO: ############ "
     write( 6,*)
     write(17,*)" ############ HAMILTONIAN INFO: ############ "
     write(17,*)' Extracted hmin, hmax        ',real(hmin),real(hmax)

     havg = (hmax+hmin)/2d0
     dh   = (hmax-hmin)/2d0
     dh = dh*dhscl

     write(6,'(X,A,F16.8)')' hmin       = ',real(hmin)
     write(6,'(X,A,F16.8)')' hmax       = ',real(hmax)
     write(6,'(X,A,F16.8)')' dh scaling = ',real(dhscl)
     write(6,'(X,A,F16.8)')' dh         = ',real(dh)
     write(6,'(X,A,F16.8)')' havg       = ',real(havg)
     write(17,*)' havg=',real(havg),'; dh(scaled by ',real(dhscl),') is ',real(dh)
     call flush(6)
     call flush(17)
  endif

  call bcast_scalar_r8(havg)
  call bcast_scalar_r8(dh)

  if (nchbmx==-1) nchbmx = dh/tp*30d0  !note new formula
  if(rank==0) then
     write(17,*)' nchbmx = ',nchbmx
     call flush(17)
     write(6,'(X,A,I16)')' nchbmx     = ',nchbmx
     call flush(6)
  endif

end subroutine proc_hminmax_info

