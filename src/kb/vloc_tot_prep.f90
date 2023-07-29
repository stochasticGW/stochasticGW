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
subroutine vloc_tot_prep
  use simple_mpi, only    : rank, bcast_r8
  use kb_mod
  implicit none
  integer               :: ma,nrp,ia,ig,iz,iy,ix,nxb,nyb,nzb
  integer, external     :: is_it_power2
  complex*16, parameter :: ci = (0.d0,1.d0)

  nxb = nx * scale_vh
  nyb = ny * scale_vh
  nzb = nz * scale_vh
  
  if(rank==0) then
     allocate(cvkb(nxb,nyb,nzb),cvn(nxb,nyb,nzb), stat=st); call check0(st,' cvkb cvn ')
     cvn = 0d0
  end if
  
  atomtypeloop : do ma=1,matop
     nrp = nrpp(ma)  ! note

     if(rank==0) then
        write(17,*)' kb pseudpt formfactor. atomtype(ma)=',ma,' charge= ',mapat(ma),' nrp ',nrp
        call flush(17)
     endif
     
     if(nrp<50.or.nrp>3000) stop ' nrp values problem '

     if(rank==0) then
        call calc_cvkb(nxb,nyb,nzb,nrp,rrpp(1,ma),vpploc(1,ma),cvkb)
        call check_real(cvkb, size(cvkb))
     endif

     call add_formfactor(ma,na,nxb,nyb,nzb)     ! cvn = cvn + form*wg_big*vk_big, for rank=0
  enddo atomtypeloop
  
  if(allocated(cvkb))deallocate(cvkb)

  allocate(vloc_tot(ng),stat=st);  call check0(st,' vloc_tot ')
  vloc_tot = 0d0
  if(rank==0) then
     call fft3d_general(  cvn,nxb,nyb,nzb,1)
     cvn = cvn/(dble(nxb)*dble(nyb)*dble(nzb)*dv)
     call check_c_real(   cvn,size(cvn),' cvn ')     ! after fft v should be real
     call c_to_v_compress(cvn, nxb,nyb,nzb,vloc_tot,nx,ny,nz)
  end if
 
!!! PTMOD added (to agree with dev vsn)
  if(dim_periodic.gt.0) then
!     if(rank==0) write(*,*) '??? shift vloc_tot by ',- sum(vloc_tot)/(dble(nx)*dble(ny)*dble(nz))   
     vloc_tot = vloc_tot - sum(vloc_tot)/(dble(nx)*dble(ny)*dble(nz))
  endif
!!! END PTMOD

  if(allocated(cvn))deallocate(cvn)

!  call bcast_r8(vloc_tot, ng, 0) !!! PTMOD removed since in vloc_tuma_prep

end subroutine vloc_tot_prep

subroutine  c_to_v_compress(cv, nxb,nyb, nzb, vv, nx,ny,nz)
  implicit none
  integer nxb, nyb, nzb
  integer nx, ny, nz
  integer mx, my, mz
  integer px, py, pz
  complex*16 cv(nxb, nyb, nzb)
  real*8     vv(nx,    ny,    nz)
  
  mx= 1 + nxb/2 - nx/2
  px=     mx+nx-1
  
  my= 1 + nyb/2 - ny/2
  py=     my+ny-1
  
  mz= 1 + nzb/2 - nz/2
  pz=     mz+nz-1
  
  vv = cv(mx:px,my:py,mz:pz)
end subroutine c_to_v_compress

