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
subroutine calc_cvkb(nxb,nyb,nzb,nrp,rp,vr,cvkb)
  use kb_mod, only : dx,dy,dz,dv,dim_periodic
  implicit none
  integer, parameter  :: nf=16384 ! f=fine
  integer nxb, nyb,nzb,nrp,ix,iy,iz,ir,st
  real*8  rp(nrp),vr(nrp), rv(3),vrlg,rscalar,rlg
  real*8  pi,aa_here,charge_here
  real*8, allocatable, dimension(:) :: vr2, rplog, ka, g, g2
  complex*16 cvkb( nxb, nyb, nzb )

  pi=dacos(-1d0)

  select case(dim_periodic)
  case(2)
     call calc_cvkb_2d
  case(3,0)
     call calc_cvkb_3d_or_0d
  case default
     stop ' dim_periodic problem '
  end select

contains

  subroutine calc_cvkb_2d
    implicit none
    real*8,  parameter  :: factor_length=12d0
    integer             :: i, ir, st
    real*8              :: rmax, drf, r, v
    real*8              :: dk, k, b

    ! All vec. below
    ! v(k) = int exp(-ikr) v(r) d3r = v_s(k)+va(k)
    !
    ! trick: if r<L/2, v(r) is isotropic. But v(r)=va(r) at long dist.
    ! Here: va(r) = charge(negative) /r* theta(L/2-|z|)
    !
    ! So subtruct and add.
    !
    ! v_s(k) = v_s(|k|)= int_0^rmax (v(r)-va(r)) exp(-ik * r)d3r
    !        =           int_0^rmax (v(r)-va(r)) 2 pi db r^2 dr exp(-ikr b) ! b=cos(theta)
    !        =           int_0^rmax (v(r)-va(r)) 2 pi r^2/(kr) dr (exp(-ikr)-exp(ikr))/(-i) 
    !        =    4*pi/k int_0^rmax (v(r)-va(r)) r sin(kr) dr

    ! spline in rp, not its log
    allocate(vr2(nrp), stat=st); call check0(st,' vr2 ')
    call SPLINE_dn(rp,vr,nrp,vr2)
    charge_here = vr(nrp)*rp(nrp)

    !
    ! now we need the linear grid
    !
    allocate(g(nf), stat=st); call check0(st,' vf_nf ')
    rmax = rp(nrp)* factor_length
    dk=pi/rmax
    drf  = rmax/nf

    b = 0
    do i=1,nf
       r = (i-1)*drf

       if(r<rp(1)) then
          v = vr(1)
       elseif(r>rp(nrp)) then
          v = charge_here/r
       else
          call SPLINT_dn(rp,vr,vr2,nrp,r,v)
       end if

       g(i) = v*r - charge_here  ! = (v-va)*r
       b = b+ g(i)*r
    end do
    b = 4d0*pi*b*drf

    call sinfft_1d(g, nf)
    ! g(r)-->sum(sin(kr) g(r)), k=0,dk,..,(nf-1)*dk, and dk= pi/(nf*dr). Not 2pi!
    g = g*drf

    g(1)=b ! special care about k=0 point
    do i=1+1,nf  ! note: 1+1
       k = (i-1)*dk
       g(i) = 4d0*pi/k * g(i)
    enddo

    !
    ! make spline for g
    !
    allocate(ka(nf), g2(nf), stat=st); call check0(st,' ka g2 ')
    do i=1,nf
       ka(i) = (i-1)*dk
    enddo
    call SPLINE_dn(ka,g,nf,g2)

    call add_short_long_2d

    if(allocated(vr2))deallocate(vr2)
    if(allocated(g))  deallocate(g)
    if(allocated(g2)) deallocate(g2)
    if(allocated(ka)) deallocate(ka)

  end subroutine calc_cvkb_2d

  subroutine calc_cvkb_3d_or_0d
    implicit none

    if(dim_periodic.eq.3) then
       aa_here = max(7d0/min(nxb*dx,nyb*dy,nzb*dz), 5d0/rp(nrp))
       write(17,*)' aa_here ',aa_here
    end if

    !
    ! OK, now we have vr on a long grid.  Now we need to put in on a big 3-d grid.
    !  First, fit to spline

    allocate(vr2(nrp), rplog(nrp), stat=st); call check0(st,' vr2 ')

    do ir=1,nrp
       rplog(ir) = log(rp(ir))
    enddo

    call SPLINE_dn(rplog,vr,nrp,vr2)

    charge_here = vr(nrp)*rp(nrp)

    do iz=1,nzb
       do iy=1,nyb
          do ix=1,nxb
             rv = (/ (ix-1-nxb/2)*dx, (iy-1-nyb/2)*dy, (iz-1-nzb/2)*dz /)
             rscalar = max(1d-10,sqrt(sum(rv**2)))
             rlg = log(rscalar)

             if(rscalar<rp(1)) then
                vrlg = vr(1)
             elseif(rscalar>rp(nrp)) then
                vrlg = charge_here/rscalar
             else
                call SPLINT_dn(rplog,vr,vr2,nrp,rlg,vrlg) 
             end if
  
             ! if periodic you need to remove the long range part, since it will be 
             ! calculated for all images a few lines below.
             if(dim_periodic.eq.3) vrlg = vrlg - charge_here * erf(aa_here*rscalar)/rscalar 

             cvkb(ix,iy,iz) = vrlg
          end do
       end do
    end do

    ! now we have the 3d potential in r space. convert to k:
    call fft3d_forward_many( nxb,nyb, nzb,1,cvkb)
    cvkb = cvkb*dv
    call shift_cvkb
    if(dim_periodic.eq.3) call add_long_range_periodic_3d
    call check_real(cvkb, size(cvkb))
    deallocate(vr2,rplog)

  end subroutine calc_cvkb_3d_or_0d

  subroutine interp_lin(rplog,vr,nrp,rlg,vrlg)
    use kb_mod, only : dv
    implicit none
    integer nrp, i, j,k
    real*8    vr(nrp)
    real*8 rplog(nrp)
    real*8 rlg, vrlg
    i=1
    j=nrp 
    do while(j>i+1)
       k=(i+j)/2
       if(rlg<rplog(k)) then
          j=k
       else
          i=k
       endif
    end do
    vrlg = (vr(i)*(rplog(j)-rlg)+vr(j)*(rlg-rplog(i)))/(rplog(j)-rplog(i)) 
  end subroutine interp_lin

  subroutine shift_cvkb
    use kb_mod, only : dx, dy, dz
    implicit none
    integer ikx,iky,ikz
    real*8  dkx,dky,dkz,kx,ky,kz,r0(3)
    complex*16, parameter :: ci = (0d0,1d0)

    r0 = (/ -nxb/2d0*dx, -nyb/2d0*dy, -nzb/2d0*dz /)
    do ikz=1,nzb
       do iky=1,nyb
          do ikx=1,nxb
             dkx = 0d0 ; if(dx>1e-8) dkx = 2.d0*pi/(dble(Nxb)*dx)
             dky = 0d0 ; if(dy>1e-8) dky = 2.d0*pi/(dble(Nyb)*dy)
             dkz = 0d0 ; if(dz>1e-8) dkz = 2.d0*pi/(dble(Nzb)*dz)
             
             kx = dble(ikx-1)* dkx
             ky = dble(iky-1)* dky
             kz = dble(ikz-1)* dkz
             
             if(kx>pi/dx) kx = kx-2d0*pi/dx  
             if(ky>pi/dy) ky = ky-2d0*pi/dy
             if(kz>pi/dz) kz = kz-2d0*pi/dz
             
             cvkb(ikx,iky,ikz) = &
             cvkb(ikx,iky,ikz) * exp(-ci*( r0(1)*kx+r0(2)*ky+r0(3)*kz))
          enddo
       enddo
    end do
  end subroutine shift_cvkb

  subroutine add_long_range_periodic_3d
    use kb_mod, only : dx, dy, dz
    implicit none
    integer ikx,iky,ikz
    real*8  dkx,dky,dkz,kx,ky,kz,r0(3),k2
    complex*16, parameter :: ci = (0d0,1d0)

    r0 = (/ -nxb/2d0*dx, -nyb/2d0*dy, -nzb/2d0*dz /)
    do ikz=1,nzb
       do iky=1,nyb
          do ikx=1,nxb
             dkx = 0d0 ; if(dx>1e-8) dkx = 2.d0*pi/(dble(nxb)*dx)
             dky = 0d0 ; if(dy>1e-8) dky = 2.d0*pi/(dble(nyb)*dy)
             dkz = 0d0 ; if(dz>1e-8) dkz = 2.d0*pi/(dble(nzb)*dz)
             
             kx = dble(ikx-1)* dkx
             ky = dble(iky-1)* dky
             kz = dble(ikz-1)* dkz
             
             if(kx>pi/dx) kx = kx-2d0*pi/dx  
             if(ky>pi/dy) ky = ky-2d0*pi/dy
             if(kz>pi/dz) kz = kz-2d0*pi/dz
             
             k2 = kx**2+ky**2+kz**2
             if(k2>1d-8) then
                cvkb(ikx,iky,ikz) = &
                cvkb(ikx,iky,ikz) + &
                                  charge_here * 4d0*pi/k2*exp(-k2/4d0/aa_here**2)
             else
                cvkb(ikx,iky,ikz) = 0d0
             endif
          enddo
       enddo
    enddo

  end subroutine add_long_range_periodic_3d

  subroutine add_short_long_2d
    use simple_mpi, only : rank
    implicit none
    integer ikx,iky,ikz
    real*8  dkx,dky,dkz
    real*8  kx,ky,kz
    real*8  k2,k,zc,q, va, vd

    zc = nzb*dz/2d0
    dkx = 0d0 ; if(dx>1e-8) dkx = 2.d0*pi/(dble(nxb)*dx)
    dky = 0d0 ; if(dy>1e-8) dky = 2.d0*pi/(dble(nyb)*dy)
    dkz = 0d0 ; if(dz>1e-8) dkz = 2.d0*pi/(dble(nzb)*dz)

    cvkb=0d0
    do ikz=1,nzb
       do iky=1,nyb
          do ikx=1,nxb
             kx = dble(ikx-1)* dkx
             ky = dble(iky-1)* dky
             kz = dble(ikz-1)* dkz

             if(kx>pi/dx) kx = kx-2d0*pi/dx
             if(ky>pi/dy) ky = ky-2d0*pi/dy
             if(kz>pi/dz) kz = kz-2d0*pi/dz

             k2 = kx**2+ky**2+kz**2
             k  = sqrt(k2)

             if(k>1d-8) then
                q = sqrt(kx**2+ky**2)
                va = charge_here * 4d0*pi / k2 * ( 1d0-exp(-q*zc)*cos(kz*zc) )

                if(k.ge.ka(1).and.k.le.ka(nf)) then
                   call SPLINT_dn(ka,g,g2,nf,k,vd)
                else
                   stop ' ka problem '
                end if

                cvkb(ikx,iky,ikz)=va+vd
             endif
          enddo
       enddo
    enddo

  end subroutine add_short_long_2d

end subroutine calc_cvkb

