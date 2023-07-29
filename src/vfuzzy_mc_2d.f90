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
subroutine  vk_exch_mc_2d_fuzzy
  use gwm, only : vk_exch,dx,dy,dz,nx,ny,nz
  use simple_mpi, only : rank
  implicit none
  logical, parameter :: q0=.true.
  integer,      save :: ip=17
  integer, parameter :: nmc=1024*1024*64 
  integer, parameter :: mca_floor = 10
  real*8,  parameter :: falloff = 3d0 ! put 3, not 4

  integer nkx, nky, nkz, nk
  integer ikx, iky, ikz, ik
  integer i
  integer imc, nksamp
  integer, allocatable :: pntrk(:),mca(:)

  real*8 ax, ay, az, zc
  real*8 kx, ky, kz
  real*8 area, pi, qxm, qym, qtop
  real*8 qx, qy, q
  real*8 phi, p
  real*8 px, py
  real*8 qpx, qpy ! qprime
  real*8 c
  real*8 tn, trm, savg, ds, s2
  real*8 uv(3), up, u2
  real*8 ratiok
  real*8, external :: ran_ps_neg  ! -1 to 1
  real*8, external :: ran_ps      !  0 to 1
  real*8, allocatable :: k2a(:)

  pi = dacos(-1d0)

  if(.not.allocated(vk_exch))stop ' vk_exch not allocated '
  if(size(vk_exch)<300) stop ' change 300 below '

  call make_cnsts
  call make_k2a
  call order_k2a
  call make_mca

  vk_exch = 0d0
  call set_vk_exch_det
  s2 = 0d0
  do imc=1,maxval(mca)
     if(q0) then
        qx=0d0
        qy=0d0
     else
        qx = ran_ps_neg()* qxm
        qy = ran_ps_neg()* qym
     end if

     phi      = pi*ran_ps_neg() ! -pi to pi
     p        = ran_ps()* 2d0* qtop  ! 0 to 2qtop
     
     px = p*cos(phi)
     py = p*sin(phi)

     qpx = px+qx
     qpy = py+qy

     bzcheck : if(abs(qpx).le.qxm.and.abs(qpy).le.qym) then
        kloop: do i=1,nk
           ik=pntrk(i)
           if(imc>mca(ik)) exit kloop
           call def_kxyz
           
           !uv  = (/ kx+qpx, ky+qpy, kz /)  ! negative comb
           uv  = (/ kx+px, ky+py, kz /)  ! corrected mistake
           up  = sqrt(uv(1)**2+uv(2)**2)
           u2  = up**2+kz**2
           trm = c* (1d0-exp(-zc*up)*cos(zc*kz))/u2*p 
           
           vk_exch(ik) = vk_exch(ik)+ trm
           if(ik==1)s2 = s2 + trm**2
        end do kloop
     end if bzcheck

     call plot

  end do
  
  call nrmlz_vk_exch
  call dealloc
  if(rank==0) write(ip,*)' vk_exch(k=0)   ',vk_exch(1)
  if(rank==0) write(ip,*)' nk, n          ',dble(size(vk_exch)),nx*ny*nz
  if(rank==0) write(ip,*)' vk_exch(k=0)/n ',vk_exch(1)/dble(size(vk_exch))
  if(rank==0) write(ip,*)' vk_exch(k=0)/V ',vk_exch(1)/(nx*ny*nz*dx*dy*dz)
  if(rank==0) write(ip,*)'vkx(1,2,30,300))         imc ',&
        real(vk_exch(1)),real(vk_exch(2)),real(vk_exch(30)),real(vk_exch(300))
  if(rank==0) call flush(ip)

contains
  subroutine make_cnsts
    implicit none

    nkx=nx
    nky=ny
    nkz=nz
    nk = nkx*nky*nkz
    call check(nk,size(vk_exch),' nk, size(vk_exch) ')

    ax = nkx*dx
    ay = nky*dy
    az = nkz*dz
    
    ! defs
    
    zc = az/2d0
    qxm = pi/ax
    qym = pi/ay
    area = ax*ay
    qtop = sqrt(qxm**2+qym**2)
    c = 4d0* area * qtop
  end subroutine make_cnsts

  subroutine set_vk_exch_det
    implicit none
    real*8 kv(3), kp, k2
    do i=1,nk
       ik=pntrk(i)
       if(mca(ik)==0) then
          call def_kxyz
          kv  = (/ kx, ky, kz /)
          kp  = sqrt(kv(1)**2+kv(2)**2)
          k2  = kp**2+kz**2
          vk_exch(ik) = 4d0*pi* (1d0-exp(-zc*kp)*cos(zc*kz))/k2 
       end if
    enddo
  end subroutine set_vk_exch_det

  subroutine def_kxyz
    implicit none
    real*8 dkx, dky, dkz
    
    call unravel(ik,ikx,iky,ikz)

    dkx = 2d0*pi/ax
    dky = 2d0*pi/ay
    dkz = 2d0*pi/az

    kx = dkx*(ikx-1)
    ky = dky*(iky-1)
    kz = dkz*(ikz-1)
    
    if(kx>pi/dx) kx=kx-2d0*pi/dx
    if(ky>pi/dy) ky=ky-2d0*pi/dy
    if(kz>pi/dz) kz=kz-2d0*pi/dz
  end subroutine def_kxyz

  subroutine unravel(ik,ikx,iky,ikz)
    implicit none
    integer ik, ikx, iky, ikz, i
    !    ik-1 = ikx-1 + (iky-1)*nkx + (ikz-1)*nkx*nky
    ikx = 1+mod(ik-1,nkx) !
    i   = (ik-ikx)/nkx    ! i= iky-1 + (ikz-1)*nky
    iky = 1+mod(i,nky)
    ikz = 1+i/nky
    
    call check(ik,ikx+(iky-1)*nkx+(ikz-1)*nkx*nky,' ik, comb ')
  end subroutine unravel

  subroutine make_k2a
    implicit none
    integer st
    allocate(k2a(nk), stat=st); call check0(st,' k2a ')

    ik=0
    do ikz=1,nkz
       do iky=1,nky
          do ikx=1,nkx
             ik=ik+1
             call def_kxyz
             k2a(ik)=kx**2+ky**2+kz**2
          enddo
       enddo
    enddo

    k2a(1) = 0.4*minval(k2a(2:))
    ! for the k=0 point, <1/k^2> is somewhat (details unimportant) larger than the k=(dkx,0,0)
    ! point or another such pount.  This is used only for ordering.
  end subroutine make_k2a

  subroutine order_k2a
    implicit none
    integer st
    allocate(pntrk(nk), stat=st);    call check0(st,' pntrk ')
    call order_r_index(k2a,pntrk,size(k2a))
  end subroutine order_k2a
  
  subroutine make_mca
    implicit none
    integer ikp, st
    allocate(mca(nk),stat=st); call check0(st,'mca')
    do ik=1,nk
       ratiok = (k2a(ik)/k2a(1))**0.5d0
       mca(ik) = ratiok**(-falloff) * nmc   ! note the minus; e.g., falloff=4--> mca ~ k^-4
       if(mca(ik)<mca_floor) mca(ik) = 0
    enddo
    
    ! superflouous, but checking the correct ordering never hurts
    call check(pntrk(1),1,' pntrk(1),one ')
    do i=1+1,nk
       ik =pntrk(i)
       ikp=pntrk(i-1)
       call check_le(mca(ik),mca(ikp),'mca_ik, mca_ikp ')
    end do

    ! now find how many k points would be monte carloed
    nksamp = 0
    do i=1,nk
       ik = pntrk(i)
       if(mca(ik)>0) nksamp=nksamp+1
    enddo
    if(rank==0) write(ip,*)' # fuzzy k points sampled     ',nksamp
    if(rank==0) write(ip,*)' # fuzzy mc sampling -- k=0   ',mca(1)
    if(rank==0) write(ip,*)' # fuzzy mc sampling -- total ',sum(mca)
    if(rank==0) call flush(ip)
  end subroutine make_mca

  subroutine plot
    implicit none
    integer, external :: is_it_power2
    integer, save     :: i1=1
    real*8 array(4)
    if(rank==0.and.is_it_power2(imc)==1) then
       if(i1==1) call pre_plot
       i1=-1

       savg = vk_exch(1)/dble(imc) 
       ds = sqrt(max(1d-10,s2/dble(imc)-savg**2))/sqrt(dble(imc))
       
       array(1) = vk_exch(pntrk(1))/min(imc,mca(pntrk(1)))
       array(2) = vk_exch(pntrk(2))/min(imc,mca(pntrk(2)))
       i=min(30,nksamp)
       array(3) = vk_exch(pntrk(i))/min(imc,mca(pntrk(i)))
       i=min(300,nksamp)
       array(4) = vk_exch(pntrk(i))/min(imc,mca(pntrk(i)))

       write(ip,*)'vkx(pntr(1,2,30,300)) d1 imc ',&
                   real(array), real(ds), imc
       call flush(ip)

    end if
  end subroutine plot

  subroutine pre_plot
    implicit none
    integer ii(4)
    ii=(/1,2,30,300/)
    where(ii>nksamp)ii=nksamp
    write(ip,*)
    write(ip,*)'         ii=(1,2,30,300) (limited by nsamp) ',ii
    write(ip,*)'   mca(pntrk(1,2,30,300)  (limited) ',mca(pntrk(ii(1:4)))
    write(ip,*)'   k2( pntrk(1,2,30,300)  (limited) ',real(k2a(pntrk(ii(1:4))))
    write(ip,*)' 4pi/k2(pntr(1,2,30,300)  (limited) ',real(4d0*pi/k2a(pntrk(ii(1:4))))
    write(ip,*)
    call flush(ip)
  end subroutine pre_plot

  subroutine nrmlz_vk_exch
    implicit none
    do i=1,nk
       ik=pntrk(i)
       if(mca(ik)>0)vk_exch(ik)=vk_exch(ik)/mca(ik)
    enddo
  end subroutine nrmlz_vk_exch

  subroutine dealloc
    implicit none
    integer st
    deallocate(k2a, pntrk, mca, stat=st); call check0(st,' deallocating k2a, etc. ')
  end subroutine dealloc
end subroutine vk_exch_mc_2d_fuzzy
