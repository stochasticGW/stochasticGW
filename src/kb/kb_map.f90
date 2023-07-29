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
subroutine kb_map
  use kb_mod
  use simple_mpi, only : sync_mpi
  implicit none
  integer, save :: i1=1
  integer       :: ia, ma, nga, ig, rep
  integer       :: ix0,iy0,iz0
  integer       :: ix, iy, iz
  integer       :: ixp, iyp, izp
  integer       :: mmx, mmy, mmz
  real*8        :: d2,r0(3),rg(3),rn(3)
  logical       :: found
  
  if(i1/=1) return
  i1=-1

  r0 = (/ -dble(nx)/2d0*dx, -dble(ny)/2d0*dy, -dble(nz)/2d0*dz /)

  if(dim_periodic.gt.0 .and.2d0*maxval(racut_a)>min(nx*dx, ny*dy, nz*dz)) then
     if(rnk==0) then
        write(6,*)' ERROR: '
        write(6,*)' Periodic, & cell too small, images of nuclei in NL overlap '
        write(6,*)' Thus, cell size ',real(nx*dx),real(ny*dy),real(nz*dz)
        write(6,*)'  while max cutoff  = ', real(maxval(racut_a))
        stop
     endif
  endif
  call sync_mpi
  call check_cell()
  do ma=1,size(racut_a)
     mmx = racut_a(ma)/dx + 2;
     mmy = racut_a(ma)/dy + 2;
     mmz = racut_a(ma)/dz + 2;
     if(rnk==0)then
        write(17,*)' ma,ch,mmxyz r, cushion for  ',ma,mapat(ma),mmx,mmy,mmz,real(racut_a(ma))
     endif
  end do
  if(rnk==0) call flush(17)

  reploop : do rep = 1,2
     ialoop : do ia=1,na
        ma = mapai(ia)
        mmx = racut_a(ma)/dx + 2;
        mmy = racut_a(ma)/dy + 2;
        mmz = racut_a(ma)/dz + 2;
        nga = 0

        rn = cnt(:,ia)
               
        if(nx>1) then; ix0 = (rn(1) - r0(1))/dx; else ; ix0=1; endif
        if(ny>1) then; iy0 = (rn(2) - r0(2))/dy; else ; iy0=1; endif
        if(nz>1) then; iz0 = (rn(3) - r0(3))/dz; else ; iz0=1; endif
        izloop : do izp=iz0-mmz,iz0+mmz
           iyloop : do iyp=iy0-mmy,iy0+mmy
              ixloop: do ixp=ix0-mmx,ix0+mmx

                 found = .false.
                 rg = r0 + (/ (ixp-1)*dx, (iyp-1)*dy, (izp-1)*dz /) 
                 d2 = sum((rg-rn)**2)
   
                 if(dim_periodic.gt.0) then
                    ix = mod(ixp-1+nx,nx)+1; call check_lele(1,ix,nx,'1-ix-nx')
                    iy = mod(iyp-1+ny,ny)+1; call check_lele(1,iy,ny,'1-iy-ny')
                    iz = mod(izp-1+nz,nz)+1; call check_lele(1,iz,nz,'1-iz-nz')
                 else
                    ix=ixp
                    iy=iyp
                    iz=izp
                    if(ix<1.or.ix>nx) cycle ixloop
                    if(iy<1.or.iy>ny) cycle iyloop
                    if(iz<1.or.iz>nz) cycle izloop
                 end if

                 ig = ix+(iy-1)*nx+(iz-1)*nx*ny

                 nearby: if(d2<racut_a(ma)**2) then
                    if(found)stop ' ERROR: 2 atomic images overlapping in NL '
                    found=.true.
                    nga = nga+1
                    if(rep==2) then
                       mapkbg(  nga,ia) = ig
                       rgn(:,nga,ia) = rg(:)-rn(:) 
                    end if
                 end if nearby
              end do ixloop
           end do iyloop
        end do izloop
        ngs(ia)=nga
     end do ialoop
     
     rep1 : if(rep==1) then
        ngsmall = maxval(ngs(1:na))         
        if(rnk==0) write(17,*)' minval ngs, maxval(=ngsmall) ',minval(ngs(1:na)),ngsmall; 
        allocate(mapkbg(    ngsmall, na), stat=st); call check(st,0,' map_kb_g ' )
        allocate(rgn(     3,ngsmall, na), stat=st); call check(st,0,' rgn      ' )
        mapkbg = -999999
     end if rep1

  end do reploop

  call flush(17)
contains
    subroutine check_cell()
    use simple_mpi, only : sync_mpi, rank
    implicit none
    integer i
    real*8  rL(3),amn(3),amx(3),fctr
    rL = (/ nx*dx, ny*dy, nz*dz /)
    amn = (/ minval(cnt(1,:)),minval(cnt(2,:)),minval(cnt(3,:)) /)
    amx = (/ maxval(cnt(1,:)),maxval(cnt(2,:)),maxval(cnt(3,:)) /)
    
    rnk0: if(rank==0) then
       write(17,*)' Checking atoms vs. cell '
       write(17,*)' atoms_minmax, cell length '
       do i=1,3
          write(17,*)' direction: ',i,' atoms_minmax ',real(amn(i)),real(amx(i)),&
                                      ' cell-length ',real(rl(i))
       enddo
       prdc: if(dim_periodic.gt.0) then
          write(17,*)' periodic,     atoms cell should be within -fctr*cell,fctr*cell, fctr=1.5 ' 
          fctr = 1.5d0
       else
          write(17,*)' non-periodic, atoms cell should be within -fctr*cell,fctr*cell, fctr=0.5 ' 
          fctr = 0.5d0
       end if prdc
       do i=1,3
          call check_r_le(-fctr*rL(i),amn(i), ' -fctr*rL, minval_cnt ')
          call check_r_le(amx(i), fctr*rL(i), '  maxval_cnt, fctr*rL ')
       enddo
    end if rnk0
    call sync_mpi
  end subroutine check_cell
end subroutine kb_map
