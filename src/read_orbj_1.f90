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
subroutine read_orbj_1
  use gwm, only : rdorbj, n, dv, orb_kind, nj
  use simple_mpi, only : rank,bcast_r8
  implicit none
  logical                 :: fex
  character*30            :: a, ch
  integer                 :: st, ni, ip, ij, isp
  integer                 :: nx, ny, nz, dumin
  real*8                  :: dx, dy, dz, dvrd
  complex*16, allocatable :: uc(:)

  ip = 111

  if(.not.allocated(rdorbj)) then; 
     allocate(rdorbj(n,nj),stat=st) 
     if(st/=0) stop ' rdorbj problems '
  endif

  rnk0read : if(rank==0) then
     write(6,*) " ############# READING DFT ORBITAL FROM orbj.txt #############"
     if(n/=size(rdorbj,1)) then
        write(6,*) " grid sizes (should be identical!):",n, size(rdorbj,1)
        call flush(6)
        stop
     end if

     if(nj/=size(rdorbj,2)) then
        write(6,*) " number of off diagonal orbitals :",nj, size(rdorbj,2)
        call flush(6)
        stop
     end if

     inquire(file='orbj.txt',exist=fex)
     if(.not.fex) stop 'orbj.txt missing'

     open(ip,file='orbj.txt',status='old')
     rewind(ip)
     a=' from file: orbj.txt '

     read(ip,*)ch, nx; if(ch/='nx' )then;write(6,*)' ERROR, expect nx,got ',ch,a;stop;endif
     read(ip,*)ch, ny; if(ch/='ny' )then;write(6,*)' ERROR, expect ny,got ',ch,a;stop;endif
     read(ip,*)ch, nz; if(ch/='nz' )then;write(6,*)' ERROR, expect nz,got ',ch,a;stop;endif
     read(ip,*)ch, dx; if(ch/='dx' )then;write(6,*)' ERROR, expect dx,got ',ch,a;stop;endif
     read(ip,*)ch, dy; if(ch/='dy' )then;write(6,*)' ERROR, expect dy,got ',ch,a;stop;endif
     read(ip,*)ch, dz; if(ch/='dz' )then;write(6,*)' ERROR, expect dz,got ',ch,a;stop;endif
     read(ip,*)ch, isp;if(ch/='isp')then;write(6,*)' ERROR, expect isp,got ',ch,a;stop;endif
     read(ip,*)ch;   if(ch/='orbj')then;write(6,*)' ERROR, expect orbj,got ',ch,a;stop;endif

     ni = nx*ny*nz
     dvrd = dx*dy*dz
     
     if(n.ne.ni) then
        write(6,*)' n=',n,' not equal to ni=n_read= ',ni,' in orbj.txt '
        stop
     end if
     if(dvrd /= dv) then
        write(6,*) " ERROR: volume elements not equal:"
        write(6,*) " from orbj.txt",dvrd," but from dens.txt:", dv
        call flush(6)
        stop
     end if

     select case(orb_kind)
     case(1)
        do ij=1,nj 
          read(ip,*) rdorbj(:,ij)
        enddo
     case(2)
        allocate(uc(n),stat=st);call check0(st,' uc ') 
        do ij=1,nj
          read(ip,*) uc
          rdorbj(:,ij) = dble(uc)
        enddo
        deallocate(uc)
     case default
        write(6,*)' ERROR orb_kind not 1(=real*8) nor 2(=complex*16) but ',orb_kind
        call flush(6)
        stop
     end select

     write(17,*) 
     do ij=1,nj
       write(17,*) "-> sum(psi^2)",sum(rdorbj(:,ij)**2)*dv
       if(abs(sum(rdorbj(:,ij)**2)*dv -1.d0)>1d-7) then
          write(6,'(X,A,I0,A)') " WARNING: orbital ",ij," in orbj.txt orbital is not correct -- as it isnt normalized !"
          write(6,*) " Specifcially, its integral =",sum(rdorbj(:,ij)**2)*dv
          write(6,*) " The calculation is continued but... Good luck!" 
          call flush(6)
!         stop
       end if
     enddo
     call flush(17)

     close(ip)
  end if rnk0read
  call bcast_r8(rdorbj,size(rdorbj),0)
end subroutine read_orbj_1
