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
subroutine set_theta
  use simple_mpi, only : rank, sync_mpi
  use gwm, only :   dh, havg, tp, toll,mu, nchb=>nchbmx, nchbc, co=>th_co
  use gwm, only :   gapped_flg, homo, lumo
  implicit none
  integer st,lc
  allocate(co(0:nchb),stat=st); if(st/=0)stop ' co '
 
  if(rank==0) then
     write(17,*)' nchb,tp ',nchb,tp
     write(17,*)' mu ',mu
     write(17,*)' havg dh ',havg,dh
     write(17,*)' toll ',toll
     write(17,*)' size_co   ',size(co)
     if (gapped_flg) then
        write(17,*)' scaled_gap',(lumo-homo)/dh
        write(6,'(X,A,F16.8)') ' scaled_gap = ',(lumo-homo)/dh
     endif
  end if

  if (gapped_flg) then
    if(rank==0) write(6,'(X,A,/4X,A/)') ' -> Gapped-filtering will be used',&
                           ' to prepare the occupied stochastic states'
    call cheb_coeff_theta_zero_gap_new(nchb, havg, DH, co, lc, homo, lumo)
    nchbc = lc+1  ! 0:lc == 0:nchbc-1
  else
    if(rank==0) write(6,'(X,A,/4X,A/)') ' -> Temperature-dependent Heaviside filtering will be used',&
                           ' to prepare the occupied stochastic states'
    call cheb_coeff_theta_power_r(nchb, Tp, mu, havg, DH, toll, co, lc, 0.5d0)
    nchbc = lc+1  ! 0:lc == 0:nchbc-1
  endif

  if(rank==0) then;
     write(17,*)' maximum chebyshev index: 0 to nchbc-1, where nchbc =',nchbc
     write(17,*)' several values of co'
     write(17,*)' i=0:3 ',co(0:3)
     write(17,*)' all value',co
     write(17,*)' max,min value',maxval(co), minval(co)
     write(17,*)' i=nchbc/2, nchbc*0.8, nchbc-2,nchbc-1 ',&
                co(nchbc/2), co((nchbc*8)/10),co(nchbc-2:nchbc-1)
     call flush(17)
  end if

  call plotfilterdata

contains
  subroutine plotfilterdata

     implicit none
     integer :: i

     if (rank==0) then

        open(104, file='filtercheb.dat')

        write(104,*) 'havg dh:'
        write(104,*)  havg,dh
        if (gapped_flg) then
           write(104,*) 'homo lumo'
           write(104,*)  homo,lumo
        else
           write(104,*) 'mu'
           write(104,*)  mu
        endif
        write(104,*) 'coefs:'
        do i=0,nchb
           write(104,*) co(i)
        enddo

        close(104)

     endif

   end subroutine plotfilterdata

end subroutine set_theta
