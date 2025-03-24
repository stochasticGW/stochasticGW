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
subroutine cheb_coeff_theta_power_r(nchb, Temperature, mu, havg, DL,toll, co, nchbtop, power)
  ! chebyshev coeff. for representating f(x)= theta_smooth(x) = erfc(-x/Temperature)/2.d0 as
  ! f(x)=sum_k=0,1,2...nchb  co(k) T_n( x/DL), for -DL<x<DL.
  ! nchbtop is  the maximum index subject to a tolerance on the values of toll.
  use mpi_lib_ours, only : rank
  implicit none

  integer nchb, ip, nchbtop ! nchbtop: output
  real*8  temperature, mu, havg, dl, toll, co(0:nchb), power ! co is output

  integer i,j, m,N
  real*8  pi,w,D, x, dx
  real*8  mu_s
  complex*16, allocatable :: f(:)
  real*8, allocatable :: fapprox(:), tn(:,:)

  pi = dacos(-1.d0)

! N  is about 4 or 8 times biger than nchb, and is 2**m
  m = 3 + nint( log(dble(nchb))/log(2.d0))
  N = 2**m  
  
  allocate(f(0:N-1), stat=i); call check(i,0,' falloc  ')
  allocate(fapprox(0:N-1), stat=i); call check(i,0,' fapprox alloc  ')
  allocate(tn(0:N-1,0:nchb), stat=i); call check(i,0,' tn alloc  ')

  do i=0,N-1
     w = 2d0*pi/dble(N) * i
     D = DL * cos(w)+havg
     f(i) = (erfc((D-mu)/Temperature)/2.d0)**power
     !write(103, *) f(i)
  enddo

  call fftsa(N,f,m)
  f = 2.d0*f/dble(N)
  f(0) = f(0)/2.d0
  co   = f(0:nchb)

  nchbtop = nchb 
  deallocate(f,stat=i); call check(i,0,' fdealloc  ')
end subroutine cheb_coeff_theta_power_r

subroutine cheb_coeff_theta_zero_gap(nchb, Temperature, mu, havg, DL, co0, nchbtop, homo, lumo)
  ! chebyshev coeff. for representating f(x)= theta_smooth(x) = erfc(-x/Temperature)/2.d0 as
  ! f(x)=sum_k=0,1,2...nchb  co(k) T_n( x/DL), for -DL<x<DL.
  ! nchbtop is the maximum index subject to a tolerance on the values of toll.
  use simple_mpi, only : rank
  use mat_module
  implicit none

  integer nchb, ip, nchbtop ! nchbtop: output
  real*8  temperature, mu, havg, dl, co0(0:nchb), homo, lumo! co is output

  integer i, j, m, N
  real*8  pi, th, x, D, dx, mu_s, homo_s, lumo_s, alpha, c
  real*8  gap, shift, factor
  real*8, allocatable :: fm(:), mat(:,:), minv(:,:), b(:), co(:), tn(:,:), vec(:,:), eig(:)
  complex*16, allocatable :: w(:), f0(:)

  factor = 5d-2
  alpha = 1d-6
  c = 1d-6
  shift = (lumo - homo)*factor
  pi = dacos(-1.d0)
  mu_s = (mu-havg)/dl
  homo_s = ((homo+shift)-havg)/dl 
  lumo_s = ((lumo-shift)-havg)/dl

! N  is about 4 or 8 times biger than nchb, and is 2**m

  m = 6 + nint( log(dble(nchb))/log(2.d0))
  N = 2**m  
  
  allocate(f0(0:N), stat=i); call check(i,0,' f0 alloc  ')
  allocate(w(0:N), stat=i); call check(i,0,' w alloc  ')
  allocate(fm(0:2*nchb), stat=i); call check(i,0,' fm alloc  ')
  allocate(mat(0:nchb,0:nchb), stat=i); call check(i,0,' mat alloc  ')
  allocate(minv(0:nchb,0:nchb), stat=i); call check(i,0,' mat_inv alloc  ')
  allocate(b(0:nchb), stat=i); call check(i,0,' b alloc  ')
  allocate(vec(0:nchb,0:nchb), eig(0:nchb), stat=i); call check(i,0,' vec alloc  ')

  do i=0,N-1
    th = 2d0*pi/dble(N) * i
    if ((th .ge. acos(lumo_s)) .and. (th .le. acos(homo_s))) then
      w(i) = c !(1-x^2)^{-1} implicitly in integral 
    else if ((th .ge. 2d0*pi-acos(homo_s)) .and. (th .le. 2d0*pi-acos(lumo_s))) then
      w(i) = c 
    else
      w(i) = 1d0
    endif
    D = DL * cos(th)+havg
    f0(i) =  (erfc((D-mu)*1E12)/2d0)*w(i)
  enddo
  close(203)

  call fftsa(N,w,m) ! w(n) = N/2pi int_0^2pi e^(-i*n x) w(x) dx

  w = w/dble(N)
  fm = pi*w(0:2*nchb) !Note we are taking the real part i.e cos(nx)

  call fftsa(N,f0,m)
  b = f0(0:nchb)/dble(N)*pi

  do i=0,nchb
    do j=0,nchb
      if (i .eq. j) then
        mat(i,j) = 0.5d0*(fm(abs(i-j)) + fm(i+j)) + alpha
      else
        mat(i,j) = 0.5d0*(fm(abs(i-j)) + fm(i+j)) 
      endif
    enddo
  enddo
  if (any(mat .ne. mat)) stop ' mat has nan '

  call mat_inv(mat,minv)
  co0 = matmul(minv, b)
  call mat_diag(mat,vec,eig)

  nchbtop = nchb

end subroutine cheb_coeff_theta_zero_gap

subroutine cheb_coeff_theta_zero_gap_new(nchb, havg, DL, co0, nchbtop, homo, lumo)
  use simple_mpi, only : rank
  use mat_module, only : mat_diag

  implicit none
  integer nchb, nchbtop ! nchbtop: output
  real*8  havg, dl, co0(0:nchb), homo, lumo! co is output
  integer i, j, N
  real*8  pi, homo_s, lumo_s, toll
  real*8  shift, factor, theta_l, theta_h
  real*8, allocatable :: mat(:,:), minv(:,:), b(:), vec(:,:), eig(:)

  factor = 1d-2
  homo_s = (homo-havg)/dl
  lumo_s = (lumo-havg)/dl
  shift = (lumo_s - homo_s)*factor
  pi = dacos(-1.d0)
  homo_s = homo_s + shift/2d0
  lumo_s = lumo_s - shift/2d0
  toll = 1d-12

  theta_h = acos(homo_s)
  theta_l = acos(lumo_s)

  allocate(mat(0:nchb,0:nchb), stat=i); call check(i,0,' mat alloc  ')
  allocate(b(0:nchb), stat=i); call check(i,0,' b alloc  ')
  allocate(minv(0:nchb,0:nchb), stat=i); call check(i,0,' mat_inv alloc  ')
  allocate(vec(0:nchb,0:nchb), eig(0:nchb), stat=i); call check(i,0,' vec alloc  ')

  do i=0,nchb
    do j=0,nchb
      mat(i,j) = f(pi,theta_h)+f(theta_l,0d0)
    enddo
  enddo

  b(0) = pi - theta_h
  do i=1,nchb
    b(i) = -sin(i*theta_h)/i
  enddo

  call mat_diag(mat,vec,eig)

  if(rank==0) then
    write(17,*)
    write(17,*)'eig(0:4)  ',real(eig(0:4))
    write(17,*)'eig(nchb-4:nchb)',real(eig(nchb-4:nchb))
    if(minval(eig).le.0) write(17,*) 'negative eigvalue present'
    write(17,*)' maxval(eig)/minval(eig)', maxval(eig)/minval(eig)
  endif

  do i=0,nchb
    if(eig(i)>toll) then
      vec(:,i) = vec(:,i)/sqrt(eig(i))
    else
      vec(:,i) = 0d0
    endif
  end do
  minv = matmul(vec,transpose(vec))
  co0 = matmul(minv,b)
  nchbtop = nchb

  contains
    real*8  function f(ty,tz)
      real*8 ty,tz
      f = g(ty)-g(tz)
    end function f

    real*8 function g(t)
      real*8 t
      if(i==j.and.i==0) g = t
      if(i==j.and.i>0)  g = t/2 + sin(2*i*t)/(4*i)
      if(i.ne.j)        g = sin((i-j)*t)/(2*(i-j)) + sin((i+j)*t)/(2*(i+j))
    end function g

end subroutine cheb_coeff_theta_zero_gap_new
