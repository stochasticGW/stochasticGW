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
module lda_modu
   implicit none
  save
  integer,parameter   :: ma=501000
  real*8, parameter   :: rhomin = 1d-12, rhomax=1d8
  logical             :: first_lda_exc=.true., flushlda=.false.
  real*8, allocatable :: logrhoa(:), vcna(:), elca(:)
  real*8              :: rm, logrm, logrmin
contains
  subroutine prep_lda_exc
    use simple_mpi, only : rank
    implicit none
    integer i
    real*8, external ::   elc, n_delc_dn
    real*8           ::   r
    first_lda_exc = .false.
    call check_real_le(rhomin, rhomax, ' rmin, rmax    ')
    rm = (rhomax/rhomin)**(1d0/(ma-1d0))
    logrm = log(rm)
    logrmin  = log(rhomin)
    
    call check_le(2,ma,' two, ma  ')
    allocate(logrhoa(ma),stat=i); if(i/=0) stop ' vcn logrhoa problems '
    allocate( vcna(ma),stat=i); if(i/=0) stop ' vcn vvcna problems '
    allocate( elca(ma),stat=i); if(i/=0) stop ' vcn vvcna problems '
    
    if(rank==0.and.flushlda) open(50512,file='lda_exc.txt',status='replace')
    do i=1,ma
       logrhoa(i)   = logrmin + (i-1)*logrm
       r          = exp(logrhoa(i))
       elca(i)    =  elc(r) 
       vcna(i)    =  elca(i) + n_delc_dn(r)
       if(rank==0.and.flushlda) write(50512,"(' ',4e18.8)") r,logrhoa(i),vcna(i),elca(i)
    enddo
    if(rank==0.and.flushlda)call flush(50512) 
  end subroutine prep_lda_exc
end module lda_modu

subroutine get_vx_lda(nn, ng, vx)
  implicit none

  integer ng
  real*8 nn(ng), vx(ng)

  vx = -0.738558766382022046 * 4.d0/3.d0* (max(nn,1e-20))**(1.d0/3.d0)
end subroutine get_vx_lda

subroutine get_vx_spn_lda(nn, ng, nspn, vx)
  implicit none
  integer ng, st, nspn, is
  real*8  nn(ng, nspn), vx(ng, nspn)
  real*8, allocatable :: nn1(:)

  select case(nspn)
  case(1)
     call get_vx_lda(nn, ng, vx)
  case(2)
     allocate(nn1(ng), stat=st)
     call check0(st,' nn1 ')
     do is=1,nspn
        nn1 = 2d0* nn(:,is)
        call get_vx_lda(nn1, ng, vx(1,is))
     enddo
     deallocate(nn1)
  case default
     stop ' problem with nspn '
  end select
  
end subroutine get_vx_spn_lda

subroutine get_vcn_spn_lda(nn, ng, nspn, vc)
  implicit none
  integer ng, nspn, ig, is
  real*8 nn(ng, nspn), vc(ng, nspn), n1, n2, dn, nnn, ec0a, de1, de2
  real*8, external :: ec_pw

  select case(nspn)
  case(1)
     do ig=1,ng
        n1 = max(1e-20,nn(ig,1)/2d0 )
        n2 = n1
        nnn = n1+n2
        ec0a  = ec_pw(n1,n2)
        dn = 1e-5 * nnn
        !de1 = (ec_pw(n1+dn,n2   )-ec0a)/dn
        de1 = (ec_pw(n1+dn/2d0,n2+dn/2d0)-ec0a)/dn
        vc(ig,1) = ec0a + nnn*de1
     enddo
  case(2)
     do ig=1,ng
        !
        ! vc(ig,1) = d( e(n1,n2)*(nnn) )/ dn1 = e(n1,n2) + nnn*de(n1,n2)/dn1
        !
        n1 = max(1e-20,nn(ig,1))
        n2 = max(1e-20,nn(ig,2))
        nnn = n1+n2
        ec0a  = ec_pw(n1,n2)
        dn = 1e-5 * nnn
        de1 = (ec_pw(n1+dn,n2   )-ec0a)/dn
        de2 = (ec_pw(n1,   n2+dn)-ec0a)/dn
        vc(ig,1) = ec0a + nnn*de1
        vc(ig,2) = ec0a + nnn*de2
     end do
  case default
     stop ' nspn in vcn_spn '
  end select
end subroutine get_vcn_spn_lda

subroutine get_vcn_lda(nn , ng, vc)
  implicit none
  integer ng
  real*8 vc(ng), nn(ng),  nl, vcnf
  integer ig

  do ig=1, ng
     nl = nn(ig)
     vc(ig) = vcnf(nl) 
  end do
end subroutine get_vcn_lda

function vcnf(nl)
  use lda_modu, only : rhomin,rhomax, logrhoa, ma,  vcna, first_lda_exc, logrm, prep_lda_exc
  implicit none
  integer i
  real*8 nl, vcnf, a, x, lognl

  if(first_lda_exc) then
     call prep_lda_exc
  endif

  if(nl<rhomin) then
     vcnf=vcna(1)
     return
  endif

  lognl = log(nl)
  a = (lognl-logrhoa(1))/logrm 
  i = int(a)
  x = a-i

  if(i<1) then
     vcnf = vcna(1)
  else if(i>ma-1) then
     vcnf = vcna(ma)
  else 
     vcnf = x*vcna(i+1)+(1d0-x)*vcna(i)
  endif     
end function vcnf

function n_delc_dn( ntot)
  implicit none
  real*8 n_delc_dn, ntot, elc
  
  real*8, parameter :: ddn = 1d-6

  if(ntot<1.d-10) then
     n_delc_dn = 0d0
     return
  end if
  
  n_delc_dn = &
              (elc( ntot *(1+ddn)) -   &
               elc( ntot *(1-ddn)))   &
               / 2.d0 / ddn

end function n_delc_dn
  
!
! checks elc, vs. PW 1991 (Perdew Wang 1991, PRB 13298), Table III
!
subroutine check_elc
  implicit none
  real*8 rs, pi, ntot, term, elc
  
99 write(6,*)' rs '
  read( 5,*) rs
  if(rs<0) return
  
  pi = dacos(-1.d0)
  ntot = 1.d0/ ( 4*pi*rs**3/3.d0)                      ! Jelium density
  
  term = elc(ntot)

  write(6,*)' rs, ntot ',real(rs),real(ntot),' elc(mRY) ',real(term)*(-2000)
  goto 99
end subroutine check_elc

  
function elc( ntot) 
!
! Local (spn-unpolarized) function 
!
  implicit none
  real*8 elc, ntot, ec_pw
 
  if(ntot<1.d-10) then
     elc = 0d0
     return
  end if
  elc = ec_pw(ntot/2.d0, ntot/2.d0)    ! contributions of up+ cont.of dn

end function elc

!
! Perdeow and Wang, PRB 45, p. 13244, 1992, eqs. 1,2,10
!

function ec_pw( nup, ndn)
  implicit none

  real*8 nup, ndn, rs, z, pi
  real*8 ec0, ec1, ac, fz, ec_pw, gf

  real*8, parameter :: fdd =   1.709921; 

  !
  ! Table I, pw, last 3 columns
  !
  real*8, parameter ::   &
       A0 = 0.031091,  A1 = 0.015545,  A2 = 0.016887,  &
       a10 = 0.21370,  a11 = 0.20548,  a12 = 0.11125,  &
       b10 = 7.5957,   b11 =14.1189,   b12 =10.357,    &
       b20 = 3.5876,   b21 = 6.1977,   b22 = 3.6231,   &
       b30 = 1.6382,   b31 = 3.3662,   b32 = 0.88026,  &
       b40 = 0.49294,  b41 = 0.62517,  b42 = 0.49671

  pi = dacos(-1.d0)
  
  rs = (3.d0/(4.d0*pi*(nup+ndn)))**(1.d0/3.d0)
  
  z  = (nup-ndn)/(nup+ndn)
  
  !
  !
  ! eq. 10, pw  (p=1)
  !

  ec0 = Gf(rs, A0, a10, b10, b20, b30, b40)
  ec1 = Gf(rs, A1, a11, b11, b21, b31, b41)
  ac  = Gf(rs, A2, a12, b12, b22, b32, b42)
  
  !
  ! eq. 9,8, pw
  !

  fz  = ((1d0+z)**(4.d0/3.d0) + (1d0-z)**(4.d0/3.d0) -2d0 )/  &
       ( 2**(4.d0/3.d0)-2)                           

!  ec_pw  = ec0 + ac * fz / fdd *(1-z**4) + (ec1-ec0)*fz*z**4   ! corrected
   ec_pw  = ec0 - ac * fz / fdd *(1-z**4) + (ec1-ec0)*fz*z**4   

end function ec_pw
   
function Gf(rs, A, a1, b1, b2, b3, b4)   !perdew-wang 1992, eq. 10, p=1
  implicit none

  real*8 Gf, rs, A, a1, b1, b2, b3, b4, term


  term = 2*A*(b1* rs**0.5d0 + b2 * rs + b3 * rs **1.5d0 + b4 * rs**2)

  Gf  = -2d0*A*(1d0 + a1 * rs)* log(1d0 + 1.d0/term)
end function Gf

