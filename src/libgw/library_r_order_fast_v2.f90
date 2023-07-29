!
! Copied with a few changes (to give the pointer ) from 
!    http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/sort3_f90.txt
!
! 
!***************************************************************
!*          Sorting an array with the Heapsort method          *
!* ----------------------------------------------------------- *
!* REFERENCE:                                                  *
!*      "NUMERICAL RECIPES By W.H. Press, B.P. Flannery,       *
!*       S.A. Teukolsky and W.T. Vetterling, Cambridge         *
!*       University Press, 1986" [BIBLI 08].                   *
!* ----------------------------------------------------------- *

!*                                                             *
!***************************************************************

!program sort3_driver_program
!  call sort3_driver
!end program sort3_driver_program

!*****************************************************
!*  Sorts a COPY of array RA of length N in ascending order *
!*                by the Heapsort method             *
!* AND ALSO GIVES THE POINTER TO THE OLD VALUES
!* ------------------------------------------------- *
!* INPUTS:                                           *
!     *    N  size of table w                        *
!     *          w  table -- unchanged!             *
!* OUTPUT:                                           *
!     *    i_old_ofnew    table sorted in ascending order    *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************         
! I label by expanding each number to a pair, (a(i),i)
!

SUBROUTINE order_r_index(w,i_old_ofnew,n)
  implicit none
  integer n, st
  integer i_old_ofnew(n)
  integer L,ir,i,k,j
  real*8  w(n)
  real*8  rra(2)
  real*8, allocatable :: a(:,:)
  allocate(a(2,n),stat=st); if(st/=0)stop ' a_order '

  do i=1,n
     a(1,i) = w(i)
     a(2,i) = i  
  enddo

  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
     L=L-1
     RRA =a(:,L)
  else
     RRA =a(:,ir)
     a(:,ir)=a(:,1)
     IR=IR-1
     if(IR.eq.1)then
        a(:,1)=rra
        i_old_ofnew(:)=nint(a(2,:))
        deallocate(a)
        call check_i_old_ofnew(w,i_old_ofnew,n) ! can be omitted
        return
     end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
     if(J < IR)then
        if(a(1,J) < a(1,J+1))  J=J+1
     end if
     if(rra(1) < a(1,J))then
        a(:,I)=a(:,J)
        I=J; J=J+J
     else
        J=IR+1
     end if
     goto 20
  end if
  a(:,I)=rra
  goto 10
END

subroutine SORT3_driver
  implicit none
  integer i_old_ofnew(100), n, i
  real*8    A(100)        !Table to be sorted
  real*8    MAX_VALUE, x  !Maximum value of table
  !initialize random number generator

  N=12  !initialize size of table
  if(n>100) stop ' n too big in sort3_driver '
  MAX_VALUE = 1000

  !generate random table of numbers (from 0 to 1000)
  do i=1, N
     x = abs(cos(i*34520.283988)) ! just a list from 0 to 1
     a(i)=MAX_VALUE*x
  end do

  if(n>8) then ! check equal values
     a(n/2+1)=a(n/2)
     a(n/2+2)=a(n/2)
     a(n/2+3)=a(n/2)
     a(n-1)=a(n)
  end if

  print  *,' '
  print *,'Table to be sorted:'
  call TWRIT(N,A)

  !call heapsort subroutine
  !call HPSORT(N,A)
  call order_r_index(a,i_old_ofnew,n)  ! does not reorder!!!

  print *,' '
  print *,'Sorted table (Heapsort method):'
  call TWRIT_n(N,A,i_old_ofnew)

  print *,' '
  stop

END subroutine SORT3_driver

!write table of size N to standard output
SUBROUTINE TWRIT(n,arr)
real*8 ARR(N)
  print *,' '
  WRITE(*,10) (ARR(I),I=1,N)
  return
 10    FORMAT(10F6.1)
END

!write sorted table of size N to standard output

SUBROUTINE TWRIT_n(n,arr,i_old_ofnew)
  implicit none
  integer n
  integer i_old_ofnew(n), i
  real*8  arr(n)

  print *,' array ordered'
  WRITE(*,10) (arr(i_old_ofnew(i)),I=1,N)
  print *,' indicies '
  write(*,11) (i_old_ofnew(i),i=1,n)
  return
 10    FORMAT(10F6.1)
 11    FORMAT(10i6)
END SUBROUTINE TWRIT_n

subroutine check_i_old_ofnew(w,i_old_ofnew,n)
  implicit none
  integer n, i_old_ofnew(n), st, i
  real*8  w(n), rp, r

  integer, allocatable :: occ(:)
  allocate(occ(n), stat=st)
  if(st/=0) stop ' ERROR: occ_in_check_i_old '
  occ = 0
  do i=1,n
     occ(i_old_ofnew(i))= occ(i_old_ofnew(i))+1
  enddo
  if(minval(occ)<1.or.maxval(occ)>1) then
     write(6,*)' ERROR: problem in check_i_old_ofnew '
     write(6,*)' minval, maxval(occ) ',minval(occ),maxval(occ)
     stop
  endif
  deallocate(occ)

  do i=1+1,n
     rp = w(i_old_ofnew(i-1))
     r  = w(i_old_ofnew(i))
     if(rp>r) then
        write(6,*)' Error, stopping in check_i_old_ofnew, i,rp,r ',i,rp,r
        stop
     endif
  enddo
end subroutine check_i_old_ofnew
