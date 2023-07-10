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

! if i not equal to j, stop and print achar.

      subroutine check(i,j,achar)
      implicit none
      integer i,j
      character(*) achar
      if(i.ne.j) then 
         write(6,*)' stopping since ',i,' n.e. ', j,' in ', achar 
         stop
      endif
      end

      subroutine check0(i,achar)
      implicit none
      integer i
      character(*) achar
      if(i.ne.0) then 
         write(6,*)' stopping since,',i,' n.e.0 in ', achar
         stop
      endif
      end


! if x not close to y, stop and print achar.  Depends on tolerance.
      subroutine check_r(x,y,achar)
      implicit none
      real*8 x, y
      character(*) achar
      if(abs(x-y)/max(abs(x),abs(y),1.d-20)>1.d-12) then 
         write(6,*)x, y, achar 
         stop
      endif
      end


! if i not less or eqaul to j, stop and print achar.



      subroutine check_LE(i,j,achar)
      implicit none
      integer i,j
      character(*) achar
      if(i > j) then 
         write(6,*)' problem in check_le ',i, j, achar 
         stop
      endif
      end


      subroutine check_lele(i,j,k,achar)
      implicit none
      integer i,j,k
      character(*) achar
      if(i > j .or. j>k) then 
         write(6,*)' problem in check_lele ',i, j, k,achar 
         stop
      endif
      end

      subroutine check_real(cvec, n)
      implicit none
      integer n
      complex*16 cvec(n)
      real*8, parameter :: thresh_abs = 1.d-10
      real*8, parameter :: thresh_rel = 1.d-8
      real*8 normi, norma
      
      norma= sum(abs(cvec))
      normi= sum(abs(aimag(cvec)))
      
      if(normi>norma*thresh_rel.and.normi>thresh_abs) then
         write(6,*)' in check_real, normi = ',normi,' norma = ',norma
         stop
      endif
      end subroutine check_real

!-----------------------------------
! checks if idet=1,2,4,8,16,etc. (1 then; -1 otherwise)
!-----------------------------------
      function is_it_power2(i_in) 
      implicit none
      integer is_it_power2,i,i_in
      
      i = abs(i_in)
      if(i==0.or.i==1) then
         is_it_power2 = 1
         return
      endif

      is_it_power2= -1
      if(  2**nint((dlog(dble(i))/dlog(2.d0)))==i) 
     &                                     is_it_power2 = 1
      end function is_it_power2

      subroutine check_real_LE(x,y,achar)
      implicit none
      real*8 x,y
      character(*) achar
      if(x > y) then 
         write(6,*)x, y, achar 
         stop
      endif
      end

      function is_i_near_power_x(ii,x)
      implicit none
      integer ii, is_i_near_power_x, j, ans
      real*8 x,y
      
      j = int(log(dble(ii))/log(x))
      y = x**dble(j)
      
      if(y.le.ii.and.y+1.gt.ii) then
         ans = 1
      else
         ans = 0
      end if
      
      is_i_near_power_x = ans
      end function is_i_near_power_x

c kronecker delta; reutrns REAL*8 results.
      function delta(i,j)
      implicit none
      integer i,j
      real*8 delta
      if(i.eq.j) then
         delta = 1.d0
      else
         delta = 0.d0
      endif
      end


      subroutine check_lelele(i,j,k,l,achar)
      implicit none
      integer i,j,k,l
      character(*) achar
      if(i > j .or. j>k .or. k>l) then 
         write(6,*)' problem in check_lelele ',i, j, k,l,achar 
         stop
      endif
      end

      subroutine check_r_le(x,y,achar)
      implicit none
      real*8 x,y
      character(*) achar
      if(x > y) then 
         write(6,*)x, y, achar 
         stop
      endif
      end


      subroutine check_r_lele(x,y,z,achar)
      implicit none
      real*8 x,y,z
      character(*) achar
      if(x > y.or.y>z) then 
         write(6,*)x, y, z, achar 
         stop
      endif
      end



