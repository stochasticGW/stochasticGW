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
subroutine assign_flg_units_nuc
  use gwm, only : units_nuc,flg_units_nuc
  implicit none
  character(len=3):: c3
  integer i
  c3=units_nuc(1:3)

  i=0
  call lower_the_case(c3)
  select case(c3)
  case('boh','au')
     i=1
  case('ang')
     i=2
  case default
     write(6,*)" Error: Units of Nuclear Positions should be Angstrom or Bohr "
     write(6,*)"   while the input was ",units_nuc
     stop
  end select
  flg_units_nuc = i
end subroutine assign_flg_units_nuc

subroutine cp_character_to_3(a,c3)
  implicit none
  character a(*), c3(3)
  c3 = a(1:3)
end subroutine cp_character_to_3
