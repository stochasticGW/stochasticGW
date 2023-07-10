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
!program tmp
!  call make_atoms_list
!end program tmp

subroutine make_atoms_list
  use atoms, only : an => atom_name, un=> atom_upper_name, n118
  implicit none
  integer i
  call check(size(an),size(un),' an-vs-un ')
  call check(size(an),n118,    ' an-vs-118')

  call atoms_data
  do i=1,size(an)
     call make_upper(an(i),un(i))
  enddo
  !call print_atoms
  call check_atoms
contains

  subroutine print_atoms
    implicit none
    write(6,*)' Element and Upper Case '
    do i=1,size(an)
       write(6,*) i,an(i),',',un(i)
    enddo
  end subroutine print_atoms

  subroutine check_atoms
    implicit none
    integer j
    do i=1,size(an)
       do j=1,size(an)
          if(i/=j) then
             if(an(i)==an(j).or.an(i)==un(j).or.un(i)==un(j)) then
                write(6,*)' problem in nanes, i,j, names: ',&
                        i,j,an(i),',',an(j),',',un(i),',',un(j)
                stop
             endif
          end if
       end do
    end do
  end subroutine check_atoms
end subroutine make_atoms_list

subroutine atoms_data
  use atoms, only : an =>atom_name
  implicit none
  an( 1: 2)=&
    (/'H ',                                                                                'He'/)
  an( 3:10)=& 
    (/'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne'/)
  an(11:18)=&
    (/'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar'/)
  an(19:36)=&
    (/'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr'/)
  an(37:54)=&
    (/'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe'/)
  an(55:56)=&
    (/'Cs','Ba'/)
  an(72:86)=&
                   (/'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn'/)
  an(87:88)=&
    (/'Fr','Ra'/)
  an(104:118)=     (/'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/)

  an(57:71)=       (/'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'/)
  an(89:103)=      (/'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/)

end subroutine atoms_data
  
subroutine make_upper(nm,unm)
  implicit none
  character :: nm(2), unm(2)
  character :: lc, uc
  integer   :: ila, iua, ilz, iuz, ic, iu

  lc=nm(2)
  if(lc==' ') then
     unm(1)=nm(1)
     unm(2)=nm(2)
     return
  endif

  ila=IACHAR('a') 
  iua=IACHAR('A') 
  ilz=IACHAR('z') 
  iuz=IACHAR('Z') 
  ic=IACHAR(lc) 
  call check_lele(ila,ic,ilz,' ila, ic, ilz ')

  uc=ACHAR(ic + (iua-ila)) 
  if(IACHAR(uc)-iua.ne.ic-ila) stop ' ERROR: problem in converting to uppercase '

  iu=IACHAR(uc)
  call check_lele(iua,iu,iuz,' iuc, iu, iuz ')
  
  unm(1)=nm(1)
  unm(2)=uc
end subroutine make_upper

