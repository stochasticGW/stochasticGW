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
module rand_seed_mod
  implicit none; save
  integer, parameter :: nlseed=4,ndiffseeds=6
  integer            :: line_seed
  integer            :: mdiffseeds
  integer            :: rank_dependent(    ndiffseeds) 
  integer            :: iu
  integer(kind=8)    :: j1, j2, j3, j4, j5
  integer(kind=8)    :: seed_array(nlseed, ndiffseeds)
  integer(kind=8)    :: seed_read( nlseed, ndiffseeds)
  logical            :: fex
  logical            :: first=.true. 
contains
  subroutine read_initiate_seeds
    use mpi_lib_ours, only : rank,nodes, bcast_i, bcast_scalar_i, &
                             bcast_integer8, &
                             sync_mpi
    implicit none
    integer line_s,i
    if(.not.first) return
    first = .false.
    
    rnk0 : if(rank==0) then
       write(6,*)' ############ NOW READING RANDOM SEED (KISS ALGORITHM) ############ ' 
       write(6,*)
       write(6,*)'  Using 2 files: '
       write(6,*)'  (I)  counter.inp, with 1 small integer (run counter) modifying seeds '

       inquire(file='counter.inp',exist=fex); if(.not.fex) stop ' ERROR: counter.inp missing '
       open(3,file='counter.inp',status='old')   
       rewind(3)
       iu = 0
       read(3,*,end=33)iu
33     continue
       close(3)

       write(6,*)'         run counter =',iu
       write(6,*)
       write(6,*)'  (II) random.inp,  with  ',ndiffseeds,' lines '
       write(6,*)'         format: line_s,  seed(1:4)          '

       inquire(file='random.inp', exist=fex); if(.not.fex) stop ' ERROR: random.inp  missing '
       open(4, file='random.inp', status='old')    
       rewind(4)

       seed_read = 0
       rank_dependent(:) = 0
       do line_s=1,ndiffseeds
          read(4,*,end=99)     seed_read(:,line_s)
          write(6,*)line_s,seed_read(:,line_s)
       end do
99     continue
       close(4)
       write(6,*)' ' 
       call flush(6)
    end if rnk0

    call bcast_scalar_i(iu)
    call bcast_i(rank_dependent,size(rank_dependent),0)
    call bcast_integer8(seed_read,     size(seed_read),     0)
    call bcast_scalar_i(line_s)
    mdiffseeds = line_s-1
!    call check_lele(1,mdiffseeds,ndiffseeds,' 1 mdiffseeds ndiffseeds ') 
    call check(mdiffseeds,ndiffseeds,' mdiffseeds ndiffseeds ') 

!    write(60+rank,*)' rank_dependent ',rank_dependent(:)
    do line_s=1,mdiffseeds
       if(rank_dependent(line_s)==1) then
          do i=1,nlseed 
             j1 = i
             j2 = 921778613
             j3 = rank
             j4 = line_s
             j5 = 6938938
             seed_array(i,line_s) = &
             seed_read( i,line_s) + j1* j2* j3 + j4 * j5       
!             write(60+rank,*)' j1*j2*j3 ',j1*j2*j3
          enddo
       else
          j1 = line_s
          j2 = 693
          seed_array(:,line_s) = seed_read(:,line_s)+ j1*j2 
       end if

       ! now use counter
       j1= iu
       j2= line_s*7 + line_s*line_s*23
       j3 = 8371637
       seed_array(:,line_s) = seed_array(:,line_s)+ j1*j2*j3
    end do
!    write(60+rank,*)' seed_array_starting ',seed_array
    line_seed = 1
    call ran_ps_putseed(seed_array(:,line_seed))
  end subroutine read_initiate_seeds
end module rand_seed_mod

function ran_ps()
  use mod_kiss,      only : kiss_uniform
  use mpi_lib_ours,  only : rank
  use rand_seed_mod, only : first, line_seed, mdiffseeds, seed_array, read_initiate_seeds
  implicit none
  integer  i
  integer, save :: cnt=0
  integer, parameter :: nrepeat = 50
  real*8   r
  real*8   ran_ps
  
  cnt=cnt+1
  if(cnt>1000)cnt=1000

  if(first) then
     call read_initiate_seeds
     first=.false. 
  end if
  
  call check_lele(1,line_seed, mdiffseeds,' 1 line_seed, mdiffseeds ')

  do i=1,nrepeat
     call ran_ps_putseed(seed_array(:,line_seed))
     call kiss_uniform(r)
     ran_ps = r
     call ran_ps_getseed(seed_array(:,line_seed))
  end do

end function ran_ps

subroutine ran_ps_getseed(seed)
  use mod_kiss, only : kiss_seed_extract_DN
  use rand_seed_mod, only : nlseed
  implicit none
  integer(kind=8) seed(nlseed)
  call check(nlseed, 4     ,' nseed, four   ')
  call kiss_seed_extract_DN(seed(1),seed(2),seed(3),seed(4))
end subroutine ran_ps_getseed

subroutine ran_ps_putseed(seed)
  use mod_kiss, only : kiss_seed
  use rand_seed_mod, only : nlseed
  implicit none
  integer(kind=8) seed(nlseed)
  call check(nlseed, 4     ,' nlseed, four   ')
  call kiss_seed(seed(1), seed(2), seed(3), seed(4))
end subroutine ran_ps_putseed

subroutine ran_ps_putseed_intoarray(seed, line)  ! sets line-seed to line, sets seed_array(:,line) too seed
  use mod_kiss, only : kiss_seed_extract_DN
  use rand_seed_mod, only : nlseed,line_seed,seed_array, mdiffseeds
  implicit none
  integer line
  integer(kind=8) seed(nlseed)
  call check(nlseed, 4     ,' nseed, four   ')
  call check_lele(1,line,mdiffseeds,' 1 line mdifseds ')
  line_seed=line
  seed_array(:,line)=seed
end subroutine ran_ps_putseed_intoarray

function ran_ps_neg()  ! uniform -1 to 1.
  implicit none
  real*8           :: ran_ps_neg
  real*8, external :: ran_ps
  ran_ps_neg = 2d0*ran_ps()-1d0
end function ran_ps_neg

