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
module mpi_lib_ours
  use mpi
  implicit none
  
  save  
  integer :: rank, nodes, merr
  ! colors
  integer :: icolor, color_size, color_rank, color_comm, i1_color
  ! see
  ! http://www.bu.edu/tech/files/text/comm_split_example.txt
  !
  ! or better yet see 
  ! usr/include/mpif-common.h   or google mpif-common.h
  !
!  integer, parameter :: mpi_comm_wd=0 !=MPI_COMM_WORLD
!  integer, parameter :: mpi_double_pn=17 !=MPI_DOUBLE_PRECISION
!  integer, parameter :: mpi_double_cx=22 !=MPI_DOUBLE_COMPLEX
!  integer, parameter :: mpi_ir=7 !=MPI_INTEGER
!  integer, parameter :: mpi_i8=7 !=MPI_INTEGER8
!  integer, parameter :: mpi_r4=14 !=MPI_REAL4
!  integer, parameter :: mpi_cr=5 !=MPI_CHARACHTER
!  integer, parameter :: mpi_sm=3 !MPI_SUM

  integer, parameter :: mpi_comm_wd=MPI_COMM_WORLD
  integer, parameter :: mpi_double_pn=MPI_DOUBLE_PRECISION
  integer, parameter :: mpi_double_cx=MPI_DOUBLE_COMPLEX
  integer, parameter :: mpi_ir=MPI_INTEGER
  integer, parameter :: mpi_i8=MPI_INTEGER8
  integer, parameter :: mpi_r4=MPI_REAL4
  integer, parameter :: mpi_cr=MPI_CHARACTER
  integer, parameter :: mpi_sm=MPI_SUM

contains
  subroutine prepare_mpi_lib
    implicit none
    call mpi_init(merr)                             ; call check(merr,0,' merr, 0 ')
    call mpi_comm_rank(mpi_comm_wd, rank,  merr) ; call check(merr,0,' merr, 0 ')
    call mpi_comm_size(mpi_comm_wd, nodes, merr) ; call check(merr,0,' merr, 0 ')
    !if(rank==0)    write(6,*)' nodes, rank, merr ',nodes,rank,merr
  end subroutine prepare_mpi_lib

  subroutine sync_mpi
    implicit none
    call mpi_barrier(mpi_comm_wd, merr); call check0(merr,' merr, barrier ')
  end subroutine sync_mpi

  subroutine finalize_mpi_lib
    implicit none
    call flush(6)
    call mpi_Finalize(merr)  ;    call check(merr,0,' merr finalize, 0 ')
  end subroutine finalize_mpi_lib

  subroutine sr_mpi_scalar_r8(scalar, src, dest)
    implicit none
    integer src, dest
    real*8 scalar
    real*8 data(1)
    if(rank==src) data = scalar
    call  sr_mpi_r8(data, 1, src, dest)
    if(rank==dest) scalar = data(1)
  end subroutine sr_mpi_scalar_r8

  subroutine sr_mpi_r8(data, nn, src, dest)
    implicit none
    integer nn, src, dest, tag
    real*8  data(nn)
    tag = 100
    if(src==dest) return
    call sync_mpi
    call send_receive_r8_in_mpi(data, nn, src, dest, tag)
  end subroutine sr_mpi_r8

  subroutine sr_mpi_c16(data, nn, src, dest)
    implicit none
    integer nn, src, dest, tag
    complex*16 data(nn)
    tag = 105
    if(src==dest) return
    call sync_mpi
    call send_receive_c16_array_in_mpi(data, nn, src, dest, tag)
  end subroutine sr_mpi_c16

  subroutine send_receive_r8_in_mpi(data, nn, src, dest, tag) ! make tag unique to src & dest !!!!
    implicit none

    integer nn,tag, ir, src,dest
    real*8  data(nn)
    integer status(mpi_status_size) 

    if (rank .eq. src) then;       
       call mpi_send(data,nn, mpi_double_pn,dest ,  tag, mpi_comm_wd, ir )
       call check(ir,0,' mpisnd' )
    else if (rank .eq. dest) then;     
       call mpi_recv(data,nn, mpi_double_pn,src, tag, mpi_comm_wd, status, ir );
       call check(ir,0,' mpisnd' )
    end if
  end subroutine send_receive_r8_in_mpi

  subroutine send_receive_c16_array_in_mpi(data, nn, src, dest, tag) 
    implicit none

    integer nn,tag, ir, src,dest
    complex*16  data(nn)
    integer status(mpi_status_size) 

    if (rank .eq. src) then;       
       call mpi_send(data,nn, mpi_double_cx,dest ,  tag, mpi_comm_wd, ir )
       call check(ir,0,' mpisnd' )
    else if (rank .eq. dest) then;     
       call mpi_recv(data,nn, mpi_double_cx,src, tag, mpi_comm_wd, status, ir );
       call check(ir,0,' mpisnd' )
    end if
  end subroutine send_receive_c16_array_in_mpi

  subroutine send_receive_iarray_in_mpi(data, nn, src, dest)
    implicit none

    integer nn,tag, ir,dest,src, data(nn)
    integer status(mpi_status_size) 
    tag = 3000
    if (rank .eq. src) then;  
       call mpi_send(data,nn, mpi_ir,dest,tag, mpi_comm_wd, ir );          
       call check(ir,0,' mpisnd' )
    elseif(rank .eq. dest) then;
       call mpi_recv(data,nn, mpi_ir,src, tag, mpi_comm_wd, status, ir ); 
       call check(ir,0,' mpircv' )
    end if
  end subroutine send_receive_iarray_in_mpi

  subroutine send_receive_i8_array_in_mpi(data, nn, src, dest)
    implicit none

    integer nn,tag, ir,dest,src
    integer(kind=8) data(nn)
    integer status(mpi_status_size) 
    tag = 3000
    if (rank .eq. src) then;  
       call mpi_send(data,nn, mpi_i8,dest,tag, mpi_comm_wd, ir );    
       call check(ir,0,' mpisnd' )
    elseif(rank .eq. dest) then;
       call mpi_recv(data,nn, mpi_i8,src, tag, mpi_comm_wd, status, ir );  
       call check(ir,0,' mpircv' )
    end if
  end subroutine send_receive_i8_array_in_mpi

  subroutine bcast_i(p, n, iin)
    implicit none
    integer n, iin, ir, p(n)
    call sync_mpi
    call mpi_bcast(p, n, mpi_ir,iin, mpi_comm_wd, ir); call check(ir,0,' mpircv' )
    call sync_mpi
  end subroutine bcast_i

  subroutine bcast_scalar_general_i(p, iin)
    implicit none
    integer iin, ir, p, p1(1), n1
    n1=1
    call sync_mpi
    if(rank==iin) p1(1)=p
    call mpi_bcast(p1, n1, mpi_ir,iin, mpi_comm_wd, ir); call check(ir,0,' mpircv' )
    call sync_mpi
    p=p1(1)
  end subroutine bcast_scalar_general_i

  subroutine bcast_integer8(p, n, iin)
    implicit none
    integer n, iin, ir
    integer(kind=8) p(n)
    call sync_mpi
    call mpi_bcast(p, n, mpi_i8,iin, mpi_comm_wd, ir); call check(ir,0,' mpircv' )
    call sync_mpi
  end subroutine bcast_integer8

  subroutine bcast_r4(p, n, iin)
    implicit none
    integer n, iin, ir
    real*4 p(n)
    call sync_mpi
    call mpi_bcast(p, n, mpi_r4,iin, mpi_comm_wd, ir); call check(ir,0,' mpibr4 ' )
    call sync_mpi
  end subroutine bcast_r4

  subroutine bcast_r8(p, n, iin)
    implicit none
    integer n, iin, ir
    real*8 p(n)
    call sync_mpi
    call mpi_bcast(p, n, mpi_double_pn,iin,mpi_comm_wd,ir); call check(ir,0,' mpibr8 ' )
    call sync_mpi
  end subroutine bcast_r8

  subroutine bcast_scalar_general_r8(p, iin)
    implicit none
    integer n1, iin, ir
    real*8  p, p1(1)
    n1=1
    call sync_mpi
    if(rank==n1) p1(1)=p
    call mpi_bcast(p1,n1,mpi_double_pn,iin,mpi_comm_wd,ir); call check(ir,0,' mpibr8 ' )
    p=p1(1)
    call sync_mpi
  end subroutine bcast_scalar_general_r8


  subroutine bcast_c16(p, n, iin)
    implicit none
    integer n, iin, ir
    complex*16 p(n)
    call sync_mpi
    call mpi_bcast(p, n, mpi_double_cx,iin, mpi_comm_wd, ir); call check(ir,0,' mpibc16' )
    call sync_mpi
  end subroutine bcast_c16

  subroutine bcast_l(l, n, iin)
    implicit none
    logical l(n)
    integer n, iin, st
    integer, allocatable :: il(:)
    allocate(il(n), stat=st); call check0(st,' il ')
    if(rank==0) then
       where(l(:))
          il(:)=1
       elsewhere
          il(:)=0
       endwhere
    endif
    call bcast_i(il, n, iin)
    where(il==1)
       l=.true.
    elsewhere
       l=.false.
    end where
  end subroutine bcast_l

  subroutine bcast_char30(a, iin)
    implicit none
    integer iin,ir
    character*30 a
    call sync_mpi
    call mpi_bcast(a, 30, mpi_cr, iin, mpi_comm_wd, ir);   call check(ir,0,' mpircv' )
    call sync_mpi
  end subroutine bcast_char30

  subroutine bcast_scalar_i(m)
    implicit none
    integer m
    integer ma(1)
    integer n
    ma=0
    if(rank==0) ma=m
    n=1
    call bcast_i(ma,n,0)
    m=ma(1)
  end subroutine bcast_scalar_i

  subroutine bcast_scalar_char(a)
    implicit none
    integer iin,ir
    character a
    character b(1)
    call sync_mpi
    iin=0
    b(1)=a
    call mpi_bcast(b, 1, mpi_cr,iin, mpi_comm_wd, ir)
    call check(ir,0,' mpircv' )
    a=b(1)
    call sync_mpi
  end subroutine bcast_scalar_char

  subroutine bcast_scalar_l(l)
    implicit none
    logical l
    integer m
    integer ma(1)
    integer n
    ma=0
    if(rank==0) then
       if(l) then
          ma = 1
       else
          ma = 0
       endif
    endif
    n=1
    call bcast_i(ma,n,0)
    m=ma(1)
    if(m==1) then
       l=.true.
    else
       l=.false.
    endif
  end subroutine bcast_scalar_l 

  subroutine bcast_scalar_i_gen(m,ir)
    implicit none
    integer m,ir
    integer ma(1)
    integer n
    ma=0
    if(rank==ir) ma=m
    n=1
    call bcast_i(ma,n,ir)
    m=ma(1)
  end subroutine bcast_scalar_i_gen

  subroutine bcast_scalar_r8(x)
    implicit none
    real*8 x
    real*8 xa(1)
    integer n
    xa=0d0
    if(rank==0) xa=x
    n=1
    call bcast_r8(xa,n,0)
    x=xa(1)
  end subroutine bcast_scalar_r8

  subroutine reduce_sum_scalar_r8(a)
    implicit none
    integer n
    real*8  a
    real*8  a1(1)
    a1=a
    call reduce_sum_r8(a1,1,0)
    if(rank==0) a=a1(1)
  end subroutine reduce_sum_scalar_r8

  subroutine reduce_sum_r8(a,n,root)
    implicit none
    integer n,root,ir
    real*8  a(n)
    real*8, allocatable :: b(:)
    allocate(b(n),stat=ir);  if(ir/=0) stop ' breds '
    call sync_mpi
    call mpi_reduce(a, b, n, mpi_double_pn, mpi_sm, root, mpi_comm_wd, ir);if(ir/=0) stop ' alsm '
    call sync_mpi
    if(rank==root) a=b
    deallocate(b)
  end subroutine reduce_sum_r8

  subroutine allsum_scalar_r8(a)
    implicit none
    real*8 a
    real*8 a1(1)
    integer n
    !write(6,*)' stage c rank ',rank
    n=1
    a1(1)=a
    call allsum_r8(a1,n)
    a = a1(1)
  end subroutine allsum_scalar_r8
    
  subroutine allsum_scalar_i(i)
    implicit none
    integer i
    integer i1(1)
    i1=i
    call allsum_i(i1,1)
    i = i1(1)
  end subroutine allsum_scalar_i
    
  subroutine allsum_i(a, n)
    implicit none
    integer n, ir
    integer  a(n)
    integer, allocatable :: b(:)
    allocate(b(n),stat=ir);  if(ir/=0) stop ' bals '
    call sync_mpi
    call mpi_allreduce(a, b, n, mpi_ir, mpi_sm, mpi_comm_wd, ir);if(ir/=0) stop ' alsm '
    call sync_mpi
    a=b
    deallocate(b)
  end subroutine allsum_i

  subroutine allsum_r8(a, n)
    implicit none
    integer n, ir,i
    real*8  a(n)
    real*8, allocatable :: b(:)
    allocate(b(n),stat=ir);  if(ir/=0) stop ' bals '
    !do i=0,max(nodes-1,0)
       !if(i==rank) then;
          !write(6,*)' checking: a,rank,n',a,rank,n; 
          !call flush(6);
       !endif
    !enddo
    call sync_mpi
    !do i=0,max(nodes-1,0)
    !   if(i==rank) then;
    !      write(6,*)' postsync-pre-checking: a,rank,n',a,rank,n; 
    !      call flush(6);
    !   endif
    !enddo
    call mpi_allreduce(a, b, n, mpi_double_pn, mpi_sm, mpi_comm_wd, ir);if(ir/=0) stop ' alsm '
    !write(6,*)' postchecking: rank=',rank; call flush(6);
    call sync_mpi
    a=b
    deallocate(b)
  end subroutine allsum_r8

  subroutine allsum_c16(a, n)
    implicit none
    integer n, ir
    complex*16  a(n)
    complex*16, allocatable :: b(:)
    allocate(b(n),stat=ir);  if(ir/=0) stop ' bals '
    call sync_mpi
    call mpi_allreduce(a, b, n, mpi_double_cx, mpi_sm, mpi_comm_wd, ir);if(ir/=0) stop ' alsm '
    call sync_mpi
    a=b
    deallocate(b)
  end subroutine allsum_c16

  subroutine scatter_r4(p, g, n, root)
    ! n is the small one!!
    implicit none
    integer n, root,ir
    real*4   p(*), g(*)
    if(nodes>1) then
       call sync_mpi
       call mpi_scatter(p, n, mpi_r4, g, n, mpi_r4, root, mpi_comm_wd, ir);
       if(ir/=0) stop ' scatter_r4 '
       call sync_mpi
    else
       g(1:n) = p(1:n)
    end if
  end subroutine scatter_r4

  subroutine scatter_r8(p, g, n, root)
    implicit none
    integer n, root,ir
    real*8   p(*), g(*)
    if(nodes>1) then
       call sync_mpi
       call mpi_scatter(p, n, mpi_double_pn, g, n, mpi_double_pn, root, mpi_comm_wd, ir);
       if(ir/=0) stop ' sct '
       call sync_mpi
    else
       g(1:n) = p(1:n)
    end if
  end subroutine scatter_r8

  subroutine scatter_c16(p, g, n, root)
    implicit none
    integer n, root,ir
    complex*16   p(*), g(*)
    if(nodes>1) then
       call sync_mpi
       call mpi_scatter(p,n,mpi_double_cx,g,n,mpi_double_cx, root, mpi_comm_wd, ir);
       if(ir/=0) stop ' sct c16 '
       call sync_mpi
    else
       g(1:n) = p(1:n)
    end if
  end subroutine scatter_c16

  subroutine gather_i(g, p, n, root)
    implicit none
    integer n, root,ir
    integer   p(*), g(*)
    if(nodes>1) then
       call sync_mpi
       call mpi_gather(g, n, mpi_ir, p, n, mpi_ir, root, mpi_comm_wd, ir);
       if(ir/=0) stop ' gather r8 '
       call sync_mpi
    else
       p(1:n) = g(1:n)
    end if
  end subroutine gather_i

  subroutine gather_r4(g, p, n, root)
    implicit none
    integer n, root,ir
    real*4   p(*), g(*)
    if(nodes>1) then
       call sync_mpi
       call mpi_gather(g, n, mpi_r4, p, n, mpi_r4, root, mpi_comm_wd, ir);
       if(ir/=0) stop ' gather r4 '
       call sync_mpi
    else
       p(1:n) = g(1:n)
    end if
  end subroutine gather_r4


  subroutine gather_r8(g, p, n, root)
    implicit none
    integer n, root,ir
    real*8   p(*), g(*)
    if(nodes>1) then
       call sync_mpi
       call mpi_gather(g, n, mpi_double_pn, &
                       p, n, mpi_double_pn, root, mpi_comm_wd, ir);
       if(ir/=0) stop ' gather r8 '
       call sync_mpi
    else
       p(1:n) = g(1:n)
    end if
  end subroutine gather_r8


  subroutine gather_c16(g, p, n, root)
    implicit none
    integer n, root,ir
    complex*16   p(*), g(*)
    if(nodes>1) then
       call sync_mpi
       call mpi_gather(g, n, mpi_double_cx, p, n, mpi_double_cx, root, mpi_comm_wd, ir);
       if(ir/=0) stop ' gather c16 '
       call sync_mpi
    else
       p(1:n) = g(1:n)
    end if
  end subroutine gather_c16


  !  for using gathev, per https://www.cac.cornell.edu/ranger/mpicc/gathervscatterv.aspx

  subroutine gatherv_r8(g, p, n, nlongmx, root)
    implicit none
    integer n,st,nds, root, nlong, nlongmx, sz1(1), ir, i
    real*8  g(n), p(*) 
    integer, allocatable :: displs(:), size_each(:)
    nds = max(nodes,1)

    !
    ! trivial case, one nodes
    if(nds.le.1) then
       call check(n, nlongmx,' n, nlongmx ')
       p(1:n) = g(1:n)
       return
    end if
    call check_lele(0,root,nds-1,' 0 root nds-1 ')

    ! general case: allocate 
    if(rank==root)then
       allocate(displs(0:nds-1), size_each(0:nds-1),stat=st)
       call check0(st,' disp0 ')
    else
       allocate(displs(0:0), size_each(0:0), stat=st)
       call check0(st,' disp ')
    end if

    !
    ! gather info. on sizes
    sz1(1)=n
    call gather_i(sz1, size_each, 1, root)

    !
    ! create displacement, ensure nlogmx is sufficient
    !
    if(rank==root) then
       !write(6,*)' size_each ',size_each
       nlong = sum(size_each)
       !write(6,*)' nlong ',nlong
       call check(nlong, nlongmx, 'nlong,mx ')
       displs(0) = 0
       do i=0,nds-2 ! not nds-1
          displs(i+1) = displs(i)+size_each(i)
       enddo
    end if

    call sync_mpi
    call mpi_gatherv(g, n, mpi_double_pn, p, size_each, displs, mpi_double_pn, root, mpi_comm_wd, ir)
    call sync_mpi
    call check0(ir,' gatherv_ir ')

    deallocate(displs, size_each)
  end subroutine gatherv_r8

  subroutine scatterv_r8(p, g, nlongmx, n, root)
    implicit none
    integer n,st,nds, root, nlong, nlongmx, sz1(1), ir, i
    real*8  g(n), p(*)
    integer, allocatable :: displs(:), size_each(:)
    nds = max(nodes,1)

    !
    ! trivial case, one nodes
    if(nds.le.1) then
       call check(n, nlongmx,' n, nlongmx ')
       g(1:n) = p(1:n)
       return
    end if
    call check_lele(0,root,nds-1,' 0 root nds-1 ')

    ! general case: allocate 
    allocate(displs(0:nds-1), size_each(0:nds-1),stat=st)
    call check0(st,' disp ')

    !
    ! gather info. on sizes
    sz1(1)=n
    call gather_i(sz1, size_each, 1, root)

    !
    ! create displacement, ensure nlogmx is sufficient
    !
    if(rank==root) then
       nlong = sum(size_each)
       !write(6,*)' scatter size_each ',size_each
       !write(6,*)' scatter nlong , nlongmx ',nlong,nlongmx
       call check(nlong, nlongmx, 'nlong,mx ')
       displs(0) = 0
       do i=0,nds-2 ! not nds-1
          displs(i+1) = displs(i)+size_each(i)
       enddo
       !write(6,*)' displs ',displs
       call flush(6)
    end if
    
    !erase
    !if(rank==0) then; write(6,*)' rank, p(1), p(nlong) ',rank,p(1), p(nlong); call flush(6); 
    !endif

    call sync_mpi
    call mpi_scatterv(p, size_each, displs, mpi_double_pn, g, n, &
                                   mpi_double_pn, root, mpi_comm_wd, ir)
    call check0(ir,' scatterv_ir ')

    call sync_mpi
    !do ir=0,nds-1
     !  call sync_mpi
       !if(ir==0.and.rank==0) then; write(6,*)' rank, p(1), p(nlong) ',rank,p(1), p(nlong); 
    !   call flush(6)
       !end if
       !if(ir==rank) then;write(6,*)' rank, g(1), g(n) ',g(1),g(n); call flush(6)
       !endif
      ! call sync_mpi
    !enddo

    deallocate(displs, size_each)
  end subroutine scatterv_r8

  subroutine prepare_mpi_colors(color_size_in)
    implicit none
    integer color_size_in, ir, key, nds, comm2d
    i1_color = -1
    color_size = color_size_in
    nds = max(nodes,1)
    if(color_size>nds)         stop ' ERROR: color_size too big '
    if(mod(nds,color_size)/=0) stop ' ERROR: mod(nds, color_size) /= 0 '
    icolor     = rank/color_size       ! color designator e.g.,         0 0 1 1 2 2 3 3 4 4 
    color_rank = mod(rank, color_size) ! key, i.e., index within color, 0 1 0 1 0 1 0 1 0 1 ...

    comm2d = mpi_comm_wd
    call mpi_comm_split(comm2d,icolor,color_rank,color_comm,ir); call check(ir,0,' mpicsplit ')
    !
    ! now checking that the key is correct
    !
    call mpi_comm_rank( color_comm, key, ir)
    call check(ir,0,' mpikeyank ')
    call check(color_rank, key, ' color_rank, key ')
  end subroutine prepare_mpi_colors

  subroutine color_allsum_scalar_i(a)
    implicit none
    integer  ir
    integer  a
    integer a1(1)
    a1 = a
    call color_allsum_i(a1,1)
    a = a1(1)
  end subroutine color_allsum_scalar_i


  subroutine color_allsum_i(a, n)
    implicit none
    integer n, ir
    integer  a(n)
    integer, allocatable :: b(:)
    call check(i1_color,-1,' i1_color vs. -1 ')
    allocate(b(n),stat=ir)
    call check0(ir,' b_in_color_allsum ')
    call sync_mpi
    call mpi_allreduce(a, b, n, mpi_ir, mpi_sm, color_comm, ir)
    call check0(ir,' color_allsum_i ')
    call sync_mpi
    a=b
    deallocate(b)
  end subroutine color_allsum_i

  subroutine color_allsum_scalar_r8(a)
    implicit none
    integer ir
    real*8   a
    real*8   a1(1)
    a1 = a
    call color_allsum_r8(a1,1)
    a = a1(1)
  end subroutine color_allsum_scalar_r8

  subroutine color_allsum_r8(a, n)
    implicit none
    integer n, ir
    real*8   a(n)
    real*8, allocatable :: b(:)
    call check(i1_color,-1,' i1_color vs. -1 ')
    allocate(b(n),stat=ir)
    call check0(ir,' b_r8n_color_allsum ')
    call sync_mpi
    call mpi_allreduce(a, b, n, mpi_double_pn, mpi_sm, color_comm, ir)
    call check0(ir,' color_allsum_r8 ')
    call sync_mpi
    a=b
    deallocate(b)
  end subroutine color_allsum_r8

  subroutine color_allsum_scalar_c16(a)
    implicit none
    integer ir
    complex*16   a
    complex*16   a1(1)
    a1 = a
    call color_allsum_c16(a1,1)
    a = a1(1)
  end subroutine color_allsum_scalar_c16

  subroutine color_allsum_c16(a, n)
    implicit none
    integer n, ir
    complex*16   a(n)
    complex*16, allocatable :: b(:)
    call check(i1_color,-1,' i1_color vs. -1 ')
    allocate(b(n),stat=ir)
    call check0(ir,' b_c16n_color_allsum ')
    call sync_mpi
    call mpi_allreduce(a, b, n, mpi_double_cx, mpi_sm, color_comm, ir)
    call check0(ir,' color_allsum_c16 ')
    call sync_mpi
    a=b
    deallocate(b)
  end subroutine color_allsum_c16

  subroutine color_reduce_sum_r8(a,n,root)
    implicit none
    integer n,root,ir
    real*8  a(n)
    real*8, allocatable :: b(:)
    call check(i1_color,-1,' i1_color vs. -1 ')
    allocate(b(n),stat=ir);  if(ir/=0) stop ' ERROR:color_breds '
    call sync_mpi
    call mpi_reduce(a, b, n, mpi_double_pn, mpi_sm, root, color_comm, ir)
    if(ir/=0) stop ' ERROR: color_alsm '
    call sync_mpi
    if(color_rank==root) a=b
    deallocate(b)
  end subroutine color_reduce_sum_r8

  subroutine color_reduce_sum_c16(a,n,root)
    implicit none
    integer n,root,ir
    complex*16  a(n)
    complex*16, allocatable :: b(:)
    call check(i1_color,-1,' i1_color vs. -1 ')
    allocate(b(n),stat=ir);  if(ir/=0) stop ' ERROR:color_bredsc '
    call sync_mpi
    call mpi_reduce(a, b, n, mpi_double_cx, mpi_sm, root, color_comm, ir)
    if(ir/=0) stop ' ERROR: color_alsmc '
    call sync_mpi
    if(color_rank==root) a=b
    deallocate(b)
  end subroutine color_reduce_sum_c16

  subroutine color_bcast_r8(p, n, iin)
    implicit none
    integer n, iin, ir
    real*8 p(n)
    call check(i1_color,-1,' i1_color vs. -1 ')
    call sync_mpi
    call mpi_bcast(p, n, mpi_double_pn,iin,color_comm,ir); call check(ir,0,' mpibr8 ' )
    call sync_mpi
  end subroutine color_bcast_r8
  
  subroutine color_bcast_c16(p, n, iin)
    implicit none
    integer n, iin, ir
    complex*16 p(n)
    call check(i1_color,-1,' i1_color vs. -1 ')
    call sync_mpi
    call mpi_bcast(p, n, mpi_double_cx,iin, color_comm, ir); call check(ir,0,' mpibc16' )
    call sync_mpi
  end subroutine color_bcast_c16

  subroutine color_bcast_c8(p, n, iin)
    implicit none
    integer n, iin, ir
    complex*8 p(n)
    call check(i1_color,-1,' i1_color vs. -1 ')
    call sync_mpi
    call mpi_bcast(p, n, mpi_complex,iin, color_comm, ir); call check(ir,0,' mpibc8' )
    call sync_mpi
  end subroutine color_bcast_c8


  subroutine color_scatter_r8(p, g, n, iin)
    implicit none
    integer n, iin, ir
    real*8 p(*), g(*)

    call check(i1_color,-1,' i1_color vs. -1 ')

    if(color_size>1) then
       call sync_mpi
       call mpi_scatter(p, n, mpi_double_pn, &
                        g, n, mpi_double_pn, iin, color_comm, ir)
       call check(ir,0,' mpisctr8_color ' )
       call sync_mpi
    else
       g(1:n)=p(1:n)
    end if

  end subroutine color_scatter_r8

  subroutine color_gather_r8(p, g, n, iin)
    implicit none
    integer n, iin, ir
    real*8 p(*), g(*)

    call check(i1_color,-1,' i1_color vs. -1 ')

    if(color_size>1) then

       call sync_mpi
       call mpi_gather(p, n, mpi_double_pn, &
                       g, n, mpi_double_pn, iin, color_comm, ir)
       call check(ir,0,' mpisctr8_color ' )
       call sync_mpi
    else
       g(1:n)=p(1:n)
    end if

  end subroutine color_gather_r8

  subroutine color_scatter_c16(p, g, n, iin)
    implicit none
    integer n, iin, ir
    complex*16 p(*), g(*)

    call check(i1_color,-1,' i1_color vs. -1 ')

    if(color_size>1) then
       call sync_mpi
       call mpi_scatter(p, n, mpi_double_cx, &
                        g, n, mpi_double_cx, iin, color_comm, ir)
       call check(ir,0,' mpisctc16_color ' )
       call sync_mpi
    else
       g(1:n)=p(1:n)
    end if

  end subroutine color_scatter_c16

  subroutine color_gather_c16(p, g, n, iin)
    implicit none
    integer n, iin, ir
    complex*16 p(*), g(*)

    call check(i1_color,-1,' i1_color vs. -1 ')

    if(color_size>1) then
       call sync_mpi
       call mpi_gather(p, n, mpi_double_cx, &
                       g, n, mpi_double_cx, iin, color_comm, ir)
       call check(ir,0,' mpisctc16_color ' )
       call sync_mpi
    else
       g(1:n)=p(1:n)
    end if

  end subroutine color_gather_c16

  subroutine color_gather_c8(p, g, n, iin)
    implicit none
    integer n, iin, ir
    complex*8 p(*), g(*)

    call check(i1_color,-1,' i1_color vs. -1 ')

    if(color_size>1) then
       call sync_mpi
       call mpi_gather(p, n, mpi_complex, &
                       g, n, mpi_complex, iin, color_comm, ir)
       call check(ir,0,' mpisctc8_color ' )
       call sync_mpi
    else
       g(1:n)=p(1:n)
    end if

  end subroutine color_gather_c8

end module mpi_lib_ours


subroutine prepare_mpi
  use mpi_lib_ours, only : prepare_mpi_lib
  implicit none
  call prepare_mpi_lib
end subroutine prepare_mpi

subroutine finalize_mpi
  use mpi_lib_ours, only : finalize_mpi_lib
  implicit none
  call finalize_mpi_lib
end subroutine finalize_mpi

