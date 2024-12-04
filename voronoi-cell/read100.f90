 subroutine calc_distan(length, lengthh, distan)

  implicit none
  real(8), intent(in) :: length(1:3), lengthh(1:3)
  real(8), intent(inout) :: distan(1:4)
  integer :: i1


  do i1 = 1, 3
   if(abs(distan(i1)) > lengthh(i1)) then
    if(distan(i1) > 0.d0) then
     distan(i1) = distan(i1) - length(i1)
    else
     distan(i1) = distan(i1) + length(i1)
    endif
   endif
  enddo
  distan(4) = distan(1) * distan(1) + distan(2) * distan(2) + distan(3) * distan(3)
  end subroutine calc_distan

 recursive subroutine quicksort_index(a, first, last, num,b)
  implicit none
  integer, intent(in) :: first, last, num
  real(8), intent(inout) ::  a(1:num)
  integer, intent(inout) :: b(1:num)
  real(8) :: x, t
  integer :: i, j, t2

  x = a( int((first+last) / 2) )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     enddo
     do while (x < a(j))
        j=j-1
     enddo
     if (i >= j) exit
      t = a(i);
      a(i) = a(j);
      a(j) = t
      t2 = b(i)
      b(i)=b(j)
      b(j)=t2
      i=i+1
      j=j-1
  enddo
  if (first < i-1) call quicksort_index(a, first, i-1, num, b)
  if (j+1 < last)  call quicksort_index(a, j+1, last, num, b)
  end subroutine quicksort_index

program rea
use te
implicit none
integer :: i0, i1, i4, isave, num_atom, index_atom, index_type, order(1:1023), m, num_o, num_h
real(8) :: length(1:3), lengthh(1:3), vertex(1:3,1:3), pos(1:4), poso(1:3,1:1024), vec(1:3), pair(1:1023)
real(8) :: vol,  posh(1:3,1:2048), tot_vol(1:9), tot_surf, coeff
integer :: num_order, initial_snapshot, input0, input5, last_snapshot, ind_snapshot, start_snapshot
integer :: io, ih, oh_pair(1:3,1:1024), c(1:8)
integer :: p0, i5, num_near(1:8), index_run
character(len=70) :: chac, cha(1:8)
logical :: same, closure
real(8) :: vol_index(1:8), tot_vol_index(1:8)
integer :: near_neighbor(1:1023), near(1:1023,1:8), num_near_neighbor
integer :: tot_num_near(1:9), tot_atom(1:9)
200 format(I8, 2X, I5, 2X, I3, 2X, I3, 2X, I3, 2X, I3, 2X, I3, 2X, I3, 2X, I3, 2X, I3)
300 format(I8, 2X, I5, 2X, F17.11, 2X, F17.11, 2X, F17.11, 2X, F17.11, 2X, F17.11, 2X, F17.11, 2X, F17.11, 2X, F17.11, 2X, F17.11 )
400 format(I8, 2X, I5, 2X, I3, 2X, I3, 2X, I3, 2X, I3, 2X, I3, 2X, I3, 2X, I3, 2X, I3, 2X, I3)
410 format(A12, 2X, F17.11, 2X, F17.11, 2X, F17.11, 2X, F17.11)
! 3.26 as cut
tot_vol=0.0
tot_vol_index=0.0
tot_num_near=0
tot_atom=0
coeff=(2*1.00794+15.994)*10/(6.02)
do isave=5, 7
 write(chac,*) isave
 index_run=0
 !!write(ch,*) index_run
 !!open(15, file="input"//trim(adjustl(chac))//"-"//trim(adjustl(ch)))
 open(15, file="input"//trim(adjustl(chac)))
 open(17, file="output"//trim(adjustl(chac)))
 open(10, file="../../save"//trim(adjustl(chac))//"/water.dump")
 open(70, file="tot_num_"//trim(adjustl(chac))//".dat")
 open(80, file="tot_volume_"//trim(adjustl(chac))//".dat")
 open(90, file="closure_"//trim(adjustl(chac))//".dat")
 read(15, *)
 read(15,*) input0
 if (input0 .gt. 0 ) then
  last_snapshot=input0
  write(17,*) "index of last snapshot ", last_snapshot
 endif
 read(15,*)
 read(15,*) input5
 if ( input5 .gt. 0) then
  initial_snapshot=input5
 endif
 close(15)
 write(17,*) "initial_snapshot ", initial_snapshot

 if ( initial_snapshot .lt. 100000 ) then
  start_snapshot=100000
  do ind_snapshot=initial_snapshot, start_snapshot-1
   do i0=1, 3081
    read(10,*)
   enddo
  enddo
 else
  start_snapshot=initial_snapshot
 endif
 write(17, *) "start_snapshot  ", start_snapshot
 write(17, *) "tot_num_snapshotsused in calculation"
 write(17, *) last_snapshot-start_snapshot+1
 
 !==begin reading file and calculation
 do ind_snapshot=start_snapshot, last_snapshot
 !print *, ind_snapshot
  if ( mod(ind_snapshot, 1000) .eq. 0 ) then
   write(17,*) ind_snapshot
  endif
  if ( mod(ind_snapshot, 100) .eq. 0 ) then
   pos=0.0
   poso=0.0
   posh=0.0
   do i0=1, 3
    read(10,*)
   enddo
   read(10,*) num_atom
   num_o=(num_atom)/3
   num_h=num_o*2
   read(10,*)
!read length
   do i0=1, 3
    read(10,*) vertex(1:3,i0)
    length(i0)=vertex(2,i0)-vertex(1,i0)
    lengthh(i0)=length(i0)/(2.0)
   enddo
   read(10,*)
  !read atom position
   io=0
   ih=0
   do i0=1, num_atom
    read(10,*) index_atom, index_type, pos(1:3), vec(1:3)
    if (index_type .eq. 1) then
     io=io+1
     poso(1:3,io)=pos(1:3)
    endif
    if (index_type .eq. 2) then
     ih=ih+1
     posh(1:3,ih)=pos(1:3)
    endif
   enddo 
   oh_pair=0
   call ohpair(poso, posh, num_o, num_h, length, lengthh, calc_distan, oh_pair)

   do i0=1, num_o
!=====sort ascending orderi
    vol=0.0
    tot_surf=0.0
    m=0
    pair=0.0
    order=0
    vol_index=0.0
    tot_atom(1)=tot_atom(1)+1
    c=0
    do i1=1, num_o
     if ( i0 .ne. i1) then
      m=m+1
      pos(1:3)=poso(1:3,i1)-poso(1:3,i0)
      call calc_distan(length, lengthh, pos)
      pair(m)=sqrt(pos(4))
      order(m)=i1
     endif
    enddo
    call quicksort_index(pair, 1, num_o-1, num_o-1, order)
!==========
    num_order=num_o-1
    num_near_neighbor=0
    near_neighbor=0
    num_near=0
    near=0
    call initialization2(length, lengthh, poso, i0, order(1:num_order), num_order, near_neighbor(1:num_order), vol, tot_surf)
    tot_vol(1)=tot_vol(1)+vol
  !==== near_neighbor part
    do i4=1, 1023
     if ( near_neighbor(i4) .ne. 0 ) then
      num_near_neighbor=num_near_neighbor+1
     endif 
    enddo
    tot_num_near(1)=tot_num_near(1)+num_near_neighbor
    do i4=1, num_near_neighbor
     pos(1:3)=poso(1:3,near_neighbor(i4))-poso(1:3,i0)
     call calc_distan(length, lengthh, pos) 
     if ( pos(4) .le. (4)**2 ) then
      num_near(1)=num_near(1)+1
      near(num_near(1),1)=near_neighbor(i4)
     endif
     if ( pos(4) .le. (4.4)**2 ) then
      num_near(2)=num_near(2)+1
      near(num_near(2),2)=near_neighbor(i4)
     endif
     if ( pos(4) .le. (4.5)**2 ) then
      num_near(3)=num_near(3)+1
      near(num_near(3),3)=near_neighbor(i4)
     endif
     if ( pos(4) .le. (5)**2 ) then
      num_near(4)=num_near(4)+1
      near(num_near(4),4)=near_neighbor(i4)
     endif
     if ( pos(4) .le. (5.5)**2 ) then
      num_near(5)=num_near(5)+1
      near(num_near(5),5)=near_neighbor(i4)
     endif
     if ( pos(4) .le. (6)**2 ) then
      num_near(6)=num_near(6)+1
      near(num_near(6),6)=near_neighbor(i4)
     endif
     if ( pos(4) .le. (6.5)**2 ) then
      num_near(7)=num_near(7)+1
      near(num_near(7),7)=near_neighbor(i4)
     endif
     if ( pos(4) .le. (7)**2 ) then
      num_near(8)=num_near(8)+1
      near(num_near(8),8)=near_neighbor(i4)
     endif
    enddo
    do i4=1, 8
     closure=.true.
     call initialization(length, lengthh, poso, i0, near(1:num_near(i4),i4), num_near(i4), vol_index(i4), tot_surf, closure)
     if ( closure .eq. .true. ) then
      c(i4)=1
     endif
     if ( c(i4) .eq. 1 ) then
      tot_num_near(i4+1)=tot_num_near(i4+1)+num_near(i4)
      tot_vol(i4+1)=tot_vol(i4+1)+vol
      tot_atom(i4+1)=tot_atom(i4+1)+1
      tot_vol_index(i4)=tot_vol_index(i4)+vol_index(i4)
     endif
    enddo
   
    
    write(70,400) ind_snapshot, i0, num_near_neighbor, num_near(1:8)
    write(80,300) ind_snapshot, i0, vol, vol_index(1:8)
    write(90,200) ind_snapshot, i0, c(1:8)
   enddo
  else
   do i4=1, 3081
    read(10,*)
   enddo
  endif
 enddo
 close(17)
 close(10)
 close(70)
 close(90)
 close(80)
enddo
cha(1)="4    "
cha(2)="4.4  "
cha(3)="4.5  "
cha(4)="5    "
cha(5)="5.5  "
cha(6)="6    "
cha(7)="6.5  "
cha(8)="7    "
open(95, file="average_volume.dat")
write(95,*) "tot  ", tot_vol(1)/(tot_atom(1)+0.0), tot_num_near(1)/(tot_atom(1)+0.0)
do i4=1, 8
 write(95,410) cha(i4), tot_vol(i4+1)/(tot_atom(i4+1)+0.0), tot_vol_index(i4)/(tot_atom(i4+1)+0.0), tot_num_near(i4+1)/(tot_atom(i4+1)+0.0), tot_atom(i4+1)/(tot_atom(1)+0.0)
enddo
close(95)
end program
