module cut
implicit none

type :: edge
integer :: f(1:2), v(1:2), stat
end type

type :: vert
real(8) :: posv(1:4)
integer :: e(1:3), stat
end type

type :: flist
integer :: e, link, v
end type

type :: face
real(8) :: dist
integer :: fptr, stat, vfar
end type
contains

!=========== calculate the distance between two oxygen atoms 
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

!======== calculate the plane which bisects two atoms
 subroutine bisectplane( pos0, pos1, length, lengthh, vert_te, edge_te, face_te , flist_te, elast, vlast, flast, flistlast, near)
 

 real(8), intent(in) :: pos1(1:3), pos0(1:3), length(1:3), lengthh(1:3)
 real(8) :: pos(1:4), pos_plus(1:4), new_vertex(1:4), a, temp0(1:4)
 integer, intent(inout) :: elast, vlast, flast, flistlast
 type(vert), intent(inout) :: vert_te(1:2000)
 type(edge), intent(inout) :: edge_te(1:2000)
 type(face), intent(inout) :: face_te(1:50)
 type(flist), intent(inout) :: flist_te(1:5000)
 logical :: vert_outside
 integer :: i0, i1, fs, fe, fv, firstlink, link_index, first_vertex, second_vertex, i2, elast0, link_index2
 integer :: first_cut_edge, second_cut_edge, second_cut_flist, vlast0, m, tempi
 integer ::  changed_face
 logical, intent(inout) :: near
 ! find the vertex, edge, face flist influced by the cutting edge
 vlast0=vlast
 changed_face=0
 near=.false.
 m=0
 do i1=1, vlast
  if ( vert_te(i1).stat .eq. 2 ) then
   m=m+1
  endif
 enddo
 m=0
 call judge_stat(pos0, pos1, length, lengthh, vert_te, edge_te, face_te, flist_te, vlast, flast) 

 !==========finding the new vert 
 !=======changed item: vlast,  new vert_te : stat, posv(1:4), e(1)   cutting edge_te: v
 m=0
 do i1=1, vlast
  if ( vert_te(i1).stat .eq. 2 ) then
   m=m+1
  endif
 enddo
 !!write(*,*) "m: ", m
 m=0
!=====whether add the new face to the current Vc after the consideration of xtra neighbors
 call find_new_vert(pos0, pos1, length, lengthh, vert_te, edge_te, elast, vlast, changed_face)
 if ( changed_face .gt. 0 ) then
  near=.true.
  elast0=elast
  do i0=1, flast
   if ( face_te(i0).stat .eq. 2 ) then
   !===set the accessory of the cutting face
   !====change face_te(i0).fptr  index of first_vertex, first_cut_edge, second_vertex, second_cut_flist, second_cut_edge   !==new edge(only 1): stat, v(1:2), f(1)  new flist(only 2): v, e, link   new flist loop
   !=== new flistlast, elast: 
    call cutting_face_accessory( face_te(i0), i0, edge_te, vert_te, flist_te, flistlast, elast, flast, first_vertex, first_cut_edge, second_cut_flist, second_vertex, second_cut_edge)
   endif
  enddo
  !======the new flist to new surface flist_te new added  new face : stat fptr 
  call new_face_flist(flist_te, vert_te, edge_te, face_te(flast+1), vlast0, vlast, elast0, flistlast )
 ! flistlast changed
  temp0(1:3)=pos1(1:3)-pos0(1:3)
  call calc_distan(length, lengthh, temp0)
 ! new face : dist
  face_te(flast+1).dist=sqrt(temp0(4))/(2.0)
 !!write(*,*) "new-face-dist: ", face_te(flast+1).dist
  flistlast=flistlast+vlast-vlast0
  flast=flast+1
 !================changw the stat of cutting face and cutting edge
  do i0=1, elast
   if ( edge_te(i0).stat .eq. 2 ) then
    edge_te(i0).stat = 3
   endif
  enddo   
  do i0=1, flast
   if ( face_te(i0).stat .eq. 2 ) then
    face_te(i0).stat = 3
   endif
  enddo
 endif
 !!do i0=1, vlast
  !!if ( vert_te(i0).stat .eq. 2 ) then
   !!write(*,*) "final vert: ", i0, vert_te(i0).posv(1:3)
  !!endif
 !!enddo
 end subroutine bisectplane 

! adding the new face 
 subroutine new_face_flist(flist_te, vert_te, edge_te, face_te_new, vlast0, vlast, elast0, flistlast )
 type(flist), intent(inout) :: flist_te(1:5000)
 type(vert), intent(in) :: vert_te(1:2000) 
 type(edge), intent(in) :: edge_te(1:2000)
 type(face), intent(inout) :: face_te_new
 integer, intent(in) :: vlast0, vlast, elast0, flistlast
 integer :: i0, i1, i2, m, tempi
  i0=vlast0+1
 !write(*,*) "new vertex: ", i0,  vert_te(i0).e
  m=1
  do i1=1, 3
   if ( vert_te(i0).e(i1) .gt. elast0 ) then
    tempi=vert_te(i0).e(i1)
   endif
  enddo
  face_te_new.fptr=flistlast+1
  face_te_new.stat=3
  do while ( m .le. (vlast-vlast0) )
   if ( m .gt. 1 ) then
    do i1=1, 3
     if ((vert_te(i0).e(i1) .gt. elast0) .and. ( vert_te(i0).e(i1) .ne. flist_te(flistlast+m-1).e ))  then
      tempi=vert_te(i0).e(i1)
     endif
    enddo
   endif 
   flist_te(flistlast+m).v=i0
   flist_te(flistlast+m).e=tempi
   if ( m .eq. (vlast-vlast0) ) then
    flist_te(flistlast+m).link=flistlast+1
   else
    flist_te(flistlast+m).link=flistlast+m+1
   endif

   do i2=1, 2
    if (edge_te(tempi).v(i2) .ne. i0 )  exit
   enddo
   i0=edge_te(tempi).v(i2)
   m=m+1
  enddo
 
 do i0=1, vlast-vlast0
  !!write(*,*) "new-face"
  !!write(*,*) "v: ", flist_te(flistlast+i0).v
  !!write(*,*) "e: ", flist_te(flistlast+i0).e
  !!write(*,*) "link: ", flist_te(flistlast+i0).link
 enddo

 end subroutine new_face_flist

! the revised face of VC due to the exisitance of new face
 subroutine cutting_face_accessory( face_te0, index_face, edge_te, vert_te, flist_te, flistlast, elast, flast, first_vertex, first_cut_edge, second_cut_flist, second_vertex, second_cut_edge ) 
 type(vert), intent(inout) :: vert_te(1:2000)
 type(edge), intent(inout) :: edge_te(1:2000)
 type(flist), intent(inout) :: flist_te(1:5000)
 type(face), intent(inout)    :: face_te0
 integer, intent(inout) :: flistlast, elast, first_vertex, first_cut_edge, second_cut_flist, second_vertex
 integer, intent(inout) :: second_cut_edge
 integer, intent(in) :: index_face, flast
 integer :: firstlink, fv, link_index, fe, i1
  firstlink=0
  firstlink=face_te0.fptr
  fv=flist_te(face_te0.fptr).v
  if ( vert_te(fv).stat .eq. 1 ) then
   do while ( vert_te(fv).stat .eq. 1 )
    firstlink=flist_te(firstlink).link
    fv=flist_te(firstlink).v
   enddo
  endif
  face_te0.fptr=firstlink
 !====finding the breaking point first_vertex
  link_index=firstlink
  fe=flist_te(link_index).e
  do while ( edge_te(fe).stat .eq. 3 )
   link_index=flist_te(link_index).link
   fe=flist_te(link_index).e
  enddo
  do i1=1, 2
   if (edge_te(fe).v(i1) .ne. flist_te(link_index).v ) then
    first_vertex=edge_te(fe).v(i1)
    first_cut_edge=fe
   endif
  enddo
 !===== second_vertex
  link_index=flist_te(link_index).link
  if ( edge_te(flist_te(link_index).e).stat .eq. 1 ) then
   do while ( edge_te(flist_te(link_index).e).stat .eq. 1 )
    link_index=flist_te(link_index).link
   enddo
  endif
  second_cut_flist=flist_te(link_index).link
  do i1=1, 2
   if ( edge_te(flist_te(link_index).e).v(i1) .ne. flist_te(flist_te(link_index).link).v ) then
    second_vertex=edge_te(flist_te(link_index).e).v(i1)
    second_cut_edge=flist_te(link_index).e
   endif
  enddo
 !!write(*,*) "cutvertex ",  first_vertex, second_vertex
 
 !=====creating new edge
   edge_te(elast+1).stat=3
   edge_te(elast+1).v(1)=first_vertex
   edge_te(elast+1).v(2)=second_vertex
   edge_te(elast+1).f(1)=index_face
   edge_te(elast+1).f(2)=flast+1
 !====add new edge of new vertex
   do i1=1, 3
    if ( vert_te(first_vertex).e(i1) .eq. 0 ) exit
   enddo
   vert_te(first_vertex).e(i1)=elast+1
   do i1=1, 3
    if ( vert_te(second_vertex).e(i1) .eq. 0) exit
   enddo
   vert_te(second_vertex).e(i1)=elast+1
 !==== create two new element of flist for face i0
   flist_te(flistlast+1).v=first_vertex
   flist_te(flistlast+1).e=elast+1
   flist_te(flistlast+1).link=flistlast+2
   flist_te(flistlast+2).v=second_vertex
   flist_te(flistlast+2).e=second_cut_edge
   flist_te(flistlast+2).link=second_cut_flist

 !========create the new loop of flist for face i0
 !==========link_index=face_te(i0).fptr
   link_index=face_te0.fptr
   do while( edge_te(flist_te(link_index).e).stat .eq. 3 )
    link_index=flist_te(link_index).link
   enddo
   flist_te(link_index).link=flistlast+1
 !==== print result
   link_index=face_te0.fptr
 if ( flist_te(link_index).link .ne. face_te0.fptr ) then
  do
   !!write(*,*) "link_index: ", link_index
   !!write(*,*) "v,e : ", flist_te(link_index).v, flist_te(link_index).e
   !!write(*,*) "link: ", flist_te(link_index).link
   if ( flist_te(link_index).link .eq. face_te0.fptr ) exit
   link_index=flist_te(link_index).link
  enddo
 endif
   flistlast=flistlast+2
   elast=elast+1
 end subroutine cutting_face_accessory

 subroutine find_new_vert(pos0, pos1, length, lengthh, vert_te, edge_te, elast, vlast, changed_face)
 real(8), intent(in) ::  pos1(1:3), pos0(1:3), length(1:3), lengthh(1:3)
 type(vert), intent(inout) :: vert_te(1:2000)
 type(edge), intent(inout) :: edge_te(1:2000)
 integer, intent(in) :: elast
 integer, intent(inout) :: vlast
 integer, intent(inout) :: changed_face
 integer :: i0, i1, i2, change
 real(8) :: new_vertex(1:4), pos(1:4), a
! change determine whether the cutting edge has changed the end point
 do i0=1, elast
  if ( edge_te(i0).stat .eq. 2 ) then
   change=0
   changed_face=1
   do i1=1, 2
    if ( (vert_te(edge_te(i0).v(i1)).stat .eq. 2 ) .and. ( change .eq. 0 ) ) then
     if ( i1 .eq. 1) then
      i2=2
     else
      i2=1
     endif
     !!write(*,*) "saved point: ", vert_te(edge_te(i0).v(i1)).posv(1:3)
     call new_vertex_coeff(pos0, pos1, vert_te(edge_te(i0).v(i1)).posv(1:3), vert_te(edge_te(i0).v(i2)).posv(1:3), length, lengthh, a)
     new_vertex(1:3)=vert_te(edge_te(i0).v(i1)).posv(1:3)+a*(vert_te(edge_te(i0).v(i2)).posv(1:3)-vert_te(edge_te(i0).v(i1)).posv(1:3))
     vert_te(vlast+1).posv(1:3)=new_vertex(1:3)
     vert_te(vlast+1).stat=2
     vert_te(vlast+1).e(1)=i0
     !!write(*,*) "new_vertex: ", vert_te(vlast+1).posv(1:3), vert_te(vlast+1).stat
     ! new added
     vert_te(vlast+1).e(2:3)=0
     pos(1:3)=vert_te(vlast+1).posv(1:3)-pos0(1:3)
     call calc_distan(length, lengthh, pos)
     vert_te(vlast+1).posv(4)=sqrt(sqrt(pos(4)))
     change=change+1
     vlast=vlast+1
     edge_te(i0).v(i2)=vlast
     !write(*,*) "edge: ", i0, i1, vlast
    endif
   enddo
  endif
 enddo
 end subroutine find_new_vert

 subroutine judge_stat(pos0, pos1, length, lengthh, vert_te, edge_te, face_te, flist_te, vlast, flast)
 type(vert), intent(inout) :: vert_te(1:2000)
 type(edge), intent(inout) :: edge_te(1:2000)
 type(face), intent(inout) :: face_te(1:50)
 type(flist), intent(inout) :: flist_te(1:5000)
 integer, intent(in) :: vlast, flast
 real(8), intent(in) ::  pos1(1:3), pos0(1:3), length(1:3), lengthh(1:3) 
 integer :: i0, i1, fs
 logical  :: vert_outside
 do i0=1, vlast
  if ( vert_te(i0).stat .gt. 1 ) then
   vert_outside=.false.
   call judge_vertex( pos0, pos1, vert_te(i0).posv(1:3), length, lengthh, vert_outside)
   if (vert_outside .eq. .true. ) then
    vert_te(i0).stat=1
    do i1=1, 3
     edge_te(vert_te(i0).e(i1)).stat=edge_te(vert_te(i0).e(i1)).stat-1
    enddo
   endif
  endif
 enddo
 
 do i0=1, flast
  fs=0
  if ( face_te(i0).stat .gt. 1 ) then
   i1=face_te(i0).fptr
   if ( flist_te(i1).link .ne. face_te(i0).fptr ) then
    do
   !write(*,*) "i0, i1 ",  i0, i1
     if ( vert_te(flist_te(i1).v).stat .lt. 2 ) then
      face_te(i0).stat=2
     !write(*,*) "stat=2 i0i1: ", i0, i1, flist_te(i1).v
     else
      fs=fs+1
     endif
     i1=flist_te(i1).link
    !write(*,*) "i1 ", i1
     if ( i1 .eq. face_te(i0).fptr )  exit
    enddo
   endif
   if ( fs .lt. 1) then
    face_te(i0).stat=1
   endif
  endif
 enddo
 end subroutine judge_stat

 
 subroutine new_vertex_coeff(pos0, pos1, pos2, pos3, length, lengthh, a) 
 implicit none
 real(8), intent(in) :: pos0(1:3), pos1(1:3), pos2(1:3), pos3(1:3), length(1:3), lengthh(1:3)
 real(8), intent(inout) :: a
 real(8) :: distan(1:4), shift1(1:4), shift2(1:3), temp(1:4), temp0, temp1, b, v1_v2(1:4), pv1(1:3), temp2, temp3
 !==== distan: alpha , pv1: pos2, pv1   pv2: pos3, shift2: ri pos0: rj  v1_v2: v2-v1
 distan(1:3)=pos1(1:3)-pos0(1:3)
 call calc_distan(length, lengthh, distan)
 shift1(1:3)=pos3(1:3)-pos2(1:3)
 !call calc_distan(length, lengthh, shift1)
 shift2(1:3)=pos0(1:3)+distan(1:3)
 call dotproduct(shift2, shift2, temp1)
 call dotproduct(pos0, pos0, temp0)
 b=(temp1-temp0)/(2.0)
 !v1_v2(1:3)=pos3(1:3)-pos2(1:3)
 !call calc_distan(length, lengthh, v1_v2)
 temp(1:3)=pos2(1:3)-pos0(1:3)
 !call calc_distan(length, lengthh, temp)
 pv1(1:3)=temp(1:3)+pos0(1:3)
 call dotproduct(distan(1:3), shift1(1:3), temp2)
 call dotproduct(distan(1:3), pv1(1:3), temp3)
 a=(b-temp3)/(temp2+0.0)
 end subroutine new_vertex_coeff

 subroutine judge_vertex( pos0, pos1, pos2, length, lengthh, vert_outside)
 implicit none
 real(8), intent(in) :: pos1(1:3), pos0(1:3), pos2(1:3), length(1:3), lengthh(1:3)
 logical, intent(inout) :: vert_outside
 real(8) :: distan(1:4), plus(1:4), plus2(1:4), tempr, temp0, temp1, distan1(1:4)
 vert_outside= .false.
 distan(1:3)=pos1(1:3)-pos0(1:3)
 distan1(1:3)=pos2(1:3)-pos0(1:3)
 call calc_distan(length, lengthh, distan)
 !call calc_distan(length, lengthh, distan1)
 plus(1:3)=distan(1:3)+pos0(1:3)
 plus2(1:3)=distan1(1:3)+pos0(1:3)
 call dotproduct(pos0(1:3), pos0(1:3), temp0)
 call dotproduct(plus(1:3), plus(1:3), temp1)
 call dotproduct(distan(1:3), plus2(1:3), tempr)
 !write(*,*) pos0
 !write(*,*) plus
 !write(*,*) plus2
 if ( tempr .gt. ((temp1-temp0)/(2.0)) ) then
  vert_outside= .true.
  !!write(*,*) "outside: ", vert_outside
 endif
 end subroutine judge_vertex
 
 
 subroutine dotproduct(a, b, c)
 implicit none
 real(8), intent(in) :: a(3), b(3);
 real(8), intent(inout) :: c;
 c=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 end subroutine dotproduct
 
 subroutine crossproduct(a, b, c)
  implicit none
  real(8), intent(in) :: a(3), b(3);
  real(8), intent(inout) :: c(3);
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
  end subroutine crossproduct 
 
!=========calculate the volume of VC if we know their faces and vertices 
  subroutine volume(face_te, flist_te, vert_te, flast, vol, tot_surf, length, lengthh)
  implicit none
  type(vert), intent(in) :: vert_te(1:2000)
  type(face), intent(in) :: face_te(1:50)
  type(flist), intent(in) :: flist_te(1:5000)
  integer, intent(in)  ::  flast
  real(8), intent(inout) :: vol, tot_surf
  real(8), intent(in) :: length(1:3), lengthh(1:3)
  real(8)  :: surf
  integer :: i0
  vol=0.0
  tot_surf=0.0
  do i0=1, flast
   if ( face_te(i0).stat .eq. 3) then
    call surface(face_te(i0), flist_te, vert_te, surf, length, lengthh)
    !!write(*,*) "surface: ", i0, surf
    tot_surf=tot_surf+surf
    vol = vol + (surf*face_te(i0).dist)/3.0
   endif
  enddo
  end subroutine

  subroutine surface(face0, flist_te, vert_te, surf, length, lengthh)
  implicit none
  type(face), intent(in) :: face0
  type(vert), intent(in) :: vert_te(1:2000)
  type(flist), intent(in) :: flist_te(1:5000)
  real(8), intent(inout) :: surf
  real(8), intent(in) :: length(1:3), lengthh(1:3)
  integer :: i0, i1, link_index, m
  real(8) :: tempr(1:4), pos0(1:4), pos1(1:4)
  surf=0.0
  link_index=flist_te(face0.fptr).link
  m=0
  do while (  flist_te(link_index).link .ne. face0.fptr )
   pos0(1:3)=vert_te(flist_te(link_index).v).posv(1:3)-vert_te(flist_te(face0.fptr).v).posv(1:3)
   
   !call calc_distan(length, lengthh, pos0)
   !!write(*,*) "vertex: ", flist_te(face0.fptr).v, flist_te(link_index).v
   !!write(*,*) "pos0: ", pos0(1:3)
   link_index=flist_te(link_index).link
   pos1(1:3)=vert_te(flist_te(link_index).v).posv(1:3)-vert_te(flist_te(face0.fptr).v).posv(1:3)
   !call calc_distan(length, lengthh, pos1)
   !!write(*,*) "vertex: ", flist_te(face0.fptr).v, flist_te(link_index).v
   !!write(*,*) "pos1: ", pos1(1:3)
   call crossproduct( pos0(1:3), pos1(1:3), tempr(1:3))
   !!write(*,*) "tempr: ", tempr(1:3)
   call dotproduct( tempr(1:3), tempr(1:3), tempr(4))
   !m=m+1
   !write(*,*) "surface: ", m, (sqrt(tempr(4)))/(2.0)
   surf=surf+(sqrt(tempr(4)))/(2.0)
  enddo  
  end subroutine

!=====================hb function
!=======================
  subroutine ohpair(pos_O, pos_H, num_O, num_H, length, lengthh, calc_distan, oh_pair)
  implicit none
  real(8), intent(in) :: lengthh(1:3), length(1:3)
  real(8), dimension(3,num_O), intent(in) :: pos_O
  real(8), dimension(3,num_H), intent(in) :: pos_H
  integer :: i0, i1, m
  real(8) :: thr=1.25**2, pos(1:4)
  integer, intent(in) :: num_O, num_H
  integer, dimension(3,num_O), intent(inout) :: oh_pair
  oh_pair=0
  do i0=1, num_O
   m=1
   do i1=1, num_H
    pos(1:3)=pos_O(1:3,i0)-pos_H(1:3,i1)
    call calc_distan(length, lengthh, pos)
    if (pos(4) < thr) then
     oh_pair(m,i0)=i1
     m=m+1
    endif
   enddo
  enddo
  end subroutine ohpair  

  subroutine O_Ohb(pos_O1, pos_O2, pos_O1H, thetaOO, length, lengthh, yes_or_no)
  implicit none
  real(8), intent(in)  :: pos_O1(1:3), pos_O2(1:3), thetaOO, pos_O1H(1:3)
  real(8), intent(in)  :: length(1:3), lengthh(1:3)
  logical, intent(inout) :: yes_or_no
  real(8) :: temp(1:4), temp0(1:4), cosbeta, pi=3.1415926
  yes_or_no= .FALSE.
  temp(1:3)=pos_O2(1:3)-pos_O1(1:3)
  call calc_distan(length, lengthh, temp)
  if ( temp(4)< thetaOO ) then
      temp0(1:3)=pos_O1H(1:3)-pos_O1(1:3)
      call calc_distan(length, lengthh, temp0)
      cosbeta=(temp0(1)*temp(1)+temp0(2)*temp(2)+temp0(3)*temp(3))/sqrt(temp(4)*temp0(4))
   if ( acos(cosbeta) < pi/6 ) then
    yes_or_no= .TRUE.
   endif
  endif
  end subroutine O_Ohb

  subroutine find_acc(pos_O, pos_H, ind_atom, oh_pair, thetaoo, acc_index, length, lengthh, num_o, num_h)
  implicit none
  real(8), intent(in) :: length(1:3), lengthh(1:3), thetaoo;
  integer, intent(in) :: ind_atom, num_o, num_h
  real(8), dimension(3,num_o), intent(in) :: pos_O
  integer, dimension(3, num_o), intent(in) :: oh_pair
  real(8), dimension(3,num_h), intent(in) :: pos_H
  integer, intent(inout) :: acc_index(1:3);
  integer :: w1, w2, w3;
  real(8)  ::  tempw(1:4);
  logical :: yes;
  acc_index=0
  w3=1;
  do w2=1, 3
   if ( oh_pair(w2,ind_atom) .ne. 0) then
    do w1=1, num_o
     if ( w1 .ne. ind_atom ) then
      yes= .FALSE.
      call O_Ohb(pos_O(1:3,ind_atom), pos_O(1:3,w1), pos_H(1:3,oh_pair(w2,ind_atom)), thetaoo, length, lengthh, yes)
      if ( yes== .TRUE.) then
       acc_index(w3)=w1
       w3=w3+1
      endif
     endif
    enddo
   endif
  enddo
  end subroutine find_acc

  subroutine find_donor(pos_O, pos_H, ind_atom, oh_pair, thetaoo, donor_index, length, lengthh, num_o, num_h)
  implicit none
  real(8), intent(in) :: length(1:3), lengthh(1:3), thetaoo;
  integer, intent(in) :: ind_atom, num_o, num_h
  real(8), dimension(3,num_o), intent(in) :: pos_O
  integer, dimension(3, num_o), intent(in) :: oh_pair
  real(8), dimension(3,num_h), intent(in) :: pos_H
  integer, intent(inout) :: donor_index(1:3);
  integer :: w1, w2, w3;
  real(8)  ::  tempw(1:4);
  logical :: yes;
  donor_index=0
  w3=1;
  do w1=1, num_o
   if (w1 .ne. ind_atom) then
    do w2 =1, 3
     yes=.FALSE.
     if ( oh_pair(w2,w1) .ne. 0) then
      call O_Ohb(pos_O(1:3,w1), pos_O(1:3,ind_atom), pos_H(1:3, oh_pair(w2,w1)), thetaoo, length, lengthh, yes)
      if ( yes== .TRUE.) then
       donor_index(w3)=w1
       w3=w3+1
      endif
     endif
    enddo
   endif
  enddo
  end subroutine find_donor

  
end module
