module te
use cut
implicit none
contains 
subroutine initialization2(length, lengthh, poso, index_origin, order, num_order, near_neighbor, vol, tot_surf)
 type(vert) :: vert_te(1:2000)
 type(edge) :: edge_te(1:2000)
 type(face) :: face_te(1:50) 
 type(flist) :: flist_te(1:5000)
 real(8) :: r2, r6, rangelim, w
! the length of link id max_flist
 integer, parameter :: max_flist=500
 integer :: m, n, ne, nf, nv, s, elast=6, vlast=4, flast=4, i0, i1, flistlast, index_order
 integer :: list_v_e(1:100), list_e_v(1:100)
 integer :: list_e_f(1:100), list_flist_e(1:100), list_flist_v(1:100)
 integer :: list_f_vfar(1:100), link, p0, p1
 real(8), intent(in) :: length(1:3), lengthh(1:3), poso(1:3, 1:1024)
 integer, intent(in) :: index_origin, num_order
 real(8), intent(inout) ::  vol, tot_surf
 integer, intent(inout) :: near_neighbor(1:num_order)
 real(8) :: origin(1:3), second(1:3)
 integer, intent(in) :: order(1:num_order)
 logical :: near
 integer, allocatable :: possible_neighbor(:)
 allocate(possible_neighbor(num_order))
 !================define the initial large regular vertex for the voronoi cell
 rangelim=31.0
 elast=6
 vlast=4 
 flast=4
 list_v_e(1:12)=[1, 3, 6, 1, 2, 5, 2, 3, 4, 4, 5, 6]
 list_e_v(1:12)=[1,2, 2, 3, 1, 3, 3, 4, 2, 4, 1,4]
 list_e_f(1:12)=[1, 4, 1, 2, 1, 3, 2, 3, 2, 4, 3, 4]
 list_flist_e(1:12)=[1, 2, 3, 2, 4, 5, 3, 4, 6, 6, 5, 1]
 list_flist_v(1:12)=[1, 2, 3, 2, 3, 4, 1, 3, 4, 1, 4, 2 ]
 list_f_vfar(1:4)=[1, 2, 1, 1]
 origin(1:3)=poso(1:3, index_origin)

!=====define vert.posv
 vert_te(1).posv(1:3)=[(sqrt(6.0)*rangelim)/(3.0), (-sqrt(2.0)*rangelim/(3.0)), -(rangelim/(3.0))]+origin(1:3)
 vert_te(2).posv(1:3)=[-(sqrt(6.0)*rangelim)/(3.0), (-sqrt(2.0)*rangelim/(3.0)), -(rangelim/(3.0))]+origin(1:3)
 vert_te(3).posv(1:3)=[0.0*rangelim, (2*sqrt(2.0)*rangelim)/(3.0), -(rangelim/(3.0))]+origin(1:3)
 vert_te(4).posv(1:3)=[0.0*rangelim, 0.0*rangelim, rangelim]+origin(1:3)
!=====square the rangelim 
 vert_te(1:4).posv(4)=sqrt(rangelim) 
 vert_te(1:4).stat=2
!=====define e
 m=0
 do i0=1, vlast
  do i1=1, 3
   m=m+1
   vert_te(i0).e(i1)=list_v_e(m)
  enddo
 enddo
 m=0

!=====define   edge's vertex, edge's face
 m=0
 do i0=1, elast
  do i1=1, 2
   m=m+1
   edge_te(i0).v(i1)=list_e_v(m)
  enddo
 enddo
 m=0

 m=0
 do i0=1, elast
  do i1=1, 2
   m=m+1
   edge_te(i0).f(i1)=list_e_f(m)
  enddo
 enddo
 m=0
 edge_te(1:elast).stat=3
!=============

!=========define link, face
 do i0=1, max_flist
  flist_te(i0).link=i0+1
 enddo

 m=0 
 do i0=1, flast
  face_te(i0).stat=3
  face_te(i0).fptr=m+1
  face_te(i0).vfar=list_f_vfar(i0)
  do i1=1, 3
   m=m+1
   flist_te(m).e=list_flist_e(m)
   flist_te(m).v=list_flist_v(m)
  enddo
  flist_te(m).link=face_te(i0).fptr
 enddo
 flistlast=m
 m=0
!===== face distance not squared
 face_te(1:flast).dist=rangelim/(3.0)


 !! can be rewritten
 p0=0
 do index_order =1, num_order
  !!write(*,*) "index_order: ", index_order
  second(1:3)=poso(1:3, order(index_order))
  !!write(*,*) "second: ", second(1:3)
  near=.false.
  call bisectplane( origin, second, length, lengthh, vert_te, edge_te, face_te , flist_te, elast, vlast, flast, flistlast, near)
  if ( near .eq. .true. ) then
   p0=p0+1
   possible_neighbor(p0)=order(index_order)
   !!near_neighbor(p0)=order(index_order)
  endif
 enddo


 p1=0 
 do i0=1, flast
  if ( face_te(i0).stat .eq. 3 ) then
   link=face_te(i0).fptr
   link=flist_te(link).link
   do while ( link .ne. face_te(i0).fptr ) 
    link=flist_te(link).link
   enddo
   p1=p1+1
   near_neighbor(p1)=possible_neighbor(i0-4)
  endif
 enddo
 call volume(face_te, flist_te, vert_te, flast, vol, tot_surf, length, lengthh)
 deallocate(possible_neighbor)
end subroutine

 !======determine whether the VC is closed given by specfic set of neighbors
subroutine initialization(length, lengthh, poso, index_origin, order, num_order, vol, tot_surf, closure)
 type(vert) :: vert_te(1:2000)
 type(edge) :: edge_te(1:2000)
 type(face) :: face_te(1:50) 
 type(flist) :: flist_te(1:5000)
 real(8) :: r2, r6, rangelim, w
! the length of link id max_flist
 integer, parameter :: max_flist=500
 integer :: m, n, ne, nf, nv, s, elast=6, vlast=4, flast=4, i0, i1, flistlast, index_order
 integer :: list_v_e(1:100), list_e_v(1:100)
 integer :: list_e_f(1:100), list_flist_e(1:100), list_flist_v(1:100)
 integer :: list_f_vfar(1:100), link
 real(8), intent(in) :: length(1:3), lengthh(1:3), poso(1:3, 1:1024)
 integer, intent(in) :: index_origin, num_order
 real(8), intent(inout) ::  vol, tot_surf
 real(8) :: origin(1:3), second(1:3)
 integer, intent(in) :: order(1:num_order)
 logical :: near
 logical, intent(inout) :: closure
 rangelim=31.0
 elast=6
 vlast=4 
 flast=4
 closure=.true.
 list_v_e(1:12)=[1, 3, 6, 1, 2, 5, 2, 3, 4, 4, 5, 6]
 list_e_v(1:12)=[1,2, 2, 3, 1, 3, 3, 4, 2, 4, 1,4]
 list_e_f(1:12)=[1, 4, 1, 2, 1, 3, 2, 3, 2, 4, 3, 4]
 list_flist_e(1:12)=[1, 2, 3, 2, 4, 5, 3, 4, 6, 6, 5, 1]
 list_flist_v(1:12)=[1, 2, 3, 2, 3, 4, 1, 3, 4, 1, 4, 2 ]
 list_f_vfar(1:4)=[1, 2, 1, 1]
 origin(1:3)=poso(1:3, index_origin)

!=====define vert.posv
 vert_te(1).posv(1:3)=[(sqrt(6.0)*rangelim)/(3.0), (-sqrt(2.0)*rangelim/(3.0)), -(rangelim/(3.0))]+origin(1:3)
 vert_te(2).posv(1:3)=[-(sqrt(6.0)*rangelim)/(3.0), (-sqrt(2.0)*rangelim/(3.0)), -(rangelim/(3.0))]+origin(1:3)
 vert_te(3).posv(1:3)=[0.0*rangelim, (2*sqrt(2.0)*rangelim)/(3.0), -(rangelim/(3.0))]+origin(1:3)
 vert_te(4).posv(1:3)=[0.0*rangelim, 0.0*rangelim, rangelim]+origin(1:3)
!=====square the rangelim 
 vert_te(1:4).posv(4)=sqrt(rangelim) 
 vert_te(1:4).stat=2
!=====define e
 m=0
 do i0=1, vlast
  do i1=1, 3
   m=m+1
   vert_te(i0).e(i1)=list_v_e(m)
  enddo
 enddo
 m=0

!=====define   edge's vertex, edge's face
 m=0
 do i0=1, elast
  do i1=1, 2
   m=m+1
   edge_te(i0).v(i1)=list_e_v(m)
  enddo
 enddo
 m=0

 m=0
 do i0=1, elast
  do i1=1, 2
   m=m+1
   edge_te(i0).f(i1)=list_e_f(m)
  enddo
 enddo
 m=0
 edge_te(1:elast).stat=3
!=============

!=========define link, face
 do i0=1, max_flist
  flist_te(i0).link=i0+1
 enddo

 m=0 
 do i0=1, flast
  face_te(i0).stat=3
  face_te(i0).fptr=m+1
  face_te(i0).vfar=list_f_vfar(i0)
  do i1=1, 3
   m=m+1
   flist_te(m).e=list_flist_e(m)
   flist_te(m).v=list_flist_v(m)
  enddo
  flist_te(m).link=face_te(i0).fptr
 enddo
 flistlast=m
 m=0
!===== face distance not squared
 face_te(1:flast).dist=rangelim/(3.0)



 !! can be rewritten
 do index_order =1, num_order
  second(1:3)=poso(1:3, order(index_order))
  call bisectplane( origin, second, length, lengthh, vert_te, edge_te, face_te , flist_te, elast, vlast, flast, flistlast, near)
 enddo
 do i0=1, flast
  if ( face_te(i0).stat .eq. 3 ) then
   link=face_te(i0).fptr
   link=flist_te(link).link
   do while ( link .ne. face_te(i0).fptr ) 
    link=flist_te(link).link
   enddo
  endif
 enddo
 do i0=1, 4
  if ( face_te(i0).stat .eq. 3 ) then
   closure= .false.
  endif
 enddo
 call volume(face_te, flist_te, vert_te, flast, vol, tot_surf, length, lengthh)
end subroutine

end module te
