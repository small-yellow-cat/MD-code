module function
implicit none
public calc_distan
contains
!=========calculate distance
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
!==============
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

 ! deternmine whether o1, h1, o2 form HBs
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
      call O_Ohb(pos_O(1:3,ind_atom), pos_O(1:3,w1), pos_H(1:3, oh_pair(w2,ind_atom)), thetaoo, length, lengthh, yes)
      if ( yes== .TRUE.) then
       acc_index(w3)=w1
       w3=w3+1
      endif
     endif
    enddo
   endif
  enddo
  end subroutine find_acc





! find HB donor
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
 
  !======== vector operation   3d vector involved
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

  recursive subroutine quicksort(a, first, last, num)
  implicit none
  integer, intent(in) :: first, last, num
  real(8), intent(inout) ::  a(1:num)
  real(8) :: x, t
  integer :: i, j

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
      i=i+1
      j=j-1
  enddo
  if (first < i-1) call quicksort(a, first, i-1, num)
  if (j+1 < last)  call quicksort(a, j+1, last, num)
  end subroutine quicksort
  
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

  !==========Qindex
  subroutine Qindex( near_O, length, lengthh, Qi)
  implicit none
  real(8), intent(in) :: near_O(1:4,1:4), length(1:3), lengthh(1:3)
  real(8), intent(out) :: Qi
  integer :: i0, i1
  real(8) :: cosbeta
  Qi=0
  do i0=1, 3
   do i1=i0+1, 4
    cosbeta=(near_O(1,i0)*near_O(1,i1)+near_O(2,i0)*near_O(2,i1)+near_O(3,i0)*near_O(3,i1))/sqrt(near_O(4,i0)*near_O(4,i1))
    Qi=Qi+(cosbeta+1.d0/3)**2
   enddo
  enddo
  Qi=(3.0*Qi)/32.0
  end subroutine Qindex

  !=========sk index
  subroutine Sindex( near_O, length, lengthh, Sk)
  implicit none
  real(8), intent(in) :: near_O(1:4,1:4), length(1:3), lengthh(1:3)
  real(8), intent(out) :: Sk
  integer :: i0, i1
  real(8) :: ave_r, sum0
  ave_r=0.0
  sum0=0.0
  do i0=1, 4
   ave_r=ave_r+sqrt(near_O(4,i0))
  enddo
  ave_r=ave_r/(4.0)
  Sk=0
  do i0=1, 4
   sum0=sum0+(sqrt(near_O(4,i0))-ave_r)**2
  enddo
  Sk=1-sum0/(3*4*(ave_r)*(ave_r))
  end subroutine Sindex
 
  !====sk index 2
  subroutine Sindex2( near_O, length, lengthh, Sk)
  implicit none
  real(8), intent(in) :: near_O(1:4,1:4), length(1:3), lengthh(1:3)
  real(8), intent(out) :: Sk
  integer :: i0, i1
  real(8) :: ave_r, sum0
  ave_r=0.0
  sum0=0.0
  do i0=1, 4
   ave_r=ave_r+sqrt(near_O(4,i0))
  enddo
  ave_r=ave_r/(4.0)
  Sk=0
  do i0=1, 4
   sum0=sum0+(sqrt(near_O(4,i0))-ave_r)**2
  enddo
  Sk=sum0/(3*4*(ave_r)*(ave_r))
  end subroutine Sindex2 
  !===lsi
  subroutine Lindex(near_o, num_near_o, length, lengthh, Lsi)
  implicit none
  real(8), intent(in) :: length(1:3), lengthh(1:3)
  integer, intent(in) :: num_near_o
  real(8), dimension(4, num_near_o), intent(in) :: near_o
  real(8), intent(out) :: Lsi
  integer :: i0, i1, time
  real(8) :: ave_r, sum0   
  ave_r=0.0
  sum0=0.0
  i0=1
  time=0
  do while ( near_o(4,i0) < (3.7)**2)
   ave_r=ave_r+sqrt(near_O(4,i0+1))-sqrt(near_O(4,i0))
   time=time+1
   i0=i0+1
  enddo
  ave_r=ave_r/(time+0.00)
  Lsi=0
  do i0=1, time
   sum0=sum0+(sqrt(near_O(4,i0+1))-sqrt(near_o(4,i0))-ave_r)**2
  enddo
  Lsi=(sum0+0.0)/(time+0.0)
  end subroutine Lindex

end module
