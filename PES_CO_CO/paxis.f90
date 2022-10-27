!*******************************************************************************
! rotate to principal axis, take Cartesian geometry as input
! return new Cartesian geometry ( rotated to PA if ipaxis=1)
! return inertia matrix
  subroutine paxis(ipaxis, natom, xmass, x, inert)
  implicit real*8  (a-h,o-z) 
  integer, intent(in) :: ipaxis, natom
  real*8, intent(in) :: xmass(natom)
  real*8, intent(inout) :: x(3,natom)
  real*8, intent(out) :: inert(3,3)

  real*8 :: dum(3), dum2(3),rotmat(3,3),II(3,3)
  real*8 :: xcom(3),norm
  CHARACTER LABEL(3)
  DATA LABEL /'a','b','c'/
  data zero/0d0/

  nfout = 21
! ================================================
! 1) shift origin to com coordinates xcom
  xmassall = zero
  do iatom = 1, natom
     xmassall = xmassall + xmass(iatom)
  end do

  xcom = zero
  do iatom = 1, natom
     xcom(1) = xcom(1) + xmass(iatom)*x(1,iatom)
     xcom(2) = xcom(2) + xmass(iatom)*x(2,iatom)
     xcom(3) = xcom(3) + xmass(iatom)*x(3,iatom)
  end do
  xcom = xcom / xmassall
  !  write(nfout,'(a5,3f15.7)')'xcom',xcom
  !  write(nfout,*)
  
  ! shift Cartesians of each atom with origin at com
  do iatom = 1, natom
     do i = 1, 3
        x(i,iatom) = x(i,iatom) - xcom(i)
     end do
  end do
!  write(322,*) xcom
!  write(322,*)
!  write(nfout,*)'reference Cartesian (au) after moving to com'
!  do iatom = 1, natom
!    write(nfout,'(i5,3f15.7)')iatom, x(1:3,iatom)
!  end do
!  write(nfout,*)

! ================================================
! 2) convert to principal axes
! 1. diagonalise inertia tensor   (borrowing rotmat)
!
  nfout = 21
  inert(1:3,1:3) = zero
  nat = natom
  do i=1,nat
     inert(1,1) = inert(1,1) + xmass(i)*(x(2,i)**2 + x(3,i)**2)
     inert(2,2) = inert(2,2) + xmass(i)*(x(1,i)**2 + x(3,i)**2)
     inert(3,3) = inert(3,3) + xmass(i)*(x(1,i)**2 + x(2,i)**2)
     inert(1,2) = inert(1,2) - xmass(i)*x(1,i)*x(2,i)
     inert(1,3) = inert(1,3) - xmass(i)*x(1,i)*x(3,i)
     inert(2,3) = inert(2,3) - xmass(i)*x(2,i)*x(3,i)
  end do
  inert(2,1) = inert(1,2)
  inert(3,1) = inert(1,3)
  inert(3,2) = inert(2,3)
  
  if (ipaxis == 0) then 
     II=0d0
     do i=1,3
        II(i,i)=1d0
     enddo
!     do i=1,nat
!        do j=1,3
!           rotvecs((i-1)*3+j,1)=(x(2,i)*II(j,3)-x(3,i)*II(j,2))*sqrt(xmass(i))
!           rotvecs((i-1)*3+j,2)=(x(3,i)*II(j,1)-x(1,i)*II(j,3))*sqrt(xmass(i))
!           rotvecs((i-1)*3+j,3)=(x(1,i)*II(j,2)-x(2,i)*II(j,1))*sqrt(xmass(i))
!        enddo
!     enddo
 !    do j=1,3
 !       norm=0d0
 !       do i=1,3*nat
 !          norm=norm+(rotvecs(i,j))**2
 !      enddo
 !       norm=sqrt(norm)
 !       do i=1,3*nat
 !          rotvecs(i,j)=rotvecs(i,j)/norm
 !       enddo
 !    enddo

     !    write(nfout,*)'No rotation to PA' 
  else if (ipaxis == 1) then
     !    write(nfout,*)'Moment of inertia before rotating to PA' 
     !    write(nfout,'(3f15.7)')inert 
     
     
     rotmat = inert
     
     ntmp =3
     call diag(ntmp, rotmat,ntmp, dum, dum2)
!!$    inert=0d0
!!$    do i=1,3
!!$       inert(i,i)=dum(i)
!!$    enddo

!    write(nfout,*)'Diagonal Moment of inertia in PA' 
!    write(nfout,*) dum(:)
!    write(nfout,*)'rotation matrix to rotate to PA' 
!    write(nfout,'(3f15.7)')rotmat 
!
!   2.rotate cartesians:  x = rotmat (3x3) * x (3xnat)
!   (note that the new x,y,z are inertial axes a,b,c in that order)
!
!   check sign of determinant to establish whether handedness of axis
!   system not changed, take remedial action if necessary
     rotdet= rotmat(1,1)*(rotmat(2,2)*rotmat(3,3)-rotmat(3,2)*rotmat(2,3)) &
          -rotmat(1,2)*(rotmat(2,1)*rotmat(3,3)-rotmat(3,1)*rotmat(2,3)) &
          +rotmat(1,3)*(rotmat(2,1)*rotmat(3,2)-rotmat(3,1)*rotmat(2,2))
     
     if(rotdet < 0.d0)then
        !      write(nfout,'(a)')'sign of second row of rotation matrix' &
        !                        //' reversed to preserve handedness'
        rotmat(1,2)=-rotmat(1,2)
        rotmat(2,2)=-rotmat(2,2)
        rotmat(3,2)=-rotmat(3,2)
     endif

!     do i=1,nat
!        do j=1,3
!           rotvecs((i-1)*3+j,1)=(x(2,i)*rotmat(j,3)-x(3,i)*rotmat(j,2))*sqrt(xmass(i))
!!           rotvecs((i-1)*3+j,2)=(x(3,i)*rotmat(j,1)-x(1,i)*rotmat(j,3))*sqrt(xmass(i))
 !          rotvecs((i-1)*3+j,3)=(x(1,i)*rotmat(j,2)-x(2,i)*rotmat(j,1))*sqrt(xmass(i))
 !       enddo
 !    enddo
!     do j=1,3
!        norm=0d0
!        do i=1,3*nat
!           norm=norm+(rotvecs(i,j))**2
!        enddo
!        norm=sqrt(norm)
!        do i=1,3*nat
!           rotvecs(i,j)=rotvecs(i,j)/norm
!        enddo
!    enddo
     !
!   read external rotation matrix if specified - this is to be selected so
!   that handedness of coordinates is not changed, ie determinant is +1.
!   note that for c3v molecules conventionally the c3 (z) axis is with the
!   threefold group at the negative end, the xz plane is a symmetry plane,
!   and the off-axis atom on that plane is for -ve x
!
     inprot = 0
     if(inprot.eq.1)then
       !      write(nfout,'(1x/'' note: external rotation matrix -->'')')
        do j=1,3
           read(ir,*)(rotmat(i,j),i=1,3)
        end do
     endif
!
!   rotate
     do i=1,nat
        xtemp=zero
        ytemp=zero
        ztemp=zero
        do j=1,3
           xtemp = xtemp+rotmat(j,1)*x(j,i)
           ytemp = ytemp+rotmat(j,2)*x(j,i)
           ztemp = ztemp+rotmat(j,3)*x(j,i)
        end do
        x(1,i)= xtemp
        x(2,i)= ytemp
        x(3,i)= ztemp
     end do
     write(322,*) rotmat(1,1:3)
     write(322,*) rotmat(2,1:3)
     write(322,*) rotmat(3,1:3)
!    write(nfout,206)
     do j=1,3
        !      write(nfout,207)label(j),(rotmat(i,j),i=1,3)
     end do
     rotdet= rotmat(1,1)*(rotmat(2,2)*rotmat(3,3)-rotmat(3,2)*rotmat(2,3)) &
          -rotmat(1,2)*(rotmat(2,1)*rotmat(3,3)-rotmat(3,1)*rotmat(2,3)) &
          +rotmat(1,3)*(rotmat(2,1)*rotmat(3,2)-rotmat(3,1)*rotmat(2,2))
     
     !    write(nfout,227)rotdet
    227 format(' determinant=',f9.6)
    206 format(/' rotation matrix to principal coordinates:'/ &
           ' -----------------------------------------'      &
        /16x,'xold',16x,'yold',16x,'zold'/)
    207 format(1x,a,3(f20.6))
   
     !    write(nfout,205)
     do k=1,nat
        !      write(nfout,204)k,(x(kk,k),kk=1,3),xmass(k)
     end do
    101 format(1x,a4,2i5,3f10.4,2i5,f10.4)
    201 format(/' atomic coordinate input (angstr):'/  &
           ' ---------------------------------'/  &
    5x,'atom',2x,'origin',2x,  &
    'polar',3x,'r',9x,'theta',9x,'phi',5x,'ll',4x,'nax',4x,'alfa'/ &
    20x,'axis')
    202 format(1x,i2,2x,a4,i5,4x,a1,f12.6,f12.6,f12.6,2i5,f12.6)
    203 format(/' cartesian coordinates (angstr):'/  &
               ' -------------------------------'   &
                /' atom no.',9x,'x',13x,'y',13x,'z',7x,'mass (amu)'/)
    205 format(/' principal coordinates (angstr):'/ &
              ' ------------------------------'  &
              /' atom no.',9x,'a',13x,'b',13x,'c',7x,'mass (amu)'/)
    204 format(1x,i2,7x,3(f12.6,2x),f12.6)
  end if ! ipaxis


  nfout = 21
  inert(1:3,1:3) = zero
  nat = natom
  do i=1,nat
     inert(1,1) = inert(1,1) + xmass(i)*(x(2,i)**2 + x(3,i)**2)
     inert(2,2) = inert(2,2) + xmass(i)*(x(1,i)**2 + x(3,i)**2)
     inert(3,3) = inert(3,3) + xmass(i)*(x(1,i)**2 + x(2,i)**2)
     inert(1,2) = inert(1,2) - xmass(i)*x(1,i)*x(2,i)
     inert(1,3) = inert(1,3) - xmass(i)*x(1,i)*x(3,i)
     inert(2,3) = inert(2,3) - xmass(i)*x(2,i)*x(3,i)
  end do
  inert(2,1) = inert(1,2)
  inert(3,1) = inert(1,3)
  inert(3,2) = inert(2,3)
  
  

  end subroutine paxis
