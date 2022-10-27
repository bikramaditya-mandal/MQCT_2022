program CO
  implicit none
  integer :: i,j,k,flag,natom,i2,i3,list(4)
  real*8 :: tampon,h2wn,V,pii,poten,poten2,qmat(3,6),inert(3,3),bohr,range,dx,coord1(200),coord2(100),temp,temp3,temp2(5),position(4)
  real*8 :: jac(4),jac2(4),djac(4),shift,hartokcl,atobohr,ugrad,ass,mass(6),ref1(6),ref2(6),cart3(12),glob_min,EGrid(100,200),EGrid2(100,50)
  character(len=3),allocatable :: symb(:)

  pii=acos(-1d0)
  h2wn = 219474.63d0
  hartokcl=627.5095d0
  atobohr=1.0d0/0.529177249d0 
  ugrad=atobohr*hartokcl
  bohr=0.529177249d0
  natom=4
  allocate(symb(natom))
   open(unit=10,file='frag_data')
   do i=1,natom
      read(10,*) symb(i)    ! element label atom 1: eg. 'H'
   enddo
   do i=1,natom
      read(10,*) mass(i)
      mass(i)=mass(i)*1837.152697d0/1.007825032d0
   enddo
   do i=1,6
      read(10,*) ref1(i)
   enddo
   do i=1,6
      read(10,*) ref2(i)
   enddo
   close(10)
!!Fitted range is 15d0 > R > 2.64d0 (dist between centers of mass in Angstroms), Jac(2) and Jac(3) are cos(theta1) and cos(theta2) and range from (-1,1). Jac(4) is the dihedral with range (0,pi).
!!! PES returns potential "V".
  !!! Energy is in hartrees.


!!!!! set asymptote energy
  jac(1)=15d0
  jac(2)=0d0
  jac(3)=0d0
  jac(4)=pii/2d0
  call IMLS(jac,V)
  ass=V
  write(*,*) ass

!test side-by-side scan
! to write cartesian coords to fort.78, uncomment lines 163-168 in "PES_4D_subroutine_vdW.f90"
! atom order is O1,C1,O2,C2

jac(2)=0d0
jac(3)=0d0
jac(4)=0d0


do i=1,20
jac(1)=3.1d0+dble(i-1)*3d0/19d0
call IMLS(jac,V)
write(77,*) jac(1),(v-ass)*h2wn
write(79,*) jac(:)
write(79,*)

enddo


end program CO

