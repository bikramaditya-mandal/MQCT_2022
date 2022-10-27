      SUBROUTINE USER_DEFINED_PES(V,R,rvib,rvib2,alpha,beta,gamma,aalpha,bbeta,ggamma)
!     INPUT:  R - distance betweenn COMs, rvib - vibrational coordinate
! 		      alpha,beta,gamma - Euler's angles of the first molecule
!   		  aalpha,bbeta, ggamma - Euler's angles of the second molecule
!     OUTPUT: V - value of the potential
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 V,R,alpha,beta,gamma,aalpha,bbeta,ggamma,t,rvib,rvib2
      REAL*8 R1,COG, phi
	  real*8, parameter :: PI = 4.d0*datan(1.d0)
	  integer :: i,j,k,flag,natom,i2,i3,list(4)
	  real*8 :: tampon,h2wn,pii,poten,poten2,qmat(3,6),inert(3,3),bohr,range,dx
	  real*8 :: coord1(200),coord2(100),temp,temp3,temp2(5),position(4)
	  real*8 :: jac(4),jac2(4),djac(4),shift,hartokcl,atobohr,ugrad,ass,mass(6)
	  real*8 :: ref1(6),ref2(6),cart3(12),glob_min,EGrid(100,200),EGrid2(100,50)
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
!	  write(*,*) ass
	  
	  !test side-by-side scan
	  ! to write cartesian coords to fort.78, uncomment lines 163-168 in "PES_4D_subroutine_vdW.f90"
	  ! atom order is O1,C1,O2,C2
	  
	  jac(1)= R
	  jac(2)= dcos(beta)
	  jac(3)= dcos(bbeta)
	  jac(4)= PI - gamma
	  call IMLS(jac,V)
	  V = V-ass
    
      END SUBROUTINE USER_DEFINED_PES


!----------!
! OPTION 2 !
!----------! 
!	  USE KEYWORD "EXPANSION=YES" TO INITIATE THIS OPTION
      SUBROUTINE USER_DEFINED_TERMS(T,I,R)
!     THIS SUBROTUNE COMPUTES RADIAL COEFFICENTS OF THE PES EXPANSION AT A GIVEN DISTANCE
!     INPUT:  R - distance between COMs of particles, I - TERM NUMBER
!     OUTPUT: T - value of coefficent 	  
      IMPLICIT NONE	  
      REAL*8 T,R
      INTEGER I
!     USER MUST INSERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:
      STOP "ERROR: USER_DEFINED_TERMS IS NOT SUPPLIED"
	  END SUBROUTINE USER_DEFINED_TERMS
!----------!
! OPTION 3 !
!----------! 
!	  USE KEYWORDS "EXPANSION=YES, TERMS_FILE=YES" TO INITIATE THIS OPTION
! 	  SIMILAR TO OPTION 2, BUT NO SUBROUTINE IS REQUIRED
!     USER SHOULD PROVIDE THE FILE EXPAN_PES_TERMS.DAT 
!     IN THE MAIN PROGRAM DIRECTORY CONTAINING THE COEFFICEINS 
!     OF POTENTIAL EXPANSION PRECOMPUTED EXTERNALLY.
! 	  SEE EXAMPLE FILES SUPPLIED WITH THE CODE.
!----------!
! OPTION 4 !
!----------! 
!	  USE KEYWORDS "EXPANSION=YES, TERMS_ONFLY=YES" TO INITIATE THIS OPTION
      SUBROUTINE USER_DEFINED_COEFFS(T,DTDR,I,R) 
!     THIS SUBROUTINE COMPUTES RADIAL COEFFICENTS OF THE PES EXPANSION 
!     AND THEIR DERIVATIVES AT A GIVEN DISTANCE R
      IMPLICIT NONE
!     INPUT : R - distance between COMs of particles, I - TERM NUMBER
!     OUTPUT: T - value of coefficent, DTDR - its radial derivative 	  
      REAL*8 T,R,DTDR 
      INTEGER I
!     USER MUST INCERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE IS SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:	
      STOP "ERROR: USER_DEFINED_COEFFS IS NOT SUPPLIED"
      END SUBROUTINE USER_DEFINED_COEFFS 
