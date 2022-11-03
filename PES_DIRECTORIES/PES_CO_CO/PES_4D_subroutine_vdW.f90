subroutine IMLS(jac3,V)  
use dynamic_parameters

implicit none
  integer :: i,j,k,initflag,ierr
  real*8 :: xi(4),tampon
  real*8 :: temp(5),temp3,jac3(4),djac3(4),cart3(12)
  real*8 :: somme,somme2
  real*8 :: h2wn,v,pii
  save initflag
  data initflag /1/
INTERFACE
   FUNCTION func_actual(xi)    !!!energy of largest basis and high-level ab initio
     USE dynamic_parameters
     IMPLICIT NONE
     REAL*8, DIMENSION(:), INTENT(IN) :: xi
     REAL*8 :: func_actual
   END FUNCTION func_actual
end interface
INTERFACE
   FUNCTION func_actual_seed(xi)  !!!energy of minimal basis and low-level ab initio
     USE dynamic_parameters
     IMPLICIT NONE
     REAL*8, DIMENSION(:), INTENT(IN) :: xi
     REAL*8 :: func_actual_seed
   END FUNCTION func_actual_seed
end interface
INTERFACE
   FUNCTION dfunc_actual_anal1(xi)    !!!energy and analytic gradient of largest basis and high-level ab initio
     USE dynamic_parameters
     IMPLICIT NONE
     REAL*8, DIMENSION(:), INTENT(IN) :: xi
     REAL*8,dimension(size(xi)+1) :: dfunc_actual_anal1
   END FUNCTION dfunc_actual_anal1
end interface
hartokcl=627.5095d0
alpha=-1d0  !!!!coefficient in R=exp(alpha*r) coordinate
order_1_min=3  !!!!!basis for "minimal fit" to high level data, also used to fit low level grid.
order_2_min=3
order_3_min=3
order_4_min=3

pii=dacos(-1d0)

if(initflag==1) then   
   ddd=4
   epss=1d-14
   zz=4
   zz_low=4
   zz4=20
   natom1=2
   natom2=2
   natom=natom1+natom2
   nbdist=natom*(natom-1)/2
   W_a=0.7d0
   dist_tol=1.3d0
   allocate(ref1(3*natom1),ref2(3*natom2),bdist(nbdist))
   allocate(jac(4),jac2(4),cart(3*(natom)))
   allocate(symb(natom),mass(natom))
   open(unit=10,file='frag_data')
   do i=1,natom
      read(10,*) symb(i)    ! element label atom 1: eg. 'H'
   enddo
   do i=1,natom
      read(10,*) mass(i)
   enddo
   do i=1,3*natom1
      read(10,*) ref1(i)
   enddo
   do i=1,3*natom2
      read(10,*) ref2(i)
   enddo
   close(10)
   OPEN (UNIT=652, FILE='PES_CO_CBSsep_core.dat', FORM='UNFORMATTED', ACCESS='SEQUENTIAL') 
   read(652) count3
!   write(*,*) count3
   read(652) order_1
   read(652) order_2
   read(652) order_3
   read(652) order_4
   read(652) maxpoints
   read(652) mass
   read(652) rmax
   read(652) rmin
   read(652) Max_E
   read(652) low_grid
   read(652) count_seed
   call basis_size(ddd,order_1,order_2,order_3,order_4,basis_1)
   call basis_size(ddd,order_1-1,order_2-1,order_3-1,order_4-1,basis_2)  
   call basis_size(ddd,order_1_min,order_2_min,order_3_min,order_4_min,basis_3) 

   allocate(b2(basis_1,count3),b2_lower(basis_2,count3),b2_minimal(basis_3,count3),d(count3),coords(count3,4)) 
   allocate(b2_seed(basis_3,count_seed),d_seed(count_seed),coords_seed(count_seed,4))
   b2=0d0
   b2_lower=0d0
   b2_minimal=0d0
   b2_seed=0d0
   d=0d0
   d_seed=0d0
   coords=0d0
   coords_seed=0d0
   do i=1,count3 
      read(652) b2(:,i)       
   enddo
   do i=1,count3      
      read(652) b2_lower(:,i)       
   enddo
   do i=1,count3      
      read(652) b2_minimal(:,i)       
   enddo  
   do i=1,count3       
      read(652) d(i) 
   enddo
   do i=1,count3        
      read(652) coords(i,:)
   enddo
   if(low_grid>0)then
      read(652) Max_E_seed
      do i=1,count_seed      
         read(652) b2_seed(:,i)       
      enddo
      do i=1,count_seed 
         read(652) d_seed(i)       
      enddo
      do i=1,count_seed        
         read(652) coords_seed(i,:)       
      enddo
   endif
   
   close(652)
   initflag=2
endif
!!!!!!!!!


xi=jac3
if(xi(1)<rmin(1))then
   v=Max_E/hartokcl
   write(*,*) 'coords outside fitted range'
   write(*,*)  xi(1),rmax(1),rmin(1)
   return
endif
if(xi(4)>pii)then
   xi(4)=xi(4)-2d0*pii
endif
if(xi(4)<-pii)then
   xi(4)=xi(4)+2d0*pii
endif
!xi(4)=abs(xi(4))
if(xi(2)>1d0)then
   xi(2)=2d0-xi(2)
endif
if(xi(3)>1d0)then
   xi(3)=2d0-xi(3)
endif
if(xi(2)<-1d0)then
   xi(2)=-2d0-xi(2)
endif
if(xi(3)<-1d0)then
   xi(3)=-2d0-xi(3)
endif

!call INT_Cart(cart3,xi,mass,natom1,natom2,ref1,ref2)
!write(78,*) cart3(1:3)
!write(78,*) cart3(4:6)
!write(78,*) cart3(7:9)
!write(78,*) cart3(10:12)
!write(78,*)


temp3=func_actual_seed(xi)
if(temp3>Max_E_seed)then
   v=Max_E/hartokcl
   djac3=0d0
!   write(*,*) 'hit ceiling on low grid'
   return
endif
temp3=func_actual(xi)
if(temp3>Max_E)then
   temp3=Max_E
!   write(*,*) 'hit ceiling'
endif
V=temp3/hartokcl
return

temp=dfunc_actual_anal1(xi)


if(temp(1)>Max_E)then
   temp(1)=Max_E
   write(*,*) 'hit ceiling'
endif
V=temp(1)/hartokcl
djac3(1:4)=temp(2:5)/hartokcl

  return
end subroutine IMLS


subroutine cart_to_bdist_inter(x,natom1,natom2,dist_tol,flag)
  implicit none
  integer :: i,j,k,flag,natom1,natom2
  real*8 :: x(3*(natom1+natom2)),summ,dist_tol
  flag=0
  do i=1,natom1
     do j=natom1+1,natom1+natom2
        summ=0d0
        do k=1,3
           summ=summ+(x(3*(i-1)+k)-x(3*(j-1)+k))**2
        enddo
        if(sqrt(summ)<dist_tol)then
           flag=1
        endif
     enddo
  enddo
  return
end subroutine cart_to_bdist_inter
subroutine basis_size(d,order_1,order_2,order_3,order_4,basis)
  
  implicit none
  
  integer :: d,order_1,order_2,order_3,order_4,count,basis,l1,l2,m
  
  
!!!!!!!!basis calc
  count=0
  !  count=d*order_1
  
  do l1=0,order_2
     do l2=0,order_3
        if((l1+l2)<order_4+1)then
           do m=0,min(l1,l2)
              count=count+1
           enddo
        endif
     enddo
  enddo
  basis=count*(order_1)+1
!  basis=count+1
  return
end subroutine basis_size

subroutine INT_Cart(cart,internal2,mass,natom1,natom2,ref1,ref2)
  implicit none
  integer :: i,j,k,kp,lab,ierr,natom1,natom2
  real*8 :: internal(6),internal2(4),cart((natom1+natom2)*3),mass(natom1+natom2),ref1(natom1*3),ref1_temp(natom1*3),ref2(natom2*3),ref2_temp(natom2*3)
  real*8 :: cart_mat1(3,natom1),cart_mat(3,natom1+natom2),cart_ref1(3,natom1),cart_ref2(3,natom2),U_rot(3,3),cm(3)
  real*8 :: cart_mat2(3,natom2),cart_frag2(natom2*3),quat(4),quat2(4),pii,vec1(3),vec2(3)
  real*8 :: gamma1,gamma2,beta1,beta2,alpha1,alpha2
  real*8 :: cart_ref1b(3,natom1),cart_ref2b(3,natom2)

 
  pii=acos(-1d0)
  internal(1)=internal2(1)
  internal(2)=0d0
  internal(3)=internal2(2)
  internal(4)=0d0
  internal(5)=internal2(3)
!  internal(6)=acos(internal2(4))
  internal(6)=internal2(4)
ref1_temp=ref1
ref2_temp=ref2




!!! set center of mass of ref1 at origin
  call rm_cmass(ref1_temp,mass(1:natom1),natom1,natom1)

!!! set center of mass of ref2 at origin
  call rm_cmass(ref2_temp,mass(natom1+1:natom1+natom2),natom2,natom2)
!!!!find c.m. for frag2
  cm(1)=0d0
  cm(2)=0d0
  cm(3)=internal(1)

alpha1=0d0
gamma1=internal(2)
!beta1=internal(3)
beta1=acos(internal(3))
gamma2=internal(4)
!beta2=internal(5)
beta2=acos(internal(5))
alpha2=-internal(6)





U_rot(1,1)=cos(alpha1)*cos(beta1)*cos(gamma1)-sin(alpha1)*sin(gamma1)
U_rot(1,2)=-cos(alpha1)*cos(beta1)*sin(gamma1)-sin(alpha1)*cos(gamma1)
U_rot(1,3)=cos(alpha1)*sin(beta1)

U_rot(2,1)=sin(alpha1)*cos(beta1)*cos(gamma1)+cos(alpha1)*sin(gamma1)
U_rot(2,2)=-sin(alpha1)*cos(beta1)*sin(gamma1)+cos(alpha1)*cos(gamma1)
U_rot(2,3)=sin(alpha1)*sin(beta1)

U_rot(3,1)=-sin(beta1)*cos(gamma1)
U_rot(3,2)=sin(beta1)*sin(gamma1)
U_rot(3,3)=cos(beta1)



call vec_to_mat2(ref1_temp,cart_ref1,natom1)
call rotmol(natom1,cart_ref1,cart_ref1b,U_rot)
cart_ref1=cart_ref1b
call mat_to_vec2(cart_ref1,ref1_temp,natom1)


U_rot(1,1)=cos(alpha2)*cos(beta2)*cos(gamma2)-sin(alpha2)*sin(gamma2)
U_rot(1,2)=-cos(alpha2)*cos(beta2)*sin(gamma2)-sin(alpha2)*cos(gamma2)
U_rot(1,3)=cos(alpha2)*sin(beta2)

U_rot(2,1)=sin(alpha2)*cos(beta2)*cos(gamma2)+cos(alpha2)*sin(gamma2)
U_rot(2,2)=-sin(alpha2)*cos(beta2)*sin(gamma2)+cos(alpha2)*cos(gamma2)
U_rot(2,3)=sin(alpha2)*sin(beta2)

U_rot(3,1)=-sin(beta2)*cos(gamma2)
U_rot(3,2)=sin(beta2)*sin(gamma2)
U_rot(3,3)=cos(beta2)


call vec_to_mat2(ref2_temp,cart_ref2,natom2)
call rotmol(natom2,cart_ref2,cart_ref2b,U_rot)
cart_ref2=cart_ref2b
call mat_to_vec2(cart_ref2,ref2_temp,natom2)


!!!!!!!!
do k=1,natom2
   do kp=1,3
      ref2_temp((k-1)*3+kp)=ref2_temp((k-1)*3+kp)+cm(kp)      
   enddo
enddo



!!!!!!!!!!!
do i=1,3*natom2
   cart(3*natom1+i)=ref2_temp(i)
enddo


!!!!!!!!!!!
do i=1,3*natom1
   cart(i)=ref1_temp(i)
enddo

  return
end subroutine INT_Cart




subroutine rm_cmass(cart,mass,natom,natom1)
integer :: k,kp,natom,natom1
real*8 :: mass(natom),cart(natom*3),mtot,cmass1(3)
mtot=0d0
do k=1,natom1
   mtot=mtot+mass(k)
enddo
cmass1=0d0
do k=1,natom1
   do kp=1,3
      cmass1(kp)=cmass1(kp)+cart((k-1)*3+kp)*mass(k)
   enddo
enddo
cmass1=cmass1/mtot

do k=1,natom
   do kp=1,3
      cart((k-1)*3+kp)=cart((k-1)*3+kp)-cmass1(kp)      
   enddo
enddo
return
end subroutine rm_cmass


subroutine vec_to_mat2(cart_perms,cart_mat,natom)
integer :: k,kp,natom
real*8 :: cart_perms(3*natom),cart_mat(3,natom)
do k=1,natom
   do kp=1,3
      cart_mat(kp,k)=cart_perms((k-1)*3+kp)
   enddo
enddo
return
end subroutine vec_to_mat2

subroutine mat_to_vec2(cart_mat,cart_perms,natom)
integer :: k,kp,natom
real*8 :: cart_perms(3*natom),cart_mat(3,natom)
do k=1,natom
   do kp=1,3
      cart_perms((k-1)*3+kp)=cart_mat(kp,k)
   enddo
enddo
return
end subroutine mat_to_vec2

!c-----------------------------------------------------------------------
function func_actual_seed(xi)
USE dynamic_parameters
implicit none
REAL*8, DIMENSION(:), INTENT(IN) :: xi
REAL*8 :: func_actual_seed
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,jp,jj,kk,R,M
integer :: count
real*8 :: temp,weight,norm,somme,jac3(4),jac4(4),temp1,temp2,diff(4),pii
real*8,allocatable :: ind7(:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)
pii=acos(-1d0)
jac3=xi
allocate(ind7(count_seed),ind8(count_seed),PM1(0:order_2_min+1,0:order_2_min+1),PM2(0:order_3_min+1,0:order_3_min+1),PD1(0:order_2_min+1,0:order_2_min+1),PD2(0:order_3_min+1,0:order_3_min+1))
count=0
do ip=1,count_seed 
   count=count+1
   Jac4=coords_seed(ip,:)
   call dist_metric(jac3,jac4,W_a,somme)
   somme=somme**2
!!$!!!!!!!!!!!diff including periodicity
!!$     Jac4=jac3-Jac4
!!$!!!!!!!!!!minimum distance to reflected sym partners
!!$      if(jac4(4)>pii)then
!!$         jac4(4)=jac4(4)-2d0*pii
!!$      endif
!!$      if(jac4(4)<-pii)then
!!$         jac4(4)=jac4(4)+2d0*pii
!!$      endif
!!$   somme=0d0
!!$   jac4(1)=W_a*jac4(1)
!!$   do k=1,4
!!$      somme=somme+jac4(k)**2
!!$   enddo
   ind7(count)=exp(-((somme)/d_seed(ip)**2))/(((somme)/d_seed(ip)**2)**(zz_low/2)+epss)
enddo

call indexxy(count_seed,ind7,ind8)



quitt=0
norm=0d0
do ip=1,count_seed
   quitt=quitt+1
   norm=norm+ind7(ind8(count_seed+1-ip))
   if(norm/ind7(ind8(count_seed+1-ip))>1d8) goto 12
enddo
!! outputs how many expansions are included in interpolation
!write(701,*) quitt 

!!!!
12 jac3=xi
Jac4=jac3
jac4(1)=exp(alpha*jac4(1))
!jac4(4)=acos(jac4(4))
call LPMN(order_2_min+1,order_2_min,order_2_min,jac4(2),PM1,PD1)
call LPMN(order_3_min+1,order_3_min,order_3_min,jac4(3),PM2,PD2)


!norm=0d0
temp=0d0
do i=1,quitt
   jj=ind8(count_seed+1-i)
!   if(pot(jj)<E_limit)then
      weight=ind7(ind8(count_seed+1-i)) 
      !   write(665,*) weight
!      norm=norm+weight
      temp=temp+weight*b2_seed(1,jj)
      count=1
      do R=1,order_1_min
         do L1=0,order_2_min
            do L2=0,order_3_min
               if((L1+L2)<order_4_min+1)then
                  do M=0,min(L1,L2)
                     count=count+1
                     temp=temp+weight*b2_seed(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
                     
                  enddo
               endif
            enddo
         enddo
      enddo
!   endif
   
enddo

func_actual_seed=temp/norm

return
end function func_actual_seed

function dfunc_actual_anal1(xi)
USE dynamic_parameters
implicit none
REAL*8, DIMENSION(:), INTENT(IN) :: xi
REAL*8, DIMENSION(size(xi)) :: x2,x3
REAL*8, DIMENSION(size(xi)+1) :: dfunc_actual_anal1
real*8 :: tampon,tampon2,scale
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,count,jp,jj,kk,R,M
real*8 :: temp(4),weight,norm,h2wn,jac3(4),jac4(4),temp2(4+1),temp3,temp4,diff(4)
real*8,allocatable :: ind7(:),weight_grad(:,:),weight_grad2(:,:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)
real*8 :: somme,somme2,norm_grad(4),grad2(4),valeur,pii
pii=acos(-1d0)
allocate(ind7(count3),ind8(count3))

jac3=xi
count=0
do ip=1,count3 
   count=count+1
   Jac4=coords(ip,:)
   call dist_metric(jac3,jac4,W_a,somme)
   somme=somme**2


!!$!!!!!!!!!!!diff including periodicity
!!$     Jac4=jac3-Jac4
!!$!!!!!!!!!!minimum distance to reflected sym partners
!!$      if(jac4(4)>pii)then
!!$         jac4(4)=jac4(4)-2d0*pii
!!$      endif
!!$      if(jac4(4)<-pii)then
!!$         jac4(4)=jac4(4)+2d0*pii
!!$      endif
!!$   somme=0d0
!!$   jac4(1)=W_a*jac4(1)
!!$   do k=1,4
!!$      somme=somme+jac4(k)**2
!!$   enddo
   ind7(count)=exp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss)
enddo

call indexxy(count3,ind7,ind8)


quitt=0
norm=0d0
do ip=1,count3
   quitt=quitt+1
   norm=norm+ind7(ind8(count3+1-ip))
   if(norm/ind7(ind8(count3+1-ip))>1d8) goto 13
!   if(pot(ind8(count3+1-ip))<E_limit)then
!      norm=norm+ind7(ind8(count3+1-ip))
!   endif
!   quitt=quitt+1
enddo
13 allocate(weight_grad(quitt,4),PM1(0:order_2+1,0:order_2+1),PM2(0:order_3+1,0:order_3+1),PD1(0:order_2+1,0:order_2+1),PD2(0:order_3+1,0:order_3+1))
weight_grad=0d0
norm_grad=0d0
!write(*,*) 'made it15'
do i=1,quitt

   jj=ind8(count3+1-i)
!   if(pot(jj)<E_limit)then
      Jac4=coords(jj,:)
      scale=sqrt((1d0-jac3(2)**2)*(1d0-jac4(2)**2)*(1d0-jac3(3)**2)*(1d0-jac4(3)**2))
      call dist_metric(jac3,jac4,W_a,somme)
      somme=somme**2

      jac4(1)=(jac3(1)-jac4(1))*(W_a**2)
      if(sqrt(1d0-jac3(2)**2)>1d-14)then
         jac4(2)=(acos(jac3(2))-acos(jac4(2)))/sqrt(1d0-jac3(2)**2)
      else
         jac4(2)=(acos(jac3(2))-acos(jac4(2)))/1d-14
      endif
      if(sqrt(1d0-jac3(3)**2)>1d-14)then
         jac4(3)=(acos(jac3(3))-acos(jac4(3)))/sqrt(1d0-jac3(3)**2)
      else
         jac4(3)=(acos(jac3(3))-acos(jac4(3)))/1d-14
      endif
      jac4(4)=jac3(4)-jac4(4)
      if(jac4(4)>pii)then
         jac4(4)=jac4(4)-2d0*pii
      endif
      if(jac4(4)<-pii)then
         jac4(4)=jac4(4)+2d0*pii
      endif
      jac4(4)=jac4(4)*scale
!!!!!!!!!!!diff including periodicity
!!$      Jac4=jac3-Jac4
      !!!!!!!!!!minimum distance to reflected sym partners
!!$      if(jac4(4)>pii)then
!!$         jac4(4)=jac4(4)-2d0*pii
!!$      endif
!!$      if(jac4(4)<-pii)then
!!$         jac4(4)=jac4(4)+2d0*pii
!!$      endif
!!$      somme=0d0

!!$      jac4(1)=W_a*jac4(1)
!!$      do k=1,4
!!$         somme=somme+jac4(k)**2
!!$      enddo

      somme2=0d0
      somme2=somme/(d(jj)**2)
      temp=0d0
      if(somme>1d-10)then
         do ip=1,4
            temp(ip)=Jac4(ip)*(-2d0)*ind7(ind8(count3+1-i))*&
                 ((1.0d0/(d(jj)**2))+(zz/2)*((somme2**(zz/2))/((somme2**(zz/2))+epss))*(1.0d0/(somme)))
         enddo
      else
         temp=0d0
      endif
      weight_grad(i,:)=temp
      do ip=1,4
         norm_grad(ip)=norm_grad(ip)+weight_grad(i,ip)
      enddo
!   endif
enddo
!write(*,*) 'made it16'
!!!!
temp2=0d0
Jac4=jac3
jac4(1)=exp(alpha*jac4(1))
! jac4(4)=acos(jac4(4))
call LPMN(order_2+1,order_2,order_2,jac4(2),PM1,PD1)
call LPMN(order_3+1,order_3,order_3,jac4(3),PM2,PD2)
do i=1,quitt
   jj=ind8(count3+1-i)
!   if(pot(jj)<E_limit)then
      
      temp=0d0
      grad2=0d0
      valeur=0d0
      weight=ind7(ind8(count3+1-i)) 
      
      
      valeur=valeur+weight*b2(1,jj)
      
      count=1
      do R=1,order_1
         do L1=0,order_2
            do L2=0,order_3
               if((L1+L2)<order_4+1)then
                  do M=0,min(L1,L2)
                     count=count+1
                     valeur=valeur+weight*b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
                     
                     grad2(1)=grad2(1)+b2(count,jj)*dble(R)*alpha*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
                     grad2(2)=grad2(2)+b2(count,jj)*(jac4(1))**(R)*PD1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
                     grad2(3)=grad2(3)+b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PD2(M,L2)*cos(dble(M)*jac4(4))
                     grad2(4)=grad2(4)+b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*(-dble(M)*sin(dble(M)*jac4(4)))
                     
                  enddo
               endif
            enddo
         enddo
      enddo
      
      
      temp2(1)=temp2(1)+valeur
      do k=1,4
         temp(k)=(weight/norm)*grad2(k)
      enddo
      temp2(2:4+1)=temp2(2:4+1)+temp
      do k=1,4
         temp2(k+1)=temp2(k+1)+&
              (valeur/(weight*norm))*weight_grad(i,k)-&
              (1.0d0/norm**2)*valeur*norm_grad(k)
      enddo
!   endif
enddo
   



dfunc_actual_anal1(1)=temp2(1)/norm
dfunc_actual_anal1(2:5)=temp2(2:5)

deallocate(ind7,ind8,weight_grad)




return
end function dfunc_actual_anal1


function func_actual(xi)
USE dynamic_parameters
implicit none
REAL*8, DIMENSION(:), INTENT(IN) :: xi
REAL*8 :: func_actual
integer :: i,j,k,ipp,jpp,ip,quitt,l1,l2,l3,l4,count2,l,jp,jj,kk,R,M
integer :: count
real*8 :: temp,weight,norm,somme,jac3(4),jac4(4),temp1,temp2,diff(4),pii
real*8,allocatable :: ind7(:),PM1(:,:),PM2(:,:),PD1(:,:),PD2(:,:)
integer,allocatable :: ind8(:)
pii=acos(-1d0)
jac3=xi
allocate(ind7(count3),ind8(count3),PM1(0:order_2+1,0:order_2+1),PM2(0:order_3+1,0:order_3+1),PD1(0:order_2+1,0:order_2+1),PD2(0:order_3+1,0:order_3+1))
count=0
do ip=1,count3 
   count=count+1
   Jac4=coords(ip,:)
   call dist_metric(jac3,jac4,W_a,somme)
   somme=somme**2

!!$!!!!!!!!!!!diff including periodicity
!!$     Jac4=jac3-Jac4
!!$!!!!!!!!!!minimum distance to reflected sym partners
!!$      if(jac4(4)>pii)then
!!$         jac4(4)=jac4(4)-2d0*pii
!!$      endif
!!$      if(jac4(4)<-pii)then
!!$         jac4(4)=jac4(4)+2d0*pii
!!$      endif
!!$   somme=0d0
!!$   jac4(1)=W_a*jac4(1)
!!$   do k=1,4
!!$      somme=somme+jac4(k)**2
!!$   enddo
   ind7(count)=exp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss)
enddo

call indexxy(count3,ind7,ind8)



quitt=0
norm=0d0
do ip=1,count3
   quitt=quitt+1
   norm=norm+ind7(ind8(count3+1-ip))
   if(norm/ind7(ind8(count3+1-ip))>1d8) goto 12
enddo
!! outputs how many expansions are included in interpolation
!write(701,*) quitt 

!!!!
12 jac3=xi
Jac4=jac3
jac4(1)=exp(alpha*jac4(1))
!jac4(4)=acos(jac4(4))
call LPMN(order_2+1,order_2,order_2,jac4(2),PM1,PD1)
call LPMN(order_3+1,order_3,order_3,jac4(3),PM2,PD2)


!norm=0d0
temp=0d0
do i=1,quitt
   jj=ind8(count3+1-i)
!   if(pot(jj)<E_limit)then
      weight=ind7(ind8(count3+1-i)) 
      !   write(665,*) weight
!      norm=norm+weight
      temp=temp+weight*b2(1,jj)
      count=1
      do R=1,order_1
         do L1=0,order_2
            do L2=0,order_3
               if((L1+L2)<order_4+1)then
                  do M=0,min(L1,L2)
                     count=count+1
                     temp=temp+weight*b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*PM2(M,L2)*cos(dble(M)*jac4(4))
                     
                  enddo
               endif
            enddo
         enddo
      enddo
!   endif
   
enddo

func_actual=temp/norm

return
end function func_actual

subroutine dist_metric(jac,jac2,scale,dist)
integer :: i,j
real*8 :: jac(4),jac2(4),scale,dist,temp(4),pii
pii=acos(-1d0)
temp(1)=((jac(1)-jac2(1))*scale)**2
temp(2)=(acos(jac(2))-acos(jac2(2)))**2
temp(3)=(acos(jac(3))-acos(jac2(3)))**2
temp(4)=jac(4)-jac2(4)
if(temp(4)>pii)then
   temp(4)=temp(4)-2d0*pii
endif
if(temp(4)<-pii)then
   temp(4)=temp(4)+2d0*pii
endif
temp(4)=(temp(4)**2)*sqrt((1d0-jac(2)**2)*(1d0-jac2(2)**2)*(1d0-jac(3)**2)*(1d0-jac2(3)**2))
dist=0d0
do i=1,4
   dist=dist+temp(i)
enddo
dist=sqrt(dist)
return
end subroutine dist_metric
