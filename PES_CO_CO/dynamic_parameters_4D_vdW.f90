MODULE dynamic_parameters
implicit none
save
public
real*8, allocatable :: b2(:,:),b2_lower(:,:),b2_minimal(:,:),b2_seed(:,:),d(:),d_seed(:),Jac(:),Jac2(:),coords(:,:),coords_seed(:,:)
real*8, allocatable :: cart(:),bdist(:)
real*8, allocatable :: pot(:),pot_seed(:),grad(:,:),grad_seed(:,:),mass(:),dcart(:),ref1(:),ref2(:)
integer :: basis_1,basis_2,basis_3,order_1,order_2,order_3,order_4,count3,zz,zz_low,zz4,support,count7,ddd,myid,natom,ab_flag,ab_flag2,lab,permfac,maxpoints,natom1,natom2,nbdist,order_1_min,order_2_min,order_3_min,order_4_min,count_seed,low_grid,subzero,dist_flag
real*8 :: epss,acc,rmax(4),rmin(4),poten,hartokcl,ugrad,Max_E,Max_E_seed,E_range,dist_tol,E_limit,W_a,alpha,Glob_min,ass
character(len=3),allocatable :: symb(:)
END MODULE dynamic_parameters
