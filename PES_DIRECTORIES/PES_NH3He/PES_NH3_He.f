      SUBROUTINE USER_DEFINED_PES(V,R,rvib,rvib2,alpha,beta,gamma,
     &										   aalpha,bbeta,ggamma)
!     INPUT:  R - distance betweenn COMs, rvib - vibrational coordinate
! 		      alpha,beta,gamma - Euler's angles of the first molecule
!   		  aalpha,bbeta, ggamma - Euler's angles of the second molecule
!     OUTPUT: V - value of the potential
      IMPLICIT NONE
      REAL*8 V,R,alpha,beta,gamma,aalpha,bbeta,ggamma,t,rvib,rvib2
	  real*8, parameter ::  pi = 4.0d0*datan(1.0d0), rho = 112.150d0
	  call potential(R,beta*180d0/pi,gamma*180d0/pi,rho,V)
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


      subroutine potential( Rcm, theta, phi, rho, pot )

C     Calculates NH_3-He potential in monomer frame coordinates.
C     Consider the monomer frame and the vector R as explained in the paper. 
C     Theta is the zenith angle of this vector, phi the azimuth angle and Rcm
C     the length of this vector. Rho is as explained in the paper.
C
C     Note that the angles theta and phi are slightly different from the angles 
C     beta^bf and gamma^bf that are also used in the paper to describe the 
C     orientation of the monomer frame in the body-fixed frame.
C
C     The two sets of coordinates are easily related and can be transformed 
C     into each other in the following way: beta^bf = theta and 
C     gamma^bf = 180 - phi, where the angles are given in degrees.
C
C     In the routine the expansion coefficients v_lmp and c_lmpn are 
C     interpolated in R by means of radial RKHS fit. The data to reconstruct 
C     this fit is stored in a data block.
C
C     Input:
C     Rcm has to be given in a0. Theta, phi and rho have to be given in degrees.
C
C     Output:
C     The potential is given in cm^{-1}.
C
C     Range of validity:
C     For the long range fit, the correct analytic terms were used in the fit,
C     so that our potential is expected to be of good quality for all large R. 
C     However, in the short range fit, an RKHS interpolation was used, 
C     whereas the true behavior is exponential. RKHS is useful for interpolation,
C     so that our fit is accurate for for R > 4.0 a_0. If smaller values for
C     R are given, a warning message is printed and the routine stops.
C
C     Test input:
C     Rcm = 6.09d0, theta = 89.1d0, phi = 61.4d0, rho = 112.15d0
C     Test output:
C     potential = -35.04088154834608

      implicit  real*8 (a-h,o-z)

C     Initial values:

      parameter ( n_R_i   = 19,   ! number of R points of ab initio grid
     .            n_lm_i  = 26,   ! number of (l,m) pairs used in expansion
     .            n_k_i   = 5,    ! number of odd or even powers in rho-expansion
     .            n_lr_i  = 9,    ! number of analytic v_lm coefs in the long range
     .            n_p_i   = 2)    ! max nr of powers R^-n for each analytic coef v_lmp

      real*8
     .     Qvec  (n_R_i),
     .     CC    (n_R_i, n_k_i, n_lm_i),
     .     CP    (n_p_i, n_k_i, n_lr_i),
     .     VLM   (n_lm_i), cm1

      integer mRKHS(n_lm_i)
      common/maxqn/n_lm, n_k, n_R, n_lr, n_p

      n_R   = n_R_i
      n_lm  = n_lm_i
      n_k   = n_k_i
      n_lr  = n_lr_i
      n_p   = n_p_i

      pi        = dacos(-1d0)
      theta_rad = theta*pi/180.0d0
      phi_rad   = phi*pi/180.0d0
      
      if (Rcm .lt. 4.0d0) then
         write(*,*) 'Warning: fit not to be trusted for R < 4 a_0'
      endif

C     Load data for expansion coefficients
      call load_VLM(CC, CP, mRKHS)

C     Compute expansion coefficients and potential
      call comp_VLM(Rcm, rho, theta_rad, phi_rad,
     &              CC, CP, pot, VLM, mRKHS)

      cm1       = 4.556335252750422d0 ! cm^{-1} in microHartree
      pot       = pot/cm1 

      end

      SUBROUTINE load_VLM(CC, CP, mRKHS)

C     Loads RKHS parameters for expansion coefficients of
C     potential obtained from R fits

      implicit  real*8 (a-h,o-z)

      integer q_tab(130)
      common/maxqn/n_lm, n_k, n_R, n_lr, n_p
      common/coefs1/clmk_rkhs(2470)
      common/coefs2/clmkp_anlr(90)
      common/qnrs/q_tab

      integer mRKHS(n_lm) 
      real*8 CC(n_R, n_k, n_lm), CP(n_p, n_k, n_lm)

      nRKHS = 2
      num  = 0
      do ilm = 1, n_lm
        mRKHS(ilm) = q_tab(5*ilm);
        do ik = 1, n_k
          do iR = 1, n_R
            num = num + 1
            CC(iR,ik,ilm) = clmk_rkhs(num)
          enddo
        enddo
      enddo

      num  = 0
      do ilm = 1, n_lr
        do ik = 1, n_k
          num = num + 1
          CP(1,ik,ilm) = clmkp_anlr(num)
          num = num + 1
          CP(2,ik,ilm) = clmkp_anlr(num)
        enddo
      enddo

      end

      SUBROUTINE comp_VLM(R0, rho0, theta0, phi0,
     &                    CC, CP, pot, VLM, mRKHS)
C
C     Returns L, M expansion coefficients of potential at R0 and rho0
C
      implicit  real*8 (a-h,o-z)

      integer mRKHS(n_lm), q_tab(130), p_lr, x
      parameter (lmax=10)

      common/maxqn/n_lm, n_k, n_R, n_lr, n_p
      common/qnrs/q_tab
      common/grid/Rs(19)

      real*8 CC(n_R, n_k, n_lm), CP(n_p, n_k, n_lm), Qvec(n_R)
      real*8 VLM(n_lm), VLMK(n_k,n_lm), tes((lmax+1)*(lmax+1))

      nRKHS = 2
      pot   = 0.0d0
      call tesser(tes,lmax,theta0,phi0)

      do ilm = 1, n_lr
        VLM(ilm) = 0.d0
        p_lr     = q_tab(5*ilm-1)
        li       = q_tab(5*ilm-3)
        mi       = q_tab(5*ilm-2)
        x        = li + mi
        Call  COMP_QVEC(Qvec, Rs, n_R, R0, nRKHS, mRKHS(ilm))
        do ik = 1, n_k
          VLMK(ik,ilm) = 0.d0
          do iR = 1, n_R
            VLMK(ik,ilm)=VLMK(ik,ilm)+CC(iR,ik,ilm)*Qvec(iR)
          enddo
          call tt_df(y_df1,R0,p_lr)
          call tt_df(y_df2,R0,p_lr+2)
          alr_part = CP(1,ik,ilm)*y_df1*(R0)**(-p_lr)
     &             + CP(2,ik,ilm)*y_df2*(R0)**(-p_lr-2)
          VLM(ilm) = VLM(ilm) + (VLMK(ik,ilm)+alr_part)
     &             *(rho0 - 90.0d0)**(2*(ik-1)+(x-2*int(x/2)))
        enddo

C       The tesseral t(l,m,theta,phi) is accessed via tes(l*(l+1)+m+1)
        ang_fun = tes(li*(li+1)+mi+1)
        pot     = pot + VLM(ilm)*ang_fun
      enddo

      do ilm = n_lr + 1, n_lm
        VLM(ilm) = 0.d0
        li       = q_tab(5*ilm-3)
        mi       = q_tab(5*ilm-2)
        x        = li + mi
        Call  COMP_QVEC(Qvec, Rs, n_R, R0, nRKHS, mRKHS(ilm))
        do ik = 1, n_k
          VLMK(ik,ilm) = 0.d0
          do iR = 1, n_R
            VLMK(ik,ilm) = VLMK(ik,ilm) + CC(iR,ik,ilm)*Qvec(iR)
          enddo
          VLM(ilm) = VLM(ilm) + VLMK(ik,ilm)
     &             *(rho0 - 90.0d0)**(2*(ik-1)+(x-2*int(x/2)))
        enddo

C       The tesseral t(l,m,theta,phi) is accessed via tes(l*(l+1)+m+1)
        ang_fun = tes(li*(li+1)+mi+1)
        pot     = pot + VLM(ilm)*ang_fun
      enddo

      end

      SUBROUTINE COMP_QVEC(Qvec, R, n_R,  X, nRKHS, mRKHS)

!     Evaluate RKHS function Q(R, X)
!     [Eq. (17) in J.C.P., Vol. 104, p. 2584, (1995)]

      implicit real*8(a-h,o-z)
      parameter (lp = 10)
      real*8 Qvec(n_R),
     .       R(n_R),
     .       P(0:lP)
      save

      if (lP .lt. nRKHS-1) then
         write(*,*) 'lP too small in subroutine COMP_QVEC'
         stop 16
      endif

C     Return polynomial coefficients in P
      call setup(P, nRKHS, mRKHS)

      do i = 1,n_R
           Qvec(i) = eval(X, R(i), P, nRKHS, mRKHS)
      enddo

      end

      REAL*8 FUNCTION EVAL(X, R, P, n, m)

C     Evaluate polynomial in x< / x>

      implicit real*8(a-h,o-z)

      real*8 P(0:n)

      xlarge = max(X, R)
      xsmall = min(X, R)
      xsl    = xsmall/xlarge

      xx  = xsl
      Pol = P(0)
      do i = 1, n-1
          Pol = Pol + P(i)*xx
          xx  = xx*xsl
      enddo

      eval = Pol/xlarge**(m+1)

      end

      SUBROUTINE SETUP(P, n, m)
      implicit real*8(a-h,o-z)
      parameter (kmax = 50)
      common /bin/  binom(0:kmax, 0:kmax), jmax
      real*8 P(0:n)

C     Fill out BINOM by Pascal's triangle relation (if necessary)
      IF (jmax .EQ. 0) CALL Pascal

      if (n+m .ge. jmax) then
          write(*,'(a, i5)')
     .        'Increase kmax in subroutine Pascal to', n+m
          stop 16
      endif

!     Get expansion coefficients of hypergeometric function {_2}F_1
      Call gauss_2f1(-n+1, m+1, n+m+1, P, n)

      obin = dble(n)/binom(n+m, n)
      do i = 0,n
        P(i) = P(i)*obin
      enddo

      end

      SUBROUTINE GAUSS_2f1(i, j, k, P, n)
!     Hypergeometric function {_2}F_1,
!     see Abramowitz & Stegun 15.1.1,  15.4.1.

C     {_2}F_1 is a polynomial of degree n,
C     if either i or j is 0, -1, -2, ...

C     GAUSS_2F1 returns the expansion coefficients of
C     this polynomial in P(0..n).

      implicit real*8(a-h,o-z)
      real*8 P(0:n)
      parameter (ihuge = 2**30)

      nmax = ihuge
      if (i .le. 0) then
        nmax = -i
      endif
      if (j .le. 0) then
         nmax = min(nmax,-j)
      endif

      if (nmax .eq. ihuge) then
        write(*,*) ' Either i or j must be 0,-1,-2,...'
      endif

      lfac=1
      do  l = 0, nmax
         Pi   = pochhammer(dble(i), l)
         Pj   = pochhammer(dble(j), l)
         Pk   = pochhammer(dble(k), l)
         p(l) = Pi*Pj/(dble(lfac)*Pk)
         lfac = lfac*(l+1)
      enddo

      end

      REAL*8 FUNCTION POCHHAMMER(z, n)
C     Pochhammer's Symbol, see Abramowitz & Stegun 6.1.22
C     POCHHAMMER(Z,N) returns P=(Z)_N.
C
C     Definition:  (Z)_0 = 1
C                  (Z)_N = Z*(Z+1)*...*(Z+N-1)

      implicit real*8(a-h,o-z)
      real*8 P(0:n)

      P(0) = 1.d0
      do i = 1, n
         P(i) = P(i-1)*(z + dble(i-1))
      enddo
      pochhammer = P(n)

      end

      SUBROUTINE PASCAL

C     Fill common bin with binomial coefficients (Pascal's triangle)

      implicit double precision(a-h,o-z)
      parameter (kmax = 50)
      common /bin/  binom(0:kmax, 0:kmax), jmax

      if ( kmax .eq. jmax ) return  ! bin has been set

C     Set jmax as sign that binom has been set:
      jmax = kmax

      binom(0,0) = 1.d0
      do  i=1,jmax
         binom(i,0) = 1.d0
         binom(i,i) = 1.d0
         do  j=1,i-1
            binom(i,j) = binom(i-1,j-1) + binom(i-1,j)
         enddo
      enddo

      end

      SUBROUTINE tt_df(y, x, p_lr)
!     Tang-Toennies damping function

      IMPLICIT NONE

      double precision x, y, a_tt, sum_val
      integer i, p_lr, fac

      a_tt    = 2.0882198434643375d0  ! Calculated for NH3-He
      sum_val = 0
      do i = 0,p_lr
        call fact(fac,i)
        sum_val = sum_val + (a_tt*x)**i/fac
      enddo

      y = 1.0d0-sum_val*exp(-a_tt*x)
      end

      SUBROUTINE fact(fac,n)
C     Factorial

      IMPLICIT NONE

      integer fac, n, i

      fac = 1
      do i = 2, n
	fac = fac * i
      enddo

      end

      SUBROUTINE TESSER(tes,n,theta,phi)
C
C     Subroutine computes tesseral harmonics as defined by J.L. Prather
C     [N.B.S. monograph 19 (1961)].
C
C     The tesserals are computed from l=m=0 upwards to l=m=n for the
C     polar angles theta and phi (radians).
C     The results are stored in the linear array tes in increasing order
C     of l and m.
C     Negative m-values refer to prather's s(l,m) (= p(l,m)*sin(m*phi))
C     positive m-values to prather's c(l,m) ( = p(l,m) * cos(m*phi) )
C     Note:
C         the tesseral t(l,m,theta,phi) is accessed via tes(l*(l+1)+m+1)
C
      implicit double precision (a-h,o-z)
      parameter ( nn = 50 )
      dimension tes( (n+1)*(n+1) )
      dimension p( (nn+1)*(nn+2)/2 ),sn(nn+1),cs(nn+1)
      data tpih/3.989422804014327d-1/
C     (   tpih = 1/dsqrt(2*pi)   )
      data twmh/7.071067811865475d-1/
C     (   twmh = 1./dsqrt(2)    )
C
      if (n .gt. nn) then
         write(*,'(''0n larger than:'',i3,'' increase nn in'',
     1              '' subroutine tesser '' )' ) nn
         stop 16
      endif
C
C     Compute associated Legendre functions
C
      cost = dcos(theta)
      Call assleg(p,n,cost)
C
      tes(1) = tpih * twmh * p(1)
      if ( n .eq. 0) return
C
C     Compute necessary sines and cosines
C
      cs(1) = 1.d0
      sn(1) = 0.d0
      cosp = dcos(phi)
      sinp = dsin(phi)
      do 10 m=1,n
      cs(m+1) = cs(m)*cosp - sn(m)*sinp
      sn(m+1) = sn(m)*cosp + cs(m)*sinp
   10 continue
C
      kp = 1
      ll = 1
      dl = 0.d0
C
      do 40 l=1,n
         ll = ll + l + l
C        ( ll = l(l+1) + 1   )
         kp = kp + 1
         dl = dl + 1.d0
         fac = dsqrt(dl+dl+1.d0) * tpih
         tes(ll) = fac*twmh * p(kp)
         kt1 = ll
         kt2 = ll
         dlm = dl
         dlm1 = dl + 1.d0
C
         do 30 m=1,l
            kt1 = kt1 + 1
            kt2 = kt2 - 1
C           ( kt1 = l(l+1) +1 + m  )
C           ( kt2 = l(l+1) +1 - m )
            kp = kp + 1
            dlm  = dlm  + 1.d0
            dlm1 = dlm1 - 1.d0
C           ( dlm = l+m      )
C           ( dlm1= l+1-m    )
            fac = fac / dsqrt(dlm*dlm1)
            tes(kt1) = fac * p(kp) * cs(m+1)
            tes(kt2) = fac * p(kp) * sn(m+1)
C           ( t(l,m) = fac * p(l,m) * cos(m*phi)  )
C           ( t(l,-m)= fac * p(l,m) * sin(m*phi)  )
C
   30    continue
   40 continue
C
      end

      SUBROUTINE ASSLEG(p,n,x)
C
C     Subroutine computes associated Legendre polynomials as defined
C     by A.R. Edmonds (angular momentum in quantum mechanics).
C     x is the usual coordinate (-1 < x < +1 ), n is the maximum
C     quantum number. the associated legendre functions are computed
C     in the order (0,0),(1,0),(1,1),(2,0),(2,1),(2,2),.... ,(n,n)
C     and returned in the array p.
C     The associated Legendre function p(l,m,x) may be  accessed via
C     p( l*(l+1)/2 + m+1 )
C
      implicit double precision (a-h,o-z)
      dimension p((n+1)*(n+2)/2 )
C
      p(1) = 1.d0
      if (n .eq. 0) return
      sint = sqrt(1.d0 - x*x)
      p(2) = x
      p(3) = sint
      if (n .eq. 1) return
C
      lm1 = 1
      lm  = 3
      dl  = 1.d0
C
      do 20 l=2,n
      dl = dl + 1.d0
      lm1 = lm1 + 1
C     (   lm1 = l*(l-1)/2 + 1 )
      lm  = lm + 1
C     (   lm = l*(l+1)/2 + 1 )
C
      p(lm) = x*p(lm1) - sint*p(lm1+1)/dl
C     ( p(l,0) = x*p(l-1,0) - dsqrt(1-x*x)*p(l-1,1)/l   )
      mmax = l-1
      dlm = dl
      do 10 m=1,mmax
      lm1 = lm1 + 1
      lm  = lm + 1

      p(lm) = dlm*sint*p(lm1-1) + x*p(lm1)
C     (  p(l,m) = (l+m-1)*dsqrt(1-x*x)*p(l-1,m-1) + x*p(l-1,m)   )
      dlm = dlm + 1.d0
   10 continue
      lm = lm + 1
      p(lm) = (dl+dl-1.d0)*sint*p(lm1)
C     (  p(l,l) = (2*l-1)*dsqrt(1-x*x)*p(l-1,l-1)    )

   20 continue

      END

      BLOCK DATA BLPAS
      implicit double precision(a-h,o-z)
      PARAMETER (kmax = 50)
      common /bin/  binom(0:kmax, 0:kmax), jmax
      data jmax/0/
      END

      BLOCK DATA DATA_POT

      implicit real*8(a-h,o-z)

      integer q_tab(130)

      common/maxqn/n_lm, n_k, n_R, n_lr, n_p
      common/coefs1/clmk_rkhs(2470)
      common/coefs2/clmkp_anlr(90)
      common/qnrs/q_tab
      common/grid/Rs(19)
      
      data Rs/
     .   4.00d0,
     .   4.50d0,
     .   5.00d0,
     .   5.50d0,
     .   6.00d0,
     .   6.50d0,
     .   7.00d0,
     .   7.50d0,
     .   8.00d0,
     .   8.50d0,
     .   9.00d0,
     .   9.50d0,
     .  10.00d0,
     .  12.00d0,
     .  14.40d0,
     .  17.30d0,
     .  20.80d0,
     .  25.00d0,
     .  30.00d0/

C     Quantum numbers and related parameters (ilm,l,m,p,m_rkhs)

      data q_tab/
     .  1,  0,  0,  6,  7,
     .  2,  1,  0,  7,  8,
     .  3,  2,  0,  6,  7,
     .  4,  3,  0,  7,  8,
     .  5,  3,  3,  7,  8,
     .  6,  4,  0,  8,  9,
     .  7,  4,  3,  8,  9,
     .  8,  5,  0,  9, 10,
     .  9,  5,  3,  9, 10,
     . 10,  6,  0, 10,  9,
     . 11,  6,  3, 10,  9,
     . 12,  6,  6, 10,  9,
     . 13,  7,  0, 11, 10,
     . 14,  7,  3, 11, 10,
     . 15,  7,  6, 11, 10,
     . 16,  8,  0, 12, 11,
     . 17,  8,  3, 12, 11,
     . 18,  8,  6, 12, 11,
     . 19,  9,  0, 13, 12,
     . 20,  9,  3, 13, 12,
     . 21,  9,  6, 13, 12,
     . 22,  9,  9, 13, 12,
     . 23, 10,  0, 14, 13,
     . 24, 10,  3, 14, 13,
     . 25, 10,  6, 14, 13,
     . 26, 10,  9, 14, 13/
     

C     The expansion coefficients for the analytic longe-range fits

      data clmkp_anlr/
     +      -3.8700307675481431d+07,      -9.0205221024120700d+08,
     +      -2.3336960953173725d+03,       1.6067500354906858d+05,
     +       1.6673166129332706d+00,      -7.8204893416583774d+00,
     +      -1.0857536761105673d-03,      -5.7158546290696687d-02,
     +       2.7978100888469890d-07,       2.1900630853212775d-05,
     +      -3.0818483973665035d+05,       2.2291595021733887d+07,
     +       1.2184400504415405d+03,       7.3931217218918828d+03,
     +      -7.4904144496737035d-01,      -2.1583061230092497d+01,
     +       3.8096527466545242d-04,       1.7816157410011729d-02,
     +      -9.0379582829250257d-08,      -4.8557515861512874d-06,
     +       5.0307863246417306d+04,      -6.1215063087038152d+07,
     +      -1.8171208173245673d+03,       8.9344037547261760d+04,
     +       6.0046877065415238d-01,      -8.8863010022982309d+01,
     +      -2.2102750989393390d-04,       2.8587131559813100d-02,
     +       4.6601969612384290d-08,      -3.9646272291836229d-06,
     +      -6.6283490545712772d+05,      -2.5904288166493963d+07,
     +       5.3785559325914733d+02,       2.6202171716298210d+04,
     +      -1.7468574425380329d-01,      -3.9778418177447534d+01,
     +       2.1480836911634478d-05,       3.5141671486117208d-02,
     +       4.3770836054218590d-09,      -1.0612003506926785d-05,
     +      -8.8420825149549916d+06,      -1.0643560385043191d+09,
     +       3.1273455717185352d+03,       4.7540809534960700d+05,
     +      -6.5998116455049563d-02,      -3.5364365419348326d+02,
     +      -4.1486859338056754d-04,       3.5924396286940291d-01,
     +       1.6143889416925463d-07,      -1.2549664064666309d-04,
     +      -2.6181046813579984d+07,             0000000000000000,
     +       7.3188182426663669d+04,             0000000000000000,
     +      -3.8168850243799831d+01,             0000000000000000,
     +       8.1594690760605684d-03,             0000000000000000,
     +      -8.5178464709396469d-07,             0000000000000000,
     +       1.8643927751594579d+06,             0000000000000000,
     +      -8.5606737946434521d+02,             0000000000000000,
     +      -8.4519352791411761d-01,             0000000000000000,
     +       1.2776722979959971d-03,             0000000000000000,
     +      -4.4026220538566993d-07,             0000000000000000,
     +       4.6497807305735676d+06,             0000000000000000,
     +       3.0637295979714377d+03,             0000000000000000,
     +      -2.5769969889066964d+01,             0000000000000000,
     +       2.7122503215929877d-02,             0000000000000000,
     +      -8.4632350028629710d-06,             0000000000000000,
     +       7.9537369000529379d+07,             0000000000000000,
     +      -1.9185917958427058d+05,             0000000000000000,
     +       1.1429714547767929d+02,             0000000000000000,
     +      -8.3616973255182400d-02,             0000000000000000,
     +       3.0530520818631014d-05,             0000000000000000/

C     The data for the RHKS fit of VLMK(R)

      data clmk_rkhs/
     +       1.0239899698113600d+12,      -1.1829791682613491d+12,
     +       4.2601581096608832d+11,      -3.4300102442881299d+11,
     +       1.3115120908917360d+10,      -9.6006489524694290d+10,
     +      -7.6149488240502853d+09,      -3.7600268522765791d+08,
     +       1.7488945507543800d+10,       2.2304810943702549d+10,
     +       1.9209715115455891d+10,       6.0259933905202599d+10,
     +      -2.0508913711621460d+10,       1.1532943973942830d+11,
     +      -6.3060660380513496d+10,       2.7922636536362011d+10,
     +      -1.5694696078087851d+10,       6.7911966204607506d+09,
     +      -4.8359322546884184d+09,      -2.7125224398384109d+07,
     +       1.8133431986392379d+07,      -1.6239571240578840d+07,
     +      -3.2901889789753729d+06,      -8.7568969134970214d+06,
     +       6.2225575850394054d+05,       1.0727774412597191d+07,
     +       1.2178455043036960d+07,       1.1217440174893470d+07,
     +       1.6521370899912180d+07,       3.2778113191998340d+07,
     +      -1.8197325294011850d+07,       2.3866931240764670d+07,
     +      -7.7696040820684418d+07,       3.4543306440278932d+07,
     +      -1.4476159104029950d+07,       7.7020063215309847d+06,
     +      -5.0475617555306526d+06,       3.9671737196068168d+06,
     +      -7.2809649735071371d+03,       1.4859020110317329d+04,
     +      -2.8809084713070142d+03,       7.6452304110213536d+03,
     +       2.5785956292494111d+03,       7.2515448081113627d+03,
     +      -1.0662650338916790d+03,      -3.9504037810710161d+03,
     +       8.7928116819692968d+03,      -2.4084590499000809d+04,
     +      -4.9456388659189754d+04,      -6.1007872338159701d+04,
     +       1.1519002953134950d+05,      -4.5633335459857653d+03,
     +      -3.8024635045821669d+03,       2.7681727851511992d+03,
     +      -1.7189854668689779d+03,       4.8261120791659523d+03,
     +      -7.4593663814581332d+03,       6.2956580796102726d+00,
     +      -8.2672431698424109d+00,       3.8197293880020791d+00,
     +      -2.2259649822293710d+00,       1.5344529839539460d+00,
     +      -6.2255033560020241d+00,      -4.2756074073573576d+00,
     +      -3.7676177859589072d-01,      -1.6571797076438479d+01,
     +       1.9471740262220461d+01,       4.0068082102424952d+01,
     +       1.0479159368136830d+02,      -1.8282960214165041d+02,
     +       6.1033587401633810d+01,      -2.0144924249072691d+01,
     +       5.8958148890938231d+00,      -2.7215061151894608d+00,
     +      -2.9839981077032820d+00,       7.1931025806248670d+00,
     +      -1.8079240272393590d-03,       1.9644364967071680d-03,
     +      -1.0979764936630361d-03,       2.4116890214954541d-04,
     +      -8.1565735882104060d-04,       1.7431415568275319d-03,
     +       1.7560538044663510d-03,       3.9216428820406372d-04,
     +       5.9573789257268434d-03,      -5.7169731985703626d-03,
     +      -1.1766849341935620d-02,      -3.9594603320659547d-02,
     +       6.6736613194902109d-02,      -2.4871222307488269d-02,
     +       8.6141546309069731d-03,      -2.6101388859362819d-03,
     +       1.2332158975012330d-03,       8.3267750704272712d-04,
     +      -2.3442235262119359d-03,      -2.5311060251045101d+10,
     +       2.4907581298130360d+10,      -9.6959752710773697d+09,
     +       9.1291455593774605d+09,       1.2499828334312921d+09,
     +       4.3808537826683578d+09,       5.9984903222833109d+08,
     +      -1.3278742691854720d+08,      -3.9444616256144363d+08,
     +      -4.1750531635703498d+08,       1.6608637402449651d+09,
     +      -3.5163870791833081d+09,       4.5485944120616407d+09,
     +      -9.5590768745733070d+09,       5.1327615181201053d+08,
     +       2.5225415576474290d+09,      -1.0837483995266190d+09,
     +       6.5897751997272003d+08,       1.3535524617923641d+08,
     +       3.0716977453607069d+06,      -3.5350261626095660d+06,
     +       9.9114627657662251d+05,      -4.8833494697067549d+05,
     +       1.0645834956272449d+06,      -1.0002078163919610d+06,
     +       1.8569198799510440d+06,       6.8894034462859416d+06,
     +       2.3547127634531111d+06,      -5.9761186819056924d+06,
     +      -1.7153670705009270d+07,      -2.3777362076269750d+07,
     +       8.0284608190525090d+04,       8.1119744708796561d+07,
     +      -6.2649542969604723d+07,       2.6538563795685899d+07,
     +      -1.5063999940897491d+07,       1.0381691089393150d+07,
     +      -7.3970205824390249d+06,       1.3719425432652531d+03,
     +      -1.0220089288094490d+03,       9.8468631234545899d+02,
     +      -8.2117769622729963d+02,      -1.6884294098184930d+03,
     +       2.8226595256575929d+03,      -5.4666435166189312d+03,
     +      -1.7185030802609770d+04,      -5.5789999586664762d+03,
     +       1.4574403439040620d+04,       3.2833384595446609d+04,
     +       9.9542293554934091d+04,      -5.7237179375963467d+04,
     +      -1.6319783150137699d+05,       1.4131804712925680d+05,
     +      -6.3920406987775619d+04,       3.6790648789616716d+04,
     +      -2.5528976663877849d+04,       1.7783090981899419d+04,
     +      -1.4040104007351790d+00,       1.1076828764076190d+00,
     +      -9.9796453386448791d-01,       3.1007348967124743d-01,
     +       1.4165604950919060d+00,      -3.2414411691881528d+00,
     +       5.7984074174475504d+00,       1.6397966768760199d+01,
     +       4.7382712899600312d+00,      -1.3678878183211181d+01,
     +      -2.6547564927279190d+01,      -1.0758757524101530d+02,
     +       7.3908325071927862d+01,       1.3736504967908880d+02,
     +      -1.2471666213371580d+02,       5.7590165882454833d+01,
     +      -3.3344191321691710d+01,       2.3144692822496811d+01,
     +      -1.5908182010256080d+01,       3.9930795315046091d-04,
     +      -3.0498663626806502d-04,       2.8930411200831332d-04,
     +      -2.6100051919846299d-05,      -4.2834917175003771d-04,
     +       1.1077970789631509d-03,      -1.8950864268602651d-03,
     +      -5.0786556219559950d-03,      -1.3243442040174721d-03,
     +       4.2085537661558864d-03,       7.4562070172753498d-03,
     +       3.4504645148461237d-02,      -2.4771682570609889d-02,
     +      -4.0144918242304130d-02,       3.7190662898859229d-02,
     +      -1.7349059412578340d-02,       1.0066250401665609d-02,
     +      -6.9637353385112906d-03,       4.7292445502375591d-03,
     +      -1.6883161568906021d+11,       2.4729101946471780d+11,
     +      -5.4916611363988083d+10,       7.3106343992295837d+10,
     +      -1.2125782378551790d+10,      -9.3864779781093102d+09,
     +      -3.2566184996104839d+10,      -2.7258856127810810d+10,
     +      -2.4124584334941261d+10,      -2.1268998755299179d+10,
     +      -1.2873189550831980d+10,      -1.5898122047574789d+10,
     +       2.8501418451574421d+10,       4.2489482158893356d+10,
     +      -2.3568224414048000d+10,       1.6881243754295050d+10,
     +      -8.3995769810467434d+09,       4.3652500992089729d+09,
     +      -1.7251487112335100d+09,       1.4338541574031749d+08,
     +      -1.9198014145975840d+08,       4.4429456904980481d+07,
     +      -6.2665225463142619d+07,      -5.4141122784422478d+05,
     +      -2.5661534264763952d+06,       1.6365169125635870d+07,
     +       1.9514788496354029d+07,       2.4037874367398862d+07,
     +       2.0442063442522489d+07,       1.0546408710297370d+07,
     +       3.5765371695188828d+07,      -3.1225794614777360d+07,
     +      -2.9821389520876881d+07,       1.4035050785006670d+07,
     +      -1.3916752759787681d+07,       6.9114050364853563d+06,
     +      -3.7226572908046148d+06,       1.0295964806737110d+06,
     +      -1.9128428069992060d+04,       2.6730450809121190d+04,
     +       1.1860288976942099d+03,       1.2499136420235831d+04,
     +       4.9954715780975521d+03,       2.5123186243948869d+01,
     +       1.8269018665513150d+02,      -1.2633114031392290d+03,
     +      -2.7617065876855679d+04,      -8.3932428270690907d+03,
     +       1.7323536285938690d+03,       3.4196305374543583d+04,
     +      -6.4153421499921868d+04,       5.0682392519175279d+04,
     +      -1.6968513312670209d+04,       7.3700138509977587d+03,
     +      -3.7452729567137412d+03,       1.9668589330434850d+03,
     +      -3.0641008156560929d+01,       9.7164995798283138d-01,
     +      -1.8912586167531591d+00,      -4.2124121884132144d+00,
     +      -2.6169143313061500d+00,      -2.5995114698449902d+00,
     +       2.7023071307286348d+00,      -2.5904945243742978d+00,
     +      -7.0448062761164554d+00,       2.6697976213123621d+01,
     +       2.5784172429373831d+00,      -6.6677341492901876d+00,
     +      -6.8685678390126540d+01,       1.0095167962462079d+02,
     +      -5.2487492270838153d+01,       1.7940756208755570d+01,
     +      -4.1961607159860277d+00,       2.2313614011593779d+00,
     +      -1.1459376740976130d+00,      -2.2129950998939221d-01,
     +      -2.3097772890458030d-04,       1.8506077681455489d-04,
     +       9.9687462984646103d-04,       6.0953087865203335d-04,
     +       6.5936522183748813d-04,      -1.2178569383425450d-03,
     +       9.5070748946486216d-04,       2.9677074533704649d-03,
     +      -8.7396953189459829d-03,      -1.3423786468604859d-04,
     +       2.7785151597792389d-03,       2.5859768796088391d-02,
     +      -3.6667258868963860d-02,       1.6950716021430640d-02,
     +      -5.6664328018188291d-03,       9.4798205101135843d-04,
     +      -5.2222177895849971d-04,       2.5762292103442792d-04,
     +       1.1391229181519770d-04,       3.7335178624469330d+10,
     +      -3.7329625832492188d+10,       1.5047551886685480d+10,
     +      -1.4014953241277010d+10,      -2.0057910712845380d+09,
     +      -6.3555387860173397d+09,      -1.6191226305325711d+09,
     +      -9.7225609058066964d+07,       9.4534563499712253d+08,
     +       9.4405345484513998d+08,      -1.0055900849329931d+09,
     +       5.3121139255935802d+09,      -9.7429519832618785d+08,
     +       4.9115618362625837d+09,       7.3471549094778800d+08,
     +      -2.6413756642315998d+09,       1.5994147977531450d+09,
     +      -1.0553988196540580d+09,       2.1996408319185480d+08,
     +      -2.0298338906758498d+07,       1.8610182135892160d+07,
     +      -8.5725505505612977d+06,       6.3552466623605946d+06,
     +      -7.3390020572172198d+05,       5.5911839675066741d+06,
     +       6.3711322373438790d+05,      -4.9599133589749411d+06,
     +       2.0963284553444050d+06,       1.1215783896620110d+07,
     +       1.0575785235569220d+07,       1.5779393023871221d+07,
     +      -3.8052040654503748d+07,       8.5871699949676692d+06,
     +      -2.0819451835103940d+07,       2.1408376948444929d+07,
     +      -1.1208481321207590d+07,       4.9223482075089663d+06,
     +      -7.2011804458003817d+05,       7.0792170659004123d+03,
     +      -4.8076914242405164d+03,       2.5604838161771190d+03,
     +      -1.1499250423486360d+03,       4.1229182657057017d+03,
     +      -6.0740952144279909d+03,       3.0583824076107630d+03,
     +       1.3715259111515121d+04,      -4.2041617050560144d+03,
     +      -2.8851745000756699d+04,      -1.6392438316088421d+04,
     +      -6.8109641018680588d+04,       1.2349131233200269d+05,
     +      -5.1464611072090840d+04,       6.0085013742576732d+04,
     +      -5.0825041681851748d+04,       2.6136072803989158d+04,
     +      -1.1539146462638280d+04,       2.9869846112417158d+03,
     +      -3.1261067151520070d+00,       1.7497704667628009d+00,
     +      -6.1404182443589694d-01,       4.9871870156207909d-01,
     +      -3.9779407256752899d+00,       4.8920112281563632d+00,
     +      -4.1736996025767539d+00,      -1.2930754589958539d+01,
     +       3.5744884921614859d+00,       2.5994950074054390d+01,
     +       1.0322634574990021d+01,       7.3724371640991734d+01,
     +      -1.2469180966672820d+02,       5.6776287203979273d+01,
     +      -5.7872802298798923d+01,       4.5948162714681729d+01,
     +      -2.3567387930303578d+01,       1.0974017743177800d+01,
     +      -4.0768947271004139d+00,       7.5438935501575893d-04,
     +      -4.0013372547030189d-04,       7.8318055310512089d-05,
     +      -1.7069891283723641d-04,       1.2053315165583649d-03,
     +      -1.4533098101821740d-03,       1.4222587501526169d-03,
     +       3.9578162280421757d-03,      -1.1001502689071011d-03,
     +      -7.7926771402069368d-03,      -2.3603081586059710d-03,
     +      -2.3704720926921639d-02,       3.8964414362466199d-02,
     +      -1.8228384697860551d-02,       1.7751187212462171d-02,
     +      -1.3742412407783089d-02,       7.0545128725467742d-03,
     +      -3.4369377461975080d-03,       1.5439512540532520d-03,
     +       1.3120304921090161d+12,      -1.3532960524465330d+12,
     +       5.1758081516423962d+11,      -5.1070653800583893d+11,
     +      -6.7603920515636124d+10,      -1.9591210676764679d+11,
     +      -3.7432620126539124d+10,       4.5472354512918110d+09,
     +       6.1991156451527428d+10,       7.5660305551095352d+10,
     +       6.2526396457903900d+10,       9.9241466938730438d+10,
     +       4.2534055867806039d+09,       2.2493418553448311d+10,
     +       1.3654775626519739d+10,      -5.1614712005963697d+09,
     +       1.0967575208378830d+09,      -5.2783292714068823d+09,
     +      -1.5194491557511251d+09,      -4.5409392467688668d+08,
     +       4.8928893917845619d+08,      -1.6902948623570469d+08,
     +       1.8253663312952641d+08,       1.9190050504583310d+07,
     +       5.4803286779371269d+07,       1.1069188959777599d+07,
     +      -3.8234411368042238d+07,      -1.0249835687379900d+08,
     +      -4.9078305744720012d+07,       3.2666801343651060d+07,
     +       1.8421350690769181d+08,      -1.4321472220515451d+08,
     +      -7.0907764193839267d+07,       1.0874922930018230d+08,
     +      -8.7863261796947554d+07,       5.3268573579067387d+07,
     +      -3.5401803319702253d+07,       2.1394030321976140d+07,
     +       2.2430107133742738d+04,      -2.1021099649367250d+04,
     +      -1.4006533359919081d+02,      -5.0521549744194162d+03,
     +      -1.7807873390920260d+04,      -9.5765782328572859d+03,
     +      -1.8431571539188499d+04,       7.8384800713096811d+03,
     +       2.9964975255650718d+05,       5.5437837777264693d+04,
     +      -2.6260115213826450d+05,      -4.0383731817690592d+05,
     +       2.1779246194586990d+05,       3.1678978132935811d+05,
     +      -2.8023187376775412d+05,       1.4927199525696540d+05,
     +      -1.0071224728298229d+05,       7.4760523820294082d+04,
     +      -2.9577082504718179d+04,       7.6652545112849788d+00,
     +      -1.4838639862359560d+01,       5.7402362880637892d+00,
     +      -7.5855776785836042d+00,       1.9119789375256801d+01,
     +       9.8468464655577606d+00,       8.9971253180751880d+00,
     +       1.2865146993591560d+01,      -3.2297288285829927d+02,
     +      -2.7420066786429789d+01,       3.2944635610308501d+02,
     +       3.5711064675250208d+02,      -2.4362206569997670d+02,
     +      -3.0030921872344987d+02,       2.2058973134185700d+02,
     +      -8.0759989321942257d+01,       6.8134979996455485d+01,
     +      -5.5500002105490061d+01,       9.9083650710643685d+00,
     +       3.6164881609297771d-04,       1.8920427720855290d-03,
     +       1.5506314179697779d-04,       1.1335855415494081d-03,
     +      -6.2701055594999401d-03,      -3.4207588007970908d-03,
     +      -1.3055074139065161d-03,      -5.2404559264181113d-03,
     +       1.0611110305489660d-01,       1.6619541771755589d-03,
     +      -1.1342737411114800d-01,      -1.1550791468367470d-01,
     +       9.7271432577261355d-02,       8.4373469124029016d-02,
     +      -5.8918705737253628d-02,       1.6077768022783719d-02,
     +      -1.7428294331067840d-02,       1.5143523307849261d-02,
     +      -1.1309677023463459d-04,       1.1588303827466389d+12,
     +      -1.2160756996665359d+12,       4.0886921991435632d+11,
     +      -5.9772561745154224d+11,      -1.6588405071554681d+11,
     +      -2.3154910992965079d+11,      -1.3803133988925890d+10,
     +       5.0607105087063961d+09,       2.0854380192742450d+11,
     +       4.9692227162647668d+11,       9.1533831050507088d+09,
     +       1.1891668310955759d+12,      -1.8025884380342041d+12,
     +       8.9476432781516553d+11,      -1.9129871773647760d+11,
     +      -2.6336217108069891d+11,       2.1028378810315881d+11,
     +      -1.5365210946509839d+11,       5.7265157307193916d+10,
     +      -2.6221760819411731d+09,       2.4389221017114272d+09,
     +      -1.0857881366331980d+09,       1.1058383944972751d+09,
     +       3.3583318397644383d+08,       6.8896392837522376d+08,
     +       1.9377390830191889d+08,       2.0890293073252550d+08,
     +       6.0971942102078423d+07,      -8.8758183429005444d+08,
     +       1.2494238688444121d+09,      -5.3131550923026390d+09,
     +       6.5294357021154642d+09,      -3.3390645723259120d+09,
     +      -1.7262298035524690d+09,       3.5361657039392958d+09,
     +      -2.3888403154834962d+09,       1.6880866145516231d+09,
     +      -9.4368768886899078d+08,       1.1183932957147460d+06,
     +      -9.6467404164481815d+05,       4.8476365705501742d+05,
     +      -3.7550995492494089d+05,      -1.0585012210119770d+05,
     +      -1.8068399771180499d+05,      -3.2348254630012170d+03,
     +      -1.5411615136752260d+05,      -2.0588266651589211d+05,
     +       6.0400504507782264d+05,      -3.5272727403011709d+06,
     +       4.0973757687712582d+06,      -3.6020438753002789d+06,
     +       1.6034478653283359d+06,       6.4660820967768589d+06,
     +      -8.4198291792822201d+06,       5.4672718144912580d+06,
     +      -3.7756627904796302d+06,       2.2138839497090941d+06,
     +      -2.1778821954349800d+02,       1.4217736623386361d+02,
     +      -9.9586115725407197d+01,       4.3173105088291173d+01,
     +       1.2506823134682669d+01,      -1.8193745524987941d+02,
     +      -1.7730120638601989d+01,       1.5510638954864530d+02,
     +      -2.2439343637762201d+02,      -4.0591622447335021d+02,
     +       3.5313982608311890d+03,       6.3361013028288561d+01,
     +      -1.9025600413615509d+03,       1.0121421178573710d+03,
     +      -6.6130458088845262d+03,       7.4966846404874732d+03,
     +      -4.8280266275772756d+03,       3.3596307937227871d+03,
     +      -2.0321764685090100d+03,       2.5137946412373401d-02,
     +      -9.9978671194593336d-03,       1.1899959579294010d-02,
     +      -6.9974923945251878d-03,      -2.2690721630361749d-03,
     +       9.4200092899255405d-02,      -4.0210274837006914d-03,
     +      -6.2150297680624522d-02,       1.5749264188706030d-01,
     +       1.4179264399027300d-01,      -1.1175947105050339d+00,
     +      -6.7757031173289006d-01,       1.4486407908275720d+00,
     +      -7.1760515504767641d-01,       2.1514185044377521d+00,
     +      -2.2781491779514309d+00,       1.4611219333241370d+00,
     +      -1.0215031421632099d+00,       6.2761608691461079d-01,
     +      -1.0265172134176440d+11,       9.0715136685512772d+10,
     +      -4.7365695026638046d+10,       3.8297678065899963d+10,
     +       1.0713996079260950d+10,       3.0107368197081879d+10,
     +       1.3977503336947430d+10,       7.4083576193561411d+09,
     +       2.3802466809597812d+09,      -6.7152280038394279d+09,
     +       5.1375975653983955d+09,      -1.0364103923260690d+11,
     +       1.1446798188221809d+11,      -4.6381317792895027d+10,
     +      -1.9701454595904861d+10,       2.7968816689434700d+10,
     +      -1.7368421420946350d+10,       1.0518896068766081d+10,
     +      -3.9429100905379019d+09,       5.0730312652907021d+07,
     +      -4.3417962821323998d+07,       2.4474027421306521d+07,
     +      -1.5156614384289831d+07,       2.3618640216186638d+06,
     +      -3.7357164369533777d+07,      -7.2817514620917113d+06,
     +       4.4309053730826654d+07,      -1.5108833180587569d+07,
     +      -4.7201278512674756d+07,      -2.1353621107528850d+07,
     +       1.5493632238454959d+08,      -1.2553082183556870d+08,
     +      -4.7980737729690149d+07,       1.4911347967683151d+08,
     +      -1.2379547239133470d+08,       4.9837004101451449d+07,
     +       6.2625100038239947d+06,      -2.2686186693305589d+07,
     +      -8.1898254880441827d+03,       3.5176462741158762d+03,
     +      -2.8158245438479421d+03,      -1.0217972079692569d+02,
     +      -1.8888712210552501d+04,       6.7000798901768168d+04,
     +      -9.4693018835481635d+03,      -1.3266212814646031d+05,
     +       3.4972313784131613d+04,       1.5811602314045621d+05,
     +       1.5901513503875160d+04,      -1.1188239543751770d+05,
     +      -4.2580387557285067d+04,       2.7184329719084268d+05,
     +      -3.7696392545946519d+05,       2.8592861473885528d+05,
     +      -9.1959475395966205d+04,      -7.5370342675241423d+04,
     +       1.2678977978080540d+05,      -1.1115919851775871d+00,
     +       3.4807003372005019d+00,      -2.7603191339331330d+00,
     +       3.2428656853666049d-01,       1.9508843183906620d+01,
     +      -6.3755740611582347d+01,       1.7735223996311401d+01,
     +       1.3457108430423051d+02,      -3.7524315422928836d+01,
     +      -1.6812415010740230d+02,       1.7547898642526980d+00,
     +       1.8395233297921600d+01,       1.5993250081277961d+02,
     +      -3.0172715731581587d+02,       3.6472821687047337d+02,
     +      -2.7258923317747400d+02,       7.9370450873832354d+01,
     +       9.9012073716871683d+01,      -1.6107534291305640d+02,
     +       6.5967927466586951d-04,      -1.1876400756248070d-03,
     +       1.1382029565048641d-03,       1.9411288973155460d-04,
     +      -6.1676925344941232d-03,       2.0008971490330551d-02,
     +      -6.6550623551158783d-03,      -4.3253203242665982d-02,
     +       1.2533719870967200d-02,       5.5333671172380750d-02,
     +      -2.7581915981067791d-03,       5.6247315093746961d-03,
     +      -6.5211348463510205d-02,       9.8227238647893234d-02,
     +      -1.1424216747088610d-01,       8.5550706012484670d-02,
     +      -2.3905270038487549d-02,      -3.5020763022420943d-02,
     +       5.7034078216254759d-02,      -1.8409519176695270d+11,
     +       1.6024365616603751d+11,      -9.5313088530156174d+10,
     +       6.9230779763144012d+10,       3.1018202510261162d+10,
     +       8.9905317408216904d+10,       4.1338818179212089d+09,
     +       7.1332357880772686d+09,       1.5676675087913401d+10,
     +      -4.2099015019027267d+10,       4.7069419426584427d+10,
     +      -1.9538195686544220d+11,       4.3740349626072540d+11,
     +      -5.5872055043598181d+11,       2.1375484580313541d+11,
     +       5.7287263472845642d+10,      -7.8089595143170685d+10,
     +       3.1600820560085217d+11,      -7.0559356956092944d+11,
     +       2.0568584283868501d+08,      -1.7648766472418571d+08,
     +       9.5305144903975397d+07,      -5.2519898949063681d+07,
     +      -1.6410083370078770d+07,      -2.2641850327161270d+08,
     +       1.6470604569616279d+08,       3.0811093730930507d+08,
     +      -1.7271989491219240d+08,      -5.4815093963064086d+08,
     +      -3.3654485598234397d+08,      -4.8544754981056070d+08,
     +       1.0876820700154430d+08,       1.5970590212318790d+09,
     +       6.0889562212163307d+07,      -1.3235149726867919d+09,
     +       9.4070444002299333d+08,      -2.8015622396831760d+09,
     +       6.3925308268955507d+09,      -4.6355473042263933d+04,
     +       3.8695319498279147d+04,       2.7465781435193989d+04,
     +      -1.8488141472405070d+04,      -4.0564807670030292d+04,
     +       3.8059263971834362d+05,      -4.6289962155621487d+05,
     +      -9.0854894952461764d+05,       2.2382116367315699d+05,
     +       1.7867267841803459d+06,       7.2236753339684347d+05,
     +       2.5555954273203360d+06,      -2.5880348644260569d+06,
     +      -1.9428087873014901d+06,      -1.2679816770710060d+06,
     +       3.6911781762249982d+06,      -2.6066899378083390d+06,
     +       7.6714644732411057d+06,      -1.7387484617971878d+07,
     +      -1.9051488452517241d+01,       1.7470416746397959d+01,
     +      -6.3406277361751350d+01,       2.6543358792739610d+01,
     +       4.7262625217187448d+01,      -3.2932648012661213d+02,
     +       4.2361673767544789d+02,       8.6844831404240665d+02,
     +      -8.3339910701613121d+01,      -1.7591415087359710d+03,
     +      -6.4022324124959141d+02,      -2.9141066378265050d+03,
     +       3.2176873472873449d+03,       1.2120498487567361d+03,
     +       1.4702298244657650d+03,      -3.5210223228956679d+03,
     +       2.5797080651731908d+03,      -7.7430352943769976d+03,
     +       1.7404243658070329d+04,       8.8864865761986451d-03,
     +      -7.8242398498195812d-03,       2.1223146136109362d-02,
     +      -6.7536279418750561d-03,      -1.4561355183517611d-02,
     +       9.7841192063414517d-02,      -1.2584754122524860d-01,
     +      -2.6679926709019680d-01,       5.6693505090348551d-03,
     +       5.4719500870374915d-01,       1.8685652037253780d-01,
     +       9.8263668873410548d-01,      -1.1129462711689111d+00,
     +      -2.8288101900185669d-01,      -4.8603129903041742d-01,
     +       1.0821006217170590d+00,      -8.1708778492382084d-01,
     +       2.4887906871357588d+00,      -5.5601497502720676d+00,
     +      -3.4501503790290669d+12,       3.1620466040249268d+12,
     +      -1.7162049009094890d+12,       1.5535656859543069d+12,
     +       7.0171374954429846d+11,       1.3262983206992051d+12,
     +       4.7029129366976477d+11,       5.8263217799162720d+11,
     +      -2.3095556586665689d+11,      -2.3959194790648262d+12,
     +      -4.3245459968568610d+11,      -5.6109383154183105d+12,
     +       1.0976792105420051d+13,      -9.3109750414882949d+12,
     +       6.7895128080129365d+12,      -3.8877882146977788d+12,
     +       2.0357636849324399d+12,      -1.3615731346353501d+12,
     +       1.7493406539080869d+12,       8.2537479089391413d+09,
     +      -7.2822007921648836d+09,       4.6574331855641441d+09,
     +      -3.0547518146902232d+09,      -1.9364077635601590d+09,
     +      -3.2516609280001001d+09,      -1.8044657630927460d+09,
     +      -1.5338694337721350d+09,       5.8376423156078708d+08,
     +       8.1476416157264690d+09,       3.0702950254511487d+08,
     +       8.3930112097128592d+09,      -2.7271054263784882d+10,
     +       3.9420788908910454d+10,      -4.6818305640246490d+10,
     +       3.7899887048452988d+10,      -2.3412525139011669d+10,
     +       1.8944698340323341d+10,      -2.0441094417270908d+10,
     +      -3.8537773798836689d+06,       3.5379041112445332d+06,
     +      -2.7244862836253871d+06,       1.2927707139417981d+06,
     +       2.3615541605266160d+06,       1.3088724744689000d+05,
     +       3.3382535389902992d+06,       1.6994517279869949d+06,
     +      -7.2261600834995378d+06,      -1.5837370884823130d+07,
     +       7.7543417863419149d+06,       2.1466106026418220d+07,
     +      -3.7337803073973679d+06,      -4.4221250854700021d+07,
     +       8.0022251564292014d+07,      -7.4654873397246897d+07,
     +       4.4048141496544182d+07,      -3.2323033923919410d+07,
     +       3.5447492929459862d+07,       7.8158078087844831d+02,
     +      -7.4251332259371702d+02,       1.0700252504048799d+03,
     +      -2.2928573392728629d+02,      -2.0580910730430910d+03,
     +       2.0512407140206542d+03,      -3.6504685152873972d+03,
     +      -2.7713451628883290d+03,       1.1425887614646170d+04,
     +       1.5858750766450870d+04,      -1.3650490453290609d+04,
     +      -3.5012239512907661d+04,       2.9255548729853159d+04,
     +       2.2660508063452729d+04,      -5.8652224193136914d+04,
     +       5.8810128649864237d+04,      -3.0653587318006681d+04,
     +       1.6828468737826832d+04,      -1.8977731713724621d+04,
     +      -8.0867646891416775d-02,       7.7066592141498039d-02,
     +      -2.5094354437037009d-01,       3.2547301355050337d-02,
     +       6.6043252474359226d-01,      -8.9538021821873170d-01,
     +       1.2707805041716180d+00,       1.1240596948119199d+00,
     +      -4.3375626300778398d+00,      -5.1401387948535593d+00,
     +       5.5518424685123664d+00,       1.2639580791736220d+01,
     +      -1.2831113769880870d+01,      -4.0113116111408171d+00,
     +       1.5012223955136850d+01,      -1.5762382810387139d+01,
     +       7.1598370004483556d+00,      -2.1586544720820071d+00,
     +       2.4209165427062831d+00,      -2.3367820128506830d+11,
     +       2.8901360312376300d+11,      -1.2225610320475999d+11,
     +       1.0645860706190790d+11,       1.5054157203346251d+10,
     +       4.1342658672389252d+10,      -2.0869575277427311d+10,
     +       2.5809814531291290d+10,      -1.9262232327124439d+10,
     +      -1.5726050116433020d+11,       1.1794444039230920d+11,
     +      -4.1692513375427148d+11,       6.6402408013540466d+11,
     +      -5.8078677849621912d+11,       4.9873702983364941d+11,
     +      -3.1510170125204242d+11,       1.8062235432300110d+11,
     +      -1.0020759452912810d+11,       6.3395250331668854d+10,
     +       1.0817434941792800d+09,      -1.2926254627179010d+09,
     +       6.2128462190211856d+08,      -4.6351160753729713d+08,
     +      -1.0175072089161240d+08,      -1.1682182507148890d+08,
     +       6.8675936569158152d+07,      -2.0631328824378639d+08,
     +      -1.9800038256033700d+07,       1.0340632516903120d+09,
     +      -1.4889396676130259d+09,       2.6341356245018010d+09,
     +      -3.5663210326939330d+09,       3.4508974334825578d+09,
     +      -2.5079961239537530d+09,       1.3200548438463359d+09,
     +      -7.4251687141005278d+08,       4.6297071356168181d+08,
     +      -5.1816236528589797d+08,      -9.3011292995226558d+05,
     +       1.1375622294021361d+06,      -5.9419015338317642d+05,
     +       4.3354190111519961d+05,       9.7986242264452041d+04,
     +      -1.4568019854852519d+05,       5.7859185989726677d+04,
     +       8.6388307587176445d+04,       1.5554245690820221d+05,
     +      -1.8956972637265881d+06,       4.4226566806523316d+06,
     +      -5.6197687663999107d+06,       5.7844110550118666d+06,
     +      -5.1369513578912793d+06,       2.6692206901199091d+06,
     +      -7.1795397112773149d+05,       3.9727615454096178d+05,
     +      -7.7062630830866483d+05,       2.1494743655467578d+06,
     +       3.2071956280738908d+02,      -4.1588515168464761d+02,
     +       2.5203790877477141d+02,      -1.9931911900201860d+02,
     +      -4.6276495833290838d+01,       2.9850432670920628d+02,
     +      -1.5436247971268679d+02,      -3.1272667117151451d+00,
     +       2.1470815063355449d+01,       1.6344124360885969d+03,
     +      -4.6114932476206168d+03,       4.8667186101729094d+03,
     +      -4.0867571588043211d+03,       3.2911786837287832d+03,
     +      -9.8389818642893533d+02,      -3.7622401287021131d+02,
     +       1.8910961104600000d+02,       7.8027119552523152d+02,
     +      -2.8654596575805199d+03,      -4.6265027919358472d-02,
     +       6.5282850460895728d-02,      -4.8364044811677422d-02,
     +       4.4416363102608561d-02,       9.8853341438880551d-03,
     +      -1.1173341702056490d-01,       6.2882210671263261d-02,
     +       1.3665731598503920d-02,      -5.6732364134261137d-02,
     +      -5.2225838321955975d-01,       1.5052576601163770d+00,
     +      -1.4182787386272091d+00,       1.0656779021658480d+00,
     +      -8.0500319240117435d-01,       9.1052947215079727d-02,
     +       2.6245344630716211d-01,      -1.3230386538825120d-01,
     +      -2.8229811451622577d-01,       1.0830514738937529d+00,
     +       2.9014561374009861d+10,      -3.3871353925254841d+10,
     +       1.6130872642664680d+10,      -1.1427562930880230d+10,
     +      -1.9524187624361949d+09,      -8.6539195832810402d+09,
     +       6.4085569796936836d+09,      -8.8218552998031998d+08,
     +      -6.6182737641932964d+09,       5.0074654393351774d+09,
     +       1.3451138037786160d+09,       3.8922763220575844d+10,
     +      -6.3000021112254333d+10,       6.2431744397663887d+10,
     +      -5.4454215180884644d+10,       3.4224509591195122d+10,
     +      -2.2129509949209808d+10,       2.5109160373755089d+10,
     +      -4.0609089573753738d+10,      -3.7384122165969491d+07,
     +       4.3105528891175970d+07,      -1.9043416165659651d+07,
     +       1.1130940403044410d+07,      -1.0066257744741270d+06,
     +       3.8056429246174969d+07,      -3.8651343536309019d+07,
     +      -3.7443813564966217d+07,       5.3168407993727408d+07,
     +       6.4030524145913847d+07,      -5.7444242803215362d+07,
     +      -5.5043506393009193d+07,       7.0261655071956307d+07,
     +      -9.4682418925741509d+07,       1.0757728921986879d+08,
     +      -7.5351304109930202d+07,       6.3270466114296667d+07,
     +      -1.6658958970437470d+08,       3.9391367497814298d+08,
     +       1.4871496082453470d+04,      -1.4645386166586901d+04,
     +       9.2157135546279778d+02,       2.8266800586724880d+03,
     +       1.2511198515656250d+04,      -8.2906738852917202d+04,
     +       8.4214161453064880d+04,       1.0589499536010600d+05,
     +      -9.7292488993963721d+04,      -2.1860776514758001d+05,
     +       1.7835258110658071d+05,      -2.0642733351734780d+04,
     +       4.7841087599034050d+04,       4.5041595589757810d+04,
     +      -1.2395382422011950d+05,       1.0883228528126820d+05,
     +      -1.1479683174306130d+05,       4.2550076103948691d+05,
     +      -1.0911956816745950d+06,      -1.7603342206447130d+00,
     +      -1.1817193210837491d+00,       6.8816773755069907d+00,
     +      -6.5997645215434053d+00,      -1.4499832599625581d+01,
     +       7.5114503222015273d+01,      -7.2129583738672991d+01,
     +      -9.8934883317906383d+01,       6.9199611605740401d+01,
     +       2.2198130600200821d+02,      -1.7552276965330080d+02,
     +       6.5799442296546502d+01,      -1.0287511936816190d+02,
     +      -1.7964835474382279d+00,       9.1064087372266400d+01,
     +      -9.1562964940605085d+01,       1.0269437039892181d+02,
     +      -4.1125714785844940d+02,       1.0776942981941711d+03,
     +      -1.8353467625992421d-04,       1.2862803643005150d-03,
     +      -2.6231495322517739d-03,       2.1278757038595518d-03,
     +       4.8150277075363363d-03,      -2.2595805753956980d-02,
     +       2.0602174650409821d-02,       3.0178169739667929d-02,
     +      -1.7329255078393582d-02,      -7.0328511281420994d-02,
     +       5.4296957986763851d-02,      -2.5783503426440611d-02,
     +       3.8435318840949549d-02,      -3.4548179746836838d-03,
     +      -2.5847505551750988d-02,       2.7489304809298239d-02,
     +      -3.1267327007512169d-02,       1.2785366211398300d-01,
     +      -3.3801701519056337d-01,       4.8484183600653198d+11,
     +      -5.6036053345656934d+11,       2.6723227613905469d+11,
     +      -1.9731550184493991d+11,      -3.2324946722468109d+10,
     +      -1.0558654685013969d+11,       4.3422887314292702d+10,
     +      -1.9252216894950211d+10,       2.3080392345594250d+10,
     +       8.2940735655721252d+10,      -4.2547607643757973d+10,
     +      -3.8604570350143356d+10,      -1.1178923080711330d+11,
     +       4.3910199897540588d+11,      -3.7278673447968939d+11,
     +       2.2295288594808929d+11,      -1.3763504234997421d+11,
     +       7.4156321315833252d+10,      -1.9152146773832008d+10,
     +      -3.5239343253700542d+08,       4.4592932259822679d+08,
     +      -2.0563394736713940d+08,       1.1952167168816800d+08,
     +       7.2752675899084806d+07,       5.2959258918429397d+07,
     +      -1.2578243012688270d+08,       3.9465096154599316d+07,
     +      -7.1299519558295161d+07,      -2.1026748959682949d+07,
     +       5.4882786014162749d+06,       2.5970611661585188d+09,
     +      -3.2796962755927391d+09,       9.5548758093391263d+08,
     +      -3.2316158175018239d+08,       1.0936859849582730d+08,
     +       9.7933525063450024d+07,      -2.8655451458178413d+08,
     +      -2.0294853516490371d+06,       6.2885458710563151d+04,
     +      -1.2109903780193950d+05,       5.2190446178881422d+04,
     +       2.1155223142101459d+04,      -1.0278072501329931d+05,
     +       3.1942556765391379d+04,       2.2779550130587039d+05,
     +      -3.7808062097173708d+05,       6.5535432916679431d+05,
     +      -3.5130384372055082d+05,      -5.6771230397573195d+05,
     +      -5.7099695430265209d+06,       8.6614425108299889d+06,
     +      -3.0691912421137439d+06,       2.8115033606826782d+05,
     +       5.9794561860017991d+05,      -8.5071575442890148d+05,
     +       1.1132432196477400d+06,       4.9858963205840980d+05,
     +       2.3099384538105792d+01,      -4.3057956714110161d-02,
     +       8.4991105714415816d-01,      -3.4144581542177960d+01,
     +       5.9298675559441378d+01,      -4.3089837608754522d+01,
     +      -1.8030168627212299d+02,       4.7117696526553982d+02,
     +      -8.7101125603335606d+02,       4.9425423762427499d+02,
     +       8.1243647169405540d+02,       4.9088302187947384d+03,
     +      -8.0830935961842752d+03,       2.8703771465887030d+03,
     +       1.2519269638222420d+02,      -9.9013040278492338d+02,
     +       1.0688167809793210d+03,      -1.1054302473820831d+03,
     +      -1.0117642426111090d+03,      -9.3849963425350607d-03,
     +       5.0577179616644987d-03,      -2.1023245010675920d-03,
     +       7.4961445735357681d-03,      -1.2453321571990320d-02,
     +       1.3562439882652041d-02,       4.9606892371508643d-02,
     +      -1.5784054849822951d-01,       3.1044392055973369d-01,
     +      -1.8917250507364089d-01,      -2.8550036038351689d-01,
     +      -1.4339150271746139d+00,       2.4690100792183021d+00,
     +      -8.7862902575123636d-01,      -9.5657023275659198d-02,
     +       3.6689709511757479d-01,      -3.6541886559806469d-01,
     +       3.3496674868995641d-01,       4.2653941060286310d-01,
     +       5.0850032557564323d+10,      -6.2924884563718369d+10,
     +       3.3076563911866711d+10,      -1.6536927465742241d+10,
     +      -9.8549521068981247d+09,      -2.7868721132554668d+10,
     +       3.3209740292465450d+10,      -1.8739590142898830d+10,
     +       4.2892212811127982d+09,       8.7181015019718857d+10,
     +      -6.1331801095646606d+10,      -1.4304747708531900d+11,
     +       3.0137621902601780d+10,       1.8335887885467459d+11,
     +      -8.6853273921456635d+10,      -2.7212678963522062d+09,
     +       5.6485839066305397d+10,      -2.5122163654348111d+11,
     +       4.3813825579677380d+11,      -1.0734724489653669d+08,
     +       1.2808916107534400d+08,      -5.6652607792300537d+07,
     +       1.4171194862180630d+07,       1.5826246713121761d+07,
     +       1.3709603355269039d+08,      -2.4044288939852649d+08,
     +      -8.3379438501277149d+07,       9.8788615862176090d+07,
     +       2.7121195082897651d+08,       1.2830270874571210d+08,
     +       8.9364133900823045d+08,      -1.0902595941624050d+09,
     +      -2.5276453977818418d+07,      -6.5614761496515989d+08,
     +       1.0089157012353270d+09,      -1.1318616006077659d+09,
     +       2.9552030318523450d+09,      -5.5376361956437511d+09,
     +       6.7591288696014191d+04,      -6.9232362487053644d+04,
     +      -1.2534239584899770d+02,       2.2986495737283170d+04,
     +       2.2732332041612761d+04,      -2.7665733490344178d+05,
     +       5.3958972969644109d+05,       3.4361939595035481d+05,
     +      -1.5344997642973429d+05,      -1.3849728132189200d+06,
     +      -1.6331637975218779d+05,      -2.3308528231731271d+06,
     +       3.8755922142629302d+06,      -1.4545699901857961d+06,
     +       2.9055650397487348d+06,      -3.3454693374582068d+06,
     +       3.3412046801427868d+06,      -7.9742165704237437d+06,
     +       1.5401591397742480d+07,      -1.5681269589274580d+01,
     +       4.2485128101094434d+00,       3.3717773675952003d+01,
     +      -2.7592439182789679d+01,      -4.4548649275732792d+01,
     +       2.3799753107879741d+02,      -4.6657004064307807d+02,
     +      -3.3351269469735848d+02,       1.6152546268810209d+01,
     +       1.4850844673685810d+03,       2.1451056100770691d+02,
     +       2.2742482525019668d+03,      -4.1079241163720153d+03,
     +       1.8292193833920901d+03,      -3.0064005485219450d+03,
     +       3.2695357247382622d+03,      -3.2148957806161029d+03,
     +       7.6986523255814336d+03,      -1.5306063207898500d+04,
     +       7.2045675364839138d-04,       4.0079718829515651d-03,
     +      -1.3114551927973299d-02,       7.7648013105706907d-03,
     +       1.7091288068269480d-02,      -6.9671312287005044d-02,
     +       1.3463570428996971d-01,       9.8900801582472761d-02,
     +       2.1519621602175760d-02,      -4.6857492659529498d-01,
     +      -8.7787968352710286d-02,      -7.1455154814805355d-01,
     +       1.3269199651241299d+00,      -6.1369821060742025d-01,
     +       9.4119713805655270d-01,      -1.0031416935816919d+00,
     +       9.8616398578259246d-01,      -2.3955709671392591d+00,
     +       4.8579310642337798d+00,       6.9595468308182275d+11,
     +      -8.9151140198465723d+11,       4.7521930833090820d+11,
     +      -3.0005858168392798d+11,      -1.0479882659036430d+11,
     +      -2.1444852688367120d+11,      -5.2765186179476753d+07,
     +      -2.4493888106849771d+11,       3.0477467145835248d+11,
     +       2.0777883464062920d+12,      -1.2345125730509089d+12,
     +      -1.6437206279158170d+11,      -2.3787289731013960d+12,
     +       4.5302564479157754d+12,      -4.6243619087084424d+12,
     +       3.3716791242533672d+12,      -2.0946064427721260d+12,
     +       8.0976115689931152d+11,       6.4845441990393274d+11,
     +      -3.6047070204570460d+09,       4.5804403874137373d+09,
     +      -2.7394600522773700d+09,       1.3522082559905379d+09,
     +       1.0422636872798001d+09,       5.5822635855751181d+08,
     +       8.5496784182229567d+08,       9.9338773615098250d+08,
     +      -1.8751077565720999d+09,      -1.7208893149607079d+10,
     +       1.3761062531472349d+10,       1.0733622616609180d+10,
     +       4.7785822244661274d+09,      -3.4661439283557259d+10,
     +       3.9926669994775520d+10,      -2.9993476530785660d+10,
     +       1.7614985464967621d+10,      -3.8963400804715209d+09,
     +      -8.0639787988512497d+09,       3.5964683312670970d+06,
     +      -4.8938320004094951d+06,       3.3512285775568732d+06,
     +      -1.3864181210506619d+06,      -2.0728571797298430d+06,
     +       1.1145866059728181d+06,      -4.1128114071847079d+06,
     +       1.1360548796237169d+06,       3.1658798935595760d+06,
     +       3.5915454359706499d+07,      -3.8193006077196017d+07,
     +      -2.4040818592102930d+07,       4.1076352320403471d+06,
     +       6.3998119462963998d+07,      -7.7965484552266210d+07,
     +       5.7990993129812591d+07,      -2.6360524729323320d+07,
     +      -2.6103060134323880d+07,       6.2619430353865243d+07,
     +      -1.4867494806084519d+03,       2.2847741545976492d+03,
     +      -1.9235796054010041d+03,       6.8616486705716443d+02,
     +       1.8315209906513651d+03,      -2.1982944916685710d+03,
     +       4.6125834842319773d+03,      -1.1582218387066041d+03,
     +      -4.0428385373492069d+03,      -3.2370645925386671d+04,
     +       4.0393066652105968d+04,       1.8999695047852889d+04,
     +      -8.8540254291578949d+03,      -5.1322165359731589d+04,
     +       6.5047166469514959d+04,      -4.7996283421606058d+04,
     +       1.5096737381252749d+04,       5.1251286939121099d+04,
     +      -9.5570580692376097d+04,       2.5333983310958819d-01,
     +      -4.4923707948233382d-01,       4.5974135905503821d-01,
     +      -1.4936540275533439d-01,      -5.6025406703430891d-01,
     +       8.2759740437796014d-01,      -1.4814267232702589d+00,
     +       1.2078013244906990d-01,       1.5553748599485719d+00,
     +       1.0149826379656799d+01,      -1.3305424270702281d+01,
     +      -5.3868233106105912d+00,       3.2988457051495752d+00,
     +       1.4770523255912900d+01,      -1.9116250764355570d+01,
     +       1.4013552734280919d+01,      -2.9667213667112531d+00,
     +      -2.1277996994802351d+01,       3.7497981798866668d+01,
     +      -5.4178848393463753d+10,       6.5044569433725899d+10,
     +      -3.8381727555538460d+10,       1.7511816304184479d+10,
     +       6.4467799988746080d+09,       3.7480494649920197d+10,
     +      -1.7403175353770229d+10,      -1.1867343748916821d+10,
     +       2.7410200749055470d+10,      -7.3328997296483658d+10,
     +      -3.2300030045996518d+09,      -1.5527141799889450d+11,
     +       4.8429757158503247d+11,      -6.8896504655474731d+11,
     +       7.4055079973844397d+11,      -5.5460147654311389d+11,
     +       3.8683307376445941d+11,      -4.5415206793669958d+11,
     +       7.7739782878846667d+11,       4.6390986612643018d+07,
     +      -5.9699580869206421d+07,       3.6621793988557577d+07,
     +      -6.0793724859191673d+06,       2.5322756786876209d+07,
     +      -1.8079355584333360d+08,       3.3284481217876799d+07,
     +       3.6851806461549813d+08,      -3.0775733047290659d+08,
     +      -3.3100087347670358d+08,       6.1296035837803257d+08,
     +       7.1536464094684327d+08,      -2.1044812107595191d+09,
     +       3.1772191597985120d+09,      -3.8471799261937790d+09,
     +       3.0208496911179991d+09,      -2.2958602429902010d+09,
     +       3.6229544805219421d+09,      -7.2783335474861250d+09,
     +      -1.2364895205866929d+04,       1.7919840027981161d+04,
     +      -1.4064714085111809d+04,      -6.4979033843072184d+03,
     +      -1.0955050830891410d+05,       4.2575420952260168d+05,
     +      -4.7877240871595983d+04,      -1.0166537991429800d+06,
     +       6.4632767451027175d+05,       1.3708852677444031d+06,
     +      -1.5808222625662889d+06,      -1.6536604486388930d+06,
     +       4.3381371438833307d+06,      -6.8572565128941406d+06,
     +       8.6299700590887405d+06,      -6.8577243693107832d+06,
     +       5.3758027829798479d+06,      -9.5402196792228743d+06,
     +       2.0574715350689039d+07,      -1.6001193773051079d+00,
     +       1.9286089101932180d+00,       1.3084779538313911d+00,
     +       2.2182197213602439d+00,       1.2156932855546100d+02,
     +      -3.9954161248275102d+02,       3.5154140084272548d+01,
     +       1.0092585045595620d+03,      -5.6902265399242185d+02,
     +      -1.5031889447713340d+03,       1.4837703960562051d+03,
     +       1.5350353101542230d+03,      -3.7876795536637378d+03,
     +       6.0998662140218494d+03,      -7.7178682424297804d+03,
     +       6.1233275765185190d+03,      -4.8578046891102986d+03,
     +       9.1950708549731862d+03,      -2.0714820752113708d+04,
     +       1.1646778117805480d-03,      -1.7888018045445430d-03,
     +       4.9700699223112564d-04,       5.0206486742546575d-04,
     +      -3.9998461804177918d-02,       1.2269463348766800d-01,
     +      -8.8643841262875042d-03,      -3.2102780793125468d-01,
     +       1.6923251293728461d-01,       4.9692863434468981d-01,
     +      -4.5184509733598388d-01,      -4.7278923497257153d-01,
     +       1.1213139343116061d+00,      -1.8189678768204041d+00,
     +       2.2980378805093311d+00,      -1.8170060332195810d+00,
     +       1.4492313500934220d+00,      -2.8520586941602679d+00,
     +       6.6023421879494606d+00,       1.0417031531397070d+12,
     +      -1.4379095071916841d+12,       8.1074923229910144d+11,
     +      -3.5235317369280463d+11,      -5.1189494305784843d+10,
     +      -7.7305610551702197d+11,       3.3180573638160580d+10,
     +      -1.0359051425487111d+12,      -6.2229151523917801d+10,
     +       1.2577893423508699d+13,      -6.6509978080216904d+12,
     +      -1.0423853406608381d+13,      -7.4261207992134888d+11,
     +       2.1925377336689961d+13,      -2.9926571112577219d+13,
     +       2.6371680059999590d+13,      -2.3926333678432121d+13,
     +       4.4732372032594141d+13,      -9.2159829300731562d+13,
     +      -8.1388651768345156d+09,       1.1310158458657890d+10,
     +      -7.0795682624238491d+09,       2.6727398551377020d+09,
     +       1.6337750082902081d+09,       3.2345020521767292d+09,
     +       7.6367392129537082d+08,       1.3695211713701611d+10,
     +       1.0059664163162260d+09,      -1.6320909976131329d+11,
     +       1.2108982008976390d+11,       1.1922892349154660d+11,
     +      -1.6975322307084040d+10,      -2.1714902136745609d+11,
     +       2.6413569496363950d+11,      -2.2148752528515750d+11,
     +       2.1470060552446481d+11,      -4.2116033535394153d+11,
     +       8.1034612228013916d+11,       1.2148899406826111d+07,
     +      -1.7743436962910101d+07,       1.1525069049711401d+07,
     +      -3.8481167489322228d+06,      -5.2843241978989691d+06,
     +       6.0207355957460159d+05,      -5.5513085494743404d+05,
     +      -1.7206801657599580d+07,      -3.1650232464467071d+07,
     +       3.9357020103854382d+08,      -3.6185001580136800d+08,
     +      -1.9073247889162111d+08,       1.7834277818996090d+07,
     +       4.9367295348781788d+08,      -5.3629846684713262d+08,
     +       4.2043056787331581d+08,      -4.5121466221878701d+08,
     +       1.0248974176109070d+09,      -2.1658744859771051d+09,
     +      -7.0775057693832841d+03,       1.1050613099090100d+04,
     +      -7.3937131900789682d+03,       1.9794916212756809d+03,
     +       6.0201135590148242d+03,      -4.1857519051852896d+03,
     +      -3.4076694049532848d+03,       7.3217094092186944d+03,
     +       4.7131882915184353d+04,      -3.5148833939710638d+05,
     +       3.6097132436687790d+05,       1.1765117266140830d+05,
     +      -1.3832484625361159d+04,      -3.9983496846366121d+05,
     +       3.9456426183036569d+05,      -2.9009904579117667d+05,
     +       3.4831506112965639d+05,      -9.1870036902794009d+05,
     +       2.0948197416584629d+06,       1.5394229965922490d+00,
     +      -2.5594923796910130d+00,       1.7553285704214130d+00,
     +      -3.3796887007916021d-01,      -2.1000608911047220d+00,
     +       1.6954236586431470d+00,       2.1782510228908389d+00,
     +      -1.3908296406911771d+00,      -1.6999939348831301d+01,
     +       1.0496105508902880d+02,      -1.1346253954573440d+02,
     +      -2.7265967925664761d+01,       4.6183828821201578d+00,
     +       1.1181021434726610d+02,      -1.0394924935356040d+02,
     +       7.2933076535283249d+01,      -9.4539206759262655d+01,
     +       2.7706932052283730d+02,      -6.6880825277940266d+02,
     +      -1.6158810946352231d+11,       2.1129861082888751d+11,
     +      -1.1821633833839630d+11,       3.2228705193924000d+10,
     +       3.7227317217454422d+10,       1.6549284771856110d+11,
     +      -2.1871928481240610d+11,       1.4722719159204489d+11,
     +       2.0008578110059460d+10,      -8.7302033515273120d+11,
     +      -4.8848223474191528d+11,       2.6321132507516870d+12,
     +      -3.0111350495089319d+11,      -2.8826447632156421d+12,
     +       3.4954528207028091d+12,      -2.9682616128590571d+12,
     +       3.3468896543593740d+12,      -1.0504889223758580d+13,
     +       2.9301017950527102d+13,       3.7829490390539593d+08,
     +      -4.7418752229732883d+08,       2.1551897057617050d+08,
     +       4.6906189978345737d+07,      -1.1648057063675950d+08,
     +      -9.6929334303785610d+08,       1.5266815570032420d+09,
     +       5.5291696532594657d+08,      -4.2260435269469649d+08,
     +      -1.9746524854570301d+09,       6.1644775783072376d+09,
     +      -8.7543106279599819d+09,      -2.6035816819208651d+09,
     +       1.5700045254225229d+10,      -1.6360068877862089d+10,
     +       1.2932349724699930d+10,      -2.0912371147381809d+10,
     +       9.8181370932583862d+10,      -3.0617165484119202d+11,
     +      -2.6334295774386491d+05,       2.8405299604893738d+05,
     +      -1.4743734652077350d+04,      -2.3903278619311389d+05,
     +       4.8227427715572812d+04,       2.1931899927300499d+06,
     +      -3.2498241761511480d+06,      -2.6728280674740388d+06,
     +      -4.5888327486573579d+05,       1.1694393569568871d+07,
     +      -1.4922492558146100d+07,       1.1873983135276791d+07,
     +       9.6450736351764053d+06,      -3.1069014656318650d+07,
     +       2.7591809451809291d+07,      -2.0115395392681140d+07,
     +       4.4113595732568488d+07,      -2.4883681851902229d+08,
     +       8.1007746052769613d+08,       6.6703395896429726d+01,
     +      -2.7774431192813559d+01,      -1.0974660358925711d+02,
     +       2.1099216151844911d+02,       5.6715188437828402d+01,
     +      -1.9717343061719109d+03,       2.6914556196246631d+03,
     +       2.6686329705526200d+03,       1.8298944683400509d+03,
     +      -1.2570449185949679d+04,       1.2744019476506381d+04,
     +      -8.4247528932371551d+03,      -9.7256877203747536d+03,
     +       2.6870080620669589d+04,      -2.1912013615433680d+04,
     +       1.5092641849717251d+04,      -3.7864354301907202d+04,
     +       2.2892405201343281d+05,      -7.6072382417437423d+05,
     +      -3.2388229676191299d-03,      -1.4115337744204740d-02,
     +       4.2829822924725168d-02,      -5.5858309648411757d-02,
     +      -3.3700361863942201d-02,       5.8649226817820488d-01,
     +      -7.5040933294340206d-01,      -8.0006961050701730d-01,
     +      -8.3223461272456578d-01,       3.9360675404012579d+00,
     +      -3.6029062362008561d+00,       2.2976523737159060d+00,
     +       3.0434611593720429d+00,      -7.9839121143944229d+00,
     +       6.2298341053705508d+00,      -4.1370783091998362d+00,
     +       1.0952322056018090d+01,      -6.8337438010848260d+01,
     +       2.2980916122231170d+02,      -1.5826830261808611d+12,
     +       2.1201702636239480d+12,      -1.3019928623900769d+12,
     +       5.4721940255339752d+11,       2.1217823222160059d+11,
     +       1.3264583589482109d+12,      -1.3500264336283270d+10,
     +       3.9611358597071033d+11,      -7.9922924192387537d+11,
     +      -1.3186358165379680d+13,       4.8573576048441562d+12,
     +       2.1371540037361781d+13,      -8.4950587781846260d+12,
     +      -1.2536148681360211d+13,       1.2810612567529119d+13,
     +      -9.9084671503957930d+12,       9.5374558729616953d+12,
     +      -1.8803442052937250d+13,       2.6663159125413148d+13,
     +       6.4644802699131069d+09,      -8.9088995850494614d+09,
     +       6.0256064651140280d+09,      -8.8599587962379122d+08,
     +      -3.8619136775967669d+09,      -5.1240967221140027d+08,
     +      -9.6340668249925995d+09,      -7.6351333662793608d+09,
     +       1.0049257912446501d+10,       9.8894452451940186d+10,
     +      -4.2926899294411034d+10,      -1.4508052246687790d+11,
     +       6.5736406622471107d+10,       1.9300147079245981d+11,
     +      -3.6520948563107343d+11,       3.4252961313087567d+11,
     +      -2.2801643595053430d+11,       2.2596386462675491d+11,
     +      -5.6601088671370618d+11,      -5.5851360247220118d+06,
     +       8.8261854355868734d+06,      -7.1305985577467382d+06,
     +      -2.1865792103102710d+06,       1.2139717132202620d+07,
     +      -1.3519374738476360d+07,       4.1412932962860152d+07,
     +       1.7842943924440801d+07,      -3.8169488089813590d+07,
     +      -2.1590280245743999d+08,       1.4069967553462100d+08,
     +       2.3650125253444129d+08,      -1.2005850639248520d+08,
     +      -4.9231062319324088d+08,       1.0275002948835230d+09,
     +      -9.8831667840298879d+08,       6.4952505209808660d+08,
     +      -6.6266199979057157d+08,       1.9098021469625061d+09,
     +       2.2020544462118619d+03,      -4.3707097766881116d+03,
     +       4.5883610661558296d+03,       3.3247234288835739d+03,
     +      -1.3294214105384841d+04,       1.9529881415726919d+04,
     +      -4.6175465913604152d+04,      -2.6365015604665110d+04,
     +       5.5727085819099259d+04,       2.0670840692269360d+05,
     +      -1.6513977709140669d+05,      -1.7799132256328309d+05,
     +       9.2502429277857227d+04,       4.8645828160846542d+05,
     +      -1.0378537402057750d+06,       1.0070246279363590d+06,
     +      -6.7575258586353308d+05,       7.8998156212901440d+05,
     +      -2.3576311232477319d+06,      -3.7473586837410350d-01,
     +       9.4774347567775252d-01,      -1.2161141251654910d+00,
     +      -1.1560801903287881d+00,       4.4110862627591443d+00,
     +      -7.0718950513361820d+00,       1.5244777982101979d+01,
     +       1.0371389164809109d+01,      -2.0522752991975651d+01,
     +      -6.5554585675066505d+01,       5.6501583410685022d+01,
     +       4.6918377865564061d+01,      -2.2810802297326369d+01,
     +      -1.5548985746974671d+02,       3.3232849369139149d+02,
     +      -3.2354673824402801d+02,       2.2133152403636890d+02,
     +      -2.8471466853806601d+02,       8.5658418006346812d+02,
     +      -2.6934134432281210d+11,       3.6939910085454559d+11,
     +      -2.2586595605861609d+11,       1.5886556891849171d+10,
     +       2.1779734303222641d+11,       3.9691821998306671d+11,
     +      -6.1689060693424939d+11,       2.4484379271773901d+11,
     +      -1.3570292934096819d+12,      -2.0199190696220359d+12,
     +       1.3529887754469470d+12,       1.1538181061284471d+13,
     +      -4.7171993097619775d+12,      -1.1305723279137580d+13,
     +       8.7794857498825781d+12,       4.0638954534787708d+11,
     +      -4.3198283842088438d+13,       2.2949357498994241d+14,
     +      -4.3912957412996131d+14,       8.2382964216975296d+08,
     +      -1.0500785935965420d+09,       6.0949160777599478d+08,
     +      -2.0512853051897921d+07,      -1.1355296321491330d+09,
     +      -1.3300266958118651d+09,       5.2913521101643858d+09,
     +       1.4125443470254691d+09,       4.6036220908458538d+09,
     +      -1.3211721172588200d+10,       1.1426595359361700d+10,
     +      -5.8839551042940903d+10,       2.1468924045911179d+10,
     +       3.0484421210983459d+10,       4.9654407363520760d+10,
     +      -1.2650471697562740d+11,       4.5841671743241638d+11,
     +      -2.1403854485423491d+12,       4.6575009747726377d+12,
     +      -7.2202517678205634d+05,       7.4930318330212624d+05,
     +      -3.9281822765392618d+05,       4.2448844444057561d+05,
     +       1.4684268268799670d+06,       1.7208241070075110d+06,
     +      -1.1609337008726809d+07,      -5.9128916818380523d+06,
     +      -1.3370496976036981d+07,       5.3372607814468458d+07,
     +      -3.2974158922387101d+07,       1.2809114402316231d+08,
     +      -5.7288695766876549d+07,      -9.5601390474132318d+06,
     +      -2.5582417204565141d+08,       4.4425370781002480d+08,
     +      -1.1917351424165349d+09,       5.2208624724414616d+09,
     +      -1.1843367037185480d+10,       2.3038329526521071d+02,
     +      -8.5693805174160488d+01,       2.9524459090665530d+01,
     +      -6.5351660759045342d+02,      -7.3643908680577579d+02,
     +      -1.3321611041509509d+03,       1.0290505916496049d+04,
     +       4.8204443074742203d+03,       1.6749513418767139d+04,
     +      -5.2231293281069004d+04,       2.6029392569632651d+04,
     +      -1.2280318533591120d+05,       6.6152661996155352d+04,
     +      -2.0122182474648911d+04,       2.9312017249948502d+05,
     +      -4.6497526582422241d+05,       1.1110868058679050d+06,
     +      -4.7223895367325880d+06,       1.1053126212146370d+07,
     +      -1.9408932116736641d-02,      -4.2851138954434848d-02,
     +       2.4187361572976750d-02,       2.4329184230204809d-01,
     +       1.3316291220411869d-01,       4.0570350407656908d-01,
     +      -3.0864765057849342d+00,      -1.1075749427118480d+00,
     +      -6.0768993346728779d+00,       1.5371770530766360d+01,
     +      -6.5121511620340549d+00,       3.9174240289267829d+01,
     +      -2.3079524178638771d+01,       1.0419559939547590d+01,
     +      -9.7178091985139147d+01,       1.4867138145454510d+02,
     +      -3.3558872969641891d+02,       1.4042614272278470d+03,
     +      -3.3611097575200829d+03,      -3.3545069507264961d+12,
     +       5.0733097518704473d+12,      -3.4420013895603950d+12,
     +       5.2751013858656238d+11,       2.0415152326773120d+12,
     +       5.3954084646337754d+12,       3.5970585536322432d+12,
     +      -3.3945769205343550d+12,      -5.2015861530005947d+12,
     +      -1.1593714939590000d+14,       9.8165406361723562d+13,
     +       1.3521590462912691d+14,      -8.7919395890969984d+13,
     +      -1.9775303203690531d+14,       3.8942446562329081d+14,
     +      -3.7779281449528812d+14,       1.3566377845112630d+14,
     +       1.0120873277955990d+15,      -3.4912814354010510d+15,
     +       2.9923197142638199d+10,      -4.6752832478448509d+10,
     +       3.4774906328373749d+10,       1.3638575675974550d+09,
     +      -3.9231089177458397d+10,      -2.8258060268106621d+10,
     +      -4.5442842623940208d+10,       1.2706307845330111d+10,
     +      -2.7671654542190289d+10,       1.6462789736748831d+12,
     +      -1.5188363454822690d+12,      -1.3704702226605510d+12,
     +       8.8246580957945264d+11,       2.5572688323552930d+12,
     +      -4.7248858722714883d+12,       4.2837853574306938d+12,
     +       5.5294019846608899d+11,      -2.8790056827967051d+13,
     +       1.0193553914157880d+14,      -5.1278116004847758d+07,
     +       8.5492688558475703d+07,      -6.7546573788967997d+07,
     +      -1.6586750609953590d+07,       1.1878726020165810d+08,
     +       7.6652701165772434d+06,       1.3056897958915401d+08,
     +      -1.2450480025793131d+08,       3.8341050638738060d+08,
     +      -4.2211000818770232d+09,       4.2842596386325660d+09,
     +       2.3482929437904649d+09,      -1.4933544330853910d+09,
     +      -6.7749789092882566d+09,       1.1836988389719780d+10,
     +      -1.0146668972661690d+10,      -5.8313544302033501d+09,
     +       1.0613546411962891d+11,      -3.7048170640816718d+11,
     +       3.3800892571010387d+04,      -6.0415652418209022d+04,
     +       5.0008871775706248d+04,       2.3162051833505149d+04,
     +      -1.1994390090302000d+05,       1.9412929547701260d+04,
     +      -1.0734062772881020d+05,       1.6226227653221361d+05,
     +      -5.4813123120440240d+05,       3.9034735705099939d+06,
     +      -4.1047607534503210d+06,      -1.6480543911343799d+06,
     +       9.9611755748741084d+05,       6.3588848184758620d+06,
     +      -1.0654496544502540d+07,       8.7565162393176015d+06,
     +       8.3169976194061320d+06,      -1.2015666959125760d+08,
     +       4.1743225550466019d+08,      -8.0536381331256646d+00,
     +       1.5245911290528079d+01,      -1.3026614119162311d+01,
     +      -8.5694382863311969d+00,       3.8426465919927743d+01,
     +      -8.7084466178171169d+00,       2.4692368683034552d+01,
     +      -5.2022046243977357d+01,       1.9859062243300610d+02,
     +      -1.1823439242955940d+03,       1.2526919786189669d+03,
     +       4.2923989954501820d+02,      -2.3425457111799551d+02,
     +      -1.9564621185556859d+03,       3.1895348195621500d+03,
     +      -2.5507170624623582d+03,      -3.0783207497036428d+03,
     +       4.0780549546547947d+04,      -1.4136032562552599d+05,
     +       3.4524979534739868d+11,      -4.8690563703054761d+11,
     +       2.7869215727998163d+11,       2.2776083165088100d+10,
     +       1.9839733221043430d+11,      -1.6681071354004500d+12,
     +       1.2043601839755181d+12,       1.2674964658024661d+12,
     +      -2.0015831879581289d+12,       2.9094200660306230d+12,
     +       4.4805004168354297d+12,      -1.5305919678759949d+13,
     +       6.7175547196606201d+12,       2.0991894810671051d+13,
     +      -4.6144538191837906d+13,       4.4230652711457703d+13,
     +      -1.6527491259041070d+13,      -1.2692958163508150d+13,
     +       8.6115261953580625d+13,      -7.1730597709884572d+08,
     +       9.6477562984659982d+08,      -3.0172008927912432d+08,
     +      -3.8432101906794930d+08,      -2.5030995046409330d+09,
     +       1.0636593620766859d+10,      -7.9571636079044323d+09,
     +      -1.4222149656599421d+10,       1.2116661585664610d+10,
     +       3.2485513893530739d+10,      -4.7560522034661400d+10,
     +       4.7620201574395050d+10,      -1.5993804577718460d+10,
     +      -1.6549086756665741d+11,       3.5802903015452173d+11,
     +      -3.3366119116681018d+11,       1.4754798826640091d+11,
     +      -1.8919516495250720d+11,      -4.4963935406334180d+11,
     +       4.5405567885257659d+05,      -4.7884061427922262d+05,
     +      -4.2782611449300102d+05,       7.2834043238672765d+05,
     +       7.1611568088154774d+06,      -2.4673139414040491d+07,
     +       1.7799398131114230d+07,       3.6410004681018867d+07,
     +      -1.6039477925067211d+07,      -1.2208221660892101d+08,
     +       1.1826034212565599d+08,      -8.3986249347872928d+07,
     +       3.6014883701695167d+07,       3.9669067794195998d+08,
     +      -8.5844945246303570d+08,       7.9563443102476394d+08,
     +      -4.2269567650454009d+08,       1.0139679617438580d+09,
     +       7.1266125342777038d+08,      -9.6492442437319781d+01,
     +      -3.1001033194509770d+01,       6.4816888635189662d+02,
     +      -4.9098540792182018d+02,      -6.9757828204682874d+03,
     +       2.2352688500444950d+04,      -1.5322225918745500d+04,
     +      -3.3561859132150217d+04,       5.1080288488160786d+03,
     +       1.2738823961334850d+05,      -1.0491624832710640d+05,
     +       7.2785867520769621d+04,      -3.8920580146931898d+04,
     +      -3.6857222938163037d+05,       7.9599064511469798d+05,
     +      -7.3795856526985357d+05,       4.4361443435831828d+05,
     +      -1.2655148157111940d+06,      -5.5710150604309759d+05,
     +      -1.2270172311996520d-03,       5.4601786932553173d-02,
     +      -2.1494729793072809d-01,       1.1149648854830320d-01,
     +       2.1603764276166531d+00,      -6.6706378165191724d+00,
     +       4.3396158676879759d+00,       1.0186354233277680d+01,
     +       2.3260572620352679d-01,      -4.0686929807407743d+01,
     +       3.0544679804555521d+01,      -2.1441788311835939d+01,
     +       1.2893104619913220d+01,       1.1253658401214680d+02,
     +      -2.4261549092187869d+02,       2.2516288967405811d+02,
     +      -1.4502462453992490d+02,       4.3841195397321508d+02,
     +       1.8763680653261110d+02,       4.4328725531657285d+12,
     +      -6.2041596406818438d+12,       3.5805244474914570d+12,
     +       6.7264756584358093d+11,      -1.8926365267184219d+11,
     +      -1.7059242309090609d+13,       1.0234780983164350d+13,
     +       1.5749772472482980d+13,      -4.1615711606977681d+12,
     +       4.8338754883841094d+13,      -1.6136415307392920d+12,
     +      -1.4082631457892400d+14,       7.4382784031930438d+13,
     +      -1.3398204211689080d+14,       3.7600512046378638d+14,
     +      -3.8404591587545881d+14,       1.4087907994155500d+14,
     +       7.8620537060618375d+14,      -3.0104341602744740d+15,
     +      -5.6783102297427778d+09,       8.8184659700661831d+09,
     +      -4.1856397890374360d+09,      -7.0909679792541695d+09,
     +       4.8743784363895388d+09,       3.5011238920721123d+10,
     +       1.7063153067668770d+10,      -8.6592819235673996d+10,
     +      -1.1176696158172810d+11,      -1.2760231485252940d+11,
     +       6.6118542851148889d+11,      -6.4542980375893457d+11,
     +       8.1363759514666150d+11,      -7.2819455389378552d+11,
     +      -1.5024800988054800d+11,       2.5552533391764050d+11,
     +       2.5835884738934268d+12,      -1.7254509367795801d+13,
     +       5.5657043820618617d+13,       3.0802294429393872d+06,
     +      -6.0282994100024393d+06,       2.1777619860571930d+06,
     +       1.5109656674149269d+07,      -1.2804579951779710d+07,
     +      -4.9084604917353272d+07,      -4.3645859421440430d+07,
     +       1.6535236278353199d+08,       4.1561627041117609d+08,
     +      -2.9295954231387001d+07,      -2.0330463824380281d+09,
     +       2.4452065409470258d+09,      -2.1872416262828951d+09,
     +       3.1271435017651801d+09,      -3.0888308441510558d+09,
     +       3.4470190810422468d+09,      -1.3513186274369101d+10,
     +       6.7079530503230911d+10,      -1.9184366965134460d+11,
     +      -8.3246455450077008d+02,       2.4500285728755211d+03,
     +      -8.1221079336653895d+02,      -1.1601116043765540d+04,
     +       1.1122958465004820d+04,       3.5546969799149600d+04,
     +       3.2227820638990670d+04,      -1.3794515077876789d+05,
     +      -4.4851085359734861d+05,       2.0975863252118949d+05,
     +       2.0184981375501191d+06,      -2.5541819272406478d+06,
     +       2.0033406630924621d+06,      -3.3906852032481511d+06,
     +       4.2923624042196209d+06,      -5.0246236702282894d+06,
     +       1.6727904844845859d+07,      -7.9228279247815266d+07,
     +       2.1625170539089060d+08,       9.0372807783604944d-02,
     +      -4.6557603762187733d-01,       1.9790237020757501d-01,
     +       3.0427671885722321d+00,      -3.2793504373907081d+00,
     +      -9.4377930421697034d+00,      -7.7565224509125592d+00,
     +       3.9905391919182883d+01,       1.4578519869198999d+02,
     +      -9.7879524889807328d+01,      -6.2531388818918322d+02,
     +       8.0679164930758964d+02,      -5.8268295756989448d+02,
     +       1.0825390555866211d+03,      -1.5084240863200059d+03,
     +       1.8145111404909480d+03,      -5.8729999780213348d+03,
     +       2.7543082748948400d+04,      -7.3496365568028457d+04,
     +      -6.3780258207819238d+12,       1.0629574101059320d+13,
     +      -6.8742277049766709d+12,      -7.1263455931075745d+11,
     +      -6.7282156766233760d+12,       3.8028643084074000d+13,
     +       1.4674095171901289d+13,       2.2766462130461941d+13,
     +       3.7469303405537047d+13,      -1.0195051728647190d+15,
     +       7.3680135741659125d+14,       6.3367424901400400d+14,
     +       5.7910844101571885d+12,      -2.0114148221247340d+15,
     +       3.6267728296854220d+15,      -4.1118929513197060d+15,
     +       6.8166726241642150d+15,      -2.5755288593228200d+16,
     +       7.3293742119344640d+16,       7.7595104515441864d+10,
     +      -1.3091787655063390d+11,       8.5531982371018005d+10,
     +       2.0823569530312222d+10,       2.6745142540144798d+10,
     +      -4.4004317747892969d+11,       5.9386301914046692d+10,
     +      -6.5506466242190186d+11,      -6.2544657836278796d+11,
     +       1.5319171168597039d+13,      -1.3661608685126920d+13,
     +      -4.5992019864384297d+12,      -4.5311854926140723d+12,
     +       3.1679782584693301d+13,      -4.9562799963349047d+13,
     +       5.1310130119224547d+13,      -5.0227049757048367d+13,
     +      -3.6768484800619789d+13,       6.5285634706438512d+14,
     +      -1.7078338911978641d+08,       2.9551193120247543d+08,
     +      -1.8414802096233511d+08,      -1.0806090914982840d+08,
     +       1.1082180599105459d+08,       8.4625282941051579d+08,
     +      -6.7038783599497712d+08,       1.7659492063337760d+09,
     +       4.0202176950989780d+09,      -4.0445512666169357d+10,
     +       3.7879495737031464d+10,       1.9671126055043681d+09,
     +       2.3707763699342850d+10,      -9.1849444252476257d+10,
     +       1.3243625242786349d+11,      -1.2984876353776340d+11,
     +       5.4964640195016167d+10,       7.9326183649967456d+11,
     +      -4.2506799549148511d+12,       1.3396238271296190d+05,
     +      -2.3722049389534959d+05,       1.3994555293615311d+05,
     +       1.2805670051879520d+05,      -2.1896736963295500d+05,
     +      -5.2482154312962468d+05,       9.0352653637291747d+05,
     +      -1.4676530939639539d+06,      -5.5771696343429117d+06,
     +       3.6919894402109303d+07,      -3.5172994990021780d+07,
     +       4.0109747087367959d+06,      -2.9158275104090970d+07,
     +       8.9061541159051314d+07,      -1.2031081145577730d+08,
     +       1.1226630274518740d+08,       1.1518257189721599d+07,
     +      -1.2561116065013199d+09,       5.7849645456666317d+09,
     +      -3.5419426368639613d+01,       6.3857528419220237d+01,
     +      -3.5970001958530148d+01,      -4.2984310329567030d+01,
     +       8.8206019101189753d+01,       1.0645050030713800d+02,
     +      -3.2912908237386961d+02,       3.9328340736956562d+02,
     +       2.0365464797036991d+03,      -1.0890577849315730d+04,
     +       1.0467574289062461d+04,      -2.1749642349664100d+03,
     +       9.9703286108223710d+03,      -2.7366337947043528d+04,
     +       3.5458964978682649d+04,      -3.1815895041589931d+04,
     +      -1.7418175916909549d+04,       4.9013131349031022d+05,
     +      -2.1341468619698151d+06,       7.9426602929339551d+11,
     +      -1.1280639505169971d+12,       6.3271585521045276d+11,
     +       2.0128898382322491d+11,      -7.3673932510955109d+10,
     +      -4.0793029058826880d+12,       4.0069080706939160d+12,
     +      -3.3804409452571738d+12,       8.3616481135499658d+12,
     +       3.0932763031220141d+13,      -4.3486538519636748d+12,
     +      -8.2891767717054781d+13,       8.6143151682943574d+12,
     +       1.2000798687561809d+14,      -1.5754881767665031d+14,
     +       1.8137815867484581d+14,      -8.5057407244967100d+14,
     +       9.3922621818947180d+15,      -4.5606678651700704d+16,
     +      -2.5870730293330731d+09,       3.3505397550331130d+09,
     +      -1.5764839285542860d+09,      -6.1153054088288033d+08,
     +      -1.2496690497214571d+08,       2.0217989265959850d+10,
     +      -3.7481606280916039d+10,       4.0076919070231762d+09,
     +      -3.7385175372714752d+10,       4.5169338709134239d+10,
     +      -1.7072482643472519d+11,       3.8526693135201282d+11,
     +       4.3332348893252548d+10,      -5.1371127208921661d+11,
     +       3.9517373859198578d+11,      -4.6607619097283209d+11,
     +       6.5044781553152617d+12,      -9.1095831533377281d+13,
     +       4.4756457766244612d+14,       2.4092484398511210d+06,
     +      -2.4902647527301749d+06,       8.0937044425879687d+05,
     +      -1.5986167626132050d+06,       4.6809067529723551d+06,
     +      -3.8861907690289170d+07,       8.0243156390665472d+07,
     +       1.7119642267602101d+07,       1.1661101374909291d+08,
     +      -3.4316856091986811d+08,       4.9500123990822142d+08,
     +      -7.2321127548131776d+08,      -1.2117055017178421d+08,
     +       6.9945060907691658d+08,       1.3089662702408610d+08,
     +      -2.8173660158540729d+07,      -1.6504274283963091d+10,
     +       2.5362142901644080d+11,      -1.2448506961064141d+12,
     +      -8.2637148500855949d+02,       3.2231169510623090d+02,
     +       1.1915485654288990d+02,       3.1440589510224690d+03,
     +      -6.9851965137350517d+03,       3.4543727658498690d+04,
     +      -6.8996017885392183d+04,      -1.7765801455346362d+04,
     +      -1.4078509423504479d+05,       3.6404404005275777d+05,
     +      -4.4670340098560508d+05,       6.6418888819548022d+05,
     +       3.2040213902658601d+04,      -3.9713539803470368d+05,
     +      -5.5316529864871199d+05,       4.6692114731614583d+05,
     +       1.5685141301353451d+07,      -2.5060112976992011d+08,
     +       1.2279531841193521d+09,       8.0815697847355789d-02,
     +       1.3386532150500910d-01,      -1.1021594142898331d-01,
     +      -1.2217263881548699d+00,       2.5556352808175928d+00,
     +      -1.0618923989662200d+01,       2.0395111330295531d+01,
     +       3.9725854927015090d+00,       4.9557720298510503d+01,
     +      -1.0938160049781720d+02,       1.2861183548796751d+02,
     +      -2.1114293327570641d+02,       1.0664992459233540d+01,
     +       8.4944750234753329d+01,       2.2284556718320681d+02,
     +      -1.9720187815387851d+02,      -4.8249282441093101d+03,
     +       7.8607705343874797d+04,      -3.8476550644650048d+05,
     +       7.7510222432228027d+12,      -1.2228756860159609d+13,
     +       8.2107065354963584d+12,      -8.9572157897159082d+11,
     +       2.1266553257001970d+12,      -2.8775451634283770d+13,
     +      -4.8623717643234297d+13,       7.2277376288547922d+13,
     +       1.0054455463782400d+13,       7.5789331717297188d+14,
     +      -8.1525824619286875d+14,      -5.1126940810433312d+14,
     +       1.5541783547757169d+14,       1.4324854161620860d+15,
     +      -2.2042585088816120d+15,       2.2678835790273380d+15,
     +      -2.9210667145174050d+15,       5.6179491894326830d+15,
     +      -2.2456230947793280d+16,      -6.4241001590602386d+10,
     +       1.0649351254201981d+11,      -8.2155320003668198d+10,
     +      -9.0798895413062401d+09,       8.7199498569213516d+10,
     +      -1.1372703148891939d+10,       7.4371531588773596d+11,
     +      -5.1623484920186261d+11,       6.6115278167468359d+11,
     +      -1.0512736257265070d+13,       9.1687328777927461d+12,
     +       1.0925219584655330d+13,      -9.1519871195245449d+12,
     +      -1.8747384923412539d+13,       4.2839476595952023d+13,
     +      -4.4738159085504883d+13,       2.2749004012247230d+13,
     +       2.1761688086508969d+13,       6.4311774377471225d+14,
     +       1.0598067276570460d+08,      -1.9009760280708781d+08,
     +       1.6458558791542369d+08,       4.8718405133479953d+07,
     +      -3.7595704159191030d+08,       6.3980720825959337d+08,
     +      -2.4142597601932092d+09,       1.4631952979620309d+09,
     +      -3.2536437742548032d+09,       2.7503244478616459d+10,
     +      -2.3353252996001221d+10,      -2.8793748933494900d+10,
     +       2.6950758854875851d+10,       5.1439186895990479d+10,
     +      -1.2421072589964880d+11,       1.2677566880596201d+11,
     +      -1.9716643281859089d+10,      -4.2031059560605042d+11,
     +      -9.1603637461210791d+11,      -6.9402231085420164d+04,
     +       1.3557633297380671d+05,      -1.2972685549600430d+05,
     +      -6.4129956812742479d+04,       4.2738671901074558d+05,
     +      -8.3589315963841637d+05,       2.3628267484830339d+06,
     +      -1.3282416097409411d+06,       3.8887239592322651d+06,
     +      -2.6098245271464702d+07,       2.1141970535328548d+07,
     +       2.8940364577347919d+07,      -2.8240345014862221d+07,
     +      -4.9271838809473768d+07,       1.2300079144104950d+08,
     +      -1.2322253339227280d+08,      -1.4197796798638981d+07,
     +       6.8774196846717465d+08,       6.2314453173045628d+07,
     +       1.6683730448070659d+01,      -3.4973215812445957d+01,
     +       3.5936876723274906d+01,       2.3602276898691940d+01,
     +      -1.4299494795048329d+02,       2.7801541581832618d+02,
     +      -6.9975223746723191d+02,       3.5656971273214720d+02,
     +      -1.3232192154349880d+03,       8.0925001371592016d+03,
     +      -6.2289601126117268d+03,      -9.4367867806163267d+03,
     +       9.2100483049300728d+03,       1.5263125013627559d+04,
     +      -3.8670449841666981d+04,       3.8206069532307833d+04,
     +       1.1332866578711570d+04,      -2.7166622094454680d+05,
     +       1.6445012153941451d+05,      -6.5769181031733301d+11,
     +       1.0056900576527220d+12,      -5.8838078884567041d+11,
     +      -3.6938202298016559d+11,      -1.3680758957870110d+12,
     +       8.6428459982791328d+12,      -3.4223487989422231d+12,
     +      -1.3235787664634689d+13,       8.6932133491306650d+12,
     +      -1.4488201956796080d+13,      -2.0651834799730830d+12,
     +       5.9120606764350320d+13,      -9.2896260559008828d+12,
     +      -1.5524339663059259d+14,       2.7972456341607159d+14,
     +      -2.9310630251022019d+14,       4.7805304241521012d+14,
     +      -4.0156890837754040d+15,       1.9631857068374300d+16,
     +       8.8326839632423937d+08,      -1.5368649420164540d+09,
     +       1.0859478907713521d+09,       4.6400456038449490d+08,
     +       9.6541887258665829d+09,      -3.7837062784846817d+10,
     +      -1.5267876554911120d+09,       1.2719077383849370d+11,
     +      -4.3665370164903191d+10,      -1.5525315041048499d+11,
     +       1.1359856375506300d+11,      -3.7149863654576459d+11,
     +       2.2181817455863629d+11,       8.3265950867487061d+11,
     +      -1.5316373602625549d+12,       1.6758121396694109d+12,
     +      -3.9701061850867759d+12,       2.7847549392458102d+13,
     +      -8.6920258233061031d+13,      -4.6887785587881628d+05,
     +       1.1553310233443240d+06,      -1.6592661959015939d+06,
     +       1.8427853321530630d+06,      -2.5679062933691829d+07,
     +       8.3475950539978683d+07,       1.5561997438899729d+07,
     +      -3.4818362352975482d+08,       4.6388403041038550d+07,
     +       6.6030040133642125d+08,      -3.1119349175159413d+08,
     +       8.4620463144327843d+08,      -7.7052959483631349d+08,
     +      -1.8922929359471941d+09,       3.8098344320803242d+09,
     +      -4.2842781834990692d+09,       9.8353695993362312d+09,
     +      -5.4659789097338242d+10,       9.5456550653366562d+10,
     +       9.5569729300171744d+01,      -5.2357102847364376d+02,
     +       1.5152929132625809d+03,      -3.2096334147602661d+03,
     +       2.5002207260060641d+04,      -7.4893120399619933d+04,
     +      -2.1435089285013561d+04,       3.4928383481557848d+05,
     +      -1.2552311888223410d+04,      -7.4912579740462429d+05,
     +       2.8045339344443422d+05,      -7.0576122017405275d+05,
     +       7.6864442953934183d+05,       1.7404836023802231d+06,
     +      -3.6189023199916729d+06,       4.1279864834072180d+06,
     +      -9.2261744116004799d+06,       4.2512160397654101d+07,
     +      -2.0759779791963831d+07,       1.1295501841159789d-04,
     +       1.0941951462161660d-01,      -4.8469182108063258d-01,
     +       1.2269797252606149d+00,      -7.7653471794090656d+00,
     +       2.2326943982424769d+01,       8.0377442294608041d+00,
     +      -1.1152604749477690d+02,      -1.3986155966931539d+00,
     +       2.5147595203181120d+02,      -8.2103877937529020d+01,
     +       1.9574654276833670d+02,      -2.3723282848124319d+02,
     +      -5.2858163114891738d+02,       1.1135234683312840d+03,
     +      -1.2838095728587959d+03,       2.8488892919861901d+03,
     +      -1.1542954096369949d+04,      -5.7844910357441149d+03/

      END

